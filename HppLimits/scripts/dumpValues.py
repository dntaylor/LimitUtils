#!/usr/bin/env python
import os
import sys
import math
import json
import pickle
import errno
import ROOT

# helper functions
def python_mkdir(dir):
    '''A function to make a unix directory as well as subdirectories'''
    try:
        os.makedirs(dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(dir):
            pass
        else: raise


def dumpResults(results,name):
    jfile = '{0}.json'.format(name)
    pfile = '{0}.pkl'.format(name)
    if os.path.dirname(jfile): python_mkdir(os.path.dirname(jfile))
    if os.path.dirname(pfile): python_mkdir(os.path.dirname(pfile))
    with open(jfile,'w') as f:
        f.write(json.dumps(results, indent=4, sort_keys=True))
    with open(pfile,'wb') as f:
        pickle.dump(results,f)


def printObjects(workspace,func):
    print func
    args = getattr(workspace,func)()
    if isinstance(args,ROOT.RooArgSet):
        args.Print()
    else:
        for arg in args:
            print arg
def dumpDict(d):
    print json.dumps(d,sort_keys=True,indent=4)

def printDict(d):
    for k,i in sorted(d.iteritems()):
        print '    ',k,i.getVal()

def getArgsetMap(workspace,func):
    allVars = {}
    args = getattr(workspace,func)()
    it = args.createIterator()
    var = it.Next()
    while var:
        allVars[var.GetName()] = var
        var = it.Next()
    return allVars

def getVals(allFuncs,doSB=False, channels=[]):
    apMap = {}
    ppMap = {}
    bgMap = {}
    hpps = {
        'll' : ['ee','em','mm'],
        'el' : ['ee','em'],
        'ml' : ['em','mm'],
        'ee' : ['ee'],
        'em' : ['em'],
        'mm' : ['mm'],
        'lt' : ['et','mt'],
        'et' : ['et'],
        'mt' : ['mt'],
        'tt' : ['tt']
    }
    hms = {
        'l' : ['e','m'],
        'e' : ['e'],
        'm' : ['m'],
        't' : ['t'],
    }
    allChannels = []
    for chan in channels:
        if len(chan) == 3:
            hpp = chan[:2]
            hm = chan[2:]
            for i in hpps[hpp]:
                for j in hms[hm]:
                    allChannels += [i+j]
        if len(chan) == 4:
            hpp = chan[:2]
            hmm = chan[2:]
            for i in hpps[hpp]:
                for j in hpps[hmm]:
                    allChannels += [i+j]
    for f,v in allFuncs.iteritems():
        if channels and not any(['_{0}_'.format(c) in f for c in allChannels]): continue
        if doSB and 'SB' not in f: continue
        if not doSB and 'SB' in f: continue
        if 'datadriven' in f:
            bgMap[f] = v.getVal()
        elif 'HppHmm' in f:
            ppMap[f] = v.getVal()
        else:
            apMap[f] = v.getVal()
    return apMap, ppMap, bgMap

def getUncertainty(expMap,shiftMap):
    total = 0.
    totalShift = 0.
    for f in expMap:
        total += expMap[f]
        totalShift += shiftMap[f]
    unc = abs(total-totalShift)/total if total else 0.
    return unc

def varyNuisances(allVars, allFuncs, doSB=False,channels=[]):
    # get unvaried expected for each channel
    apMap, ppMap, bgMap = getVals(allFuncs,doSB=doSB,channels=channels)
    ap = sum([v for f,v in apMap.iteritems()])
    pp = sum([v for f,v in ppMap.iteritems()])
    bg = sum([v for f,v in bgMap.iteritems()])
    # vary each nuisance independently and get change in expected
    apVaryMap = {'up':{}, 'down':{}}
    ppVaryMap = {'up':{}, 'down':{}}
    bgVaryMap = {'up':{}, 'down':{}}
    for f,v in allVars.iteritems():
        # remove stuff that isnt an uncertainty
        if '_In' in f: continue
        if 'CMS_fake' in f: continue
        if f in ['r','MH']: continue

        if 'alpha_13TeV80X' in f:
            # vary gmN
            start = v.getVal()
            uperr = start*(1 + 1/math.sqrt(start+1))
            downerr = start*(1 - 1/math.sqrt(start+1))
            v.setVal(uperr)
            apVaryMap['up'][f], ppVaryMap['up'][f], bgVaryMap['up'][f] = getVals(allFuncs,doSB=doSB,channels=channels)
            v.setVal(downerr)
            apVaryMap['down'][f], ppVaryMap['down'][f], bgVaryMap['down'][f] = getVals(allFuncs,doSB=doSB,channels=channels)
            v.setVal(start)
        else:
            # vary lnN
            v.setVal(1.)
            apVaryMap['up'][f], ppVaryMap['up'][f], bgVaryMap['up'][f] = getVals(allFuncs,doSB=doSB,channels=channels)
            v.setVal(-1.)
            apVaryMap['down'][f], ppVaryMap['down'][f], bgVaryMap['down'][f] = getVals(allFuncs,doSB=doSB,channels=channels)
            v.setVal(0.)
    # sum of squares the changes
    apErr2 = {'up':0.,'down':0.}
    ppErr2 = {'up':0.,'down':0.}
    bgErr2 = {'up':0.,'down':0.}
    for f in allVars:
        if f not in bgVaryMap['up']: continue
        apErr2['up'] += getUncertainty(apMap,apVaryMap['up'][f])**2
        ppErr2['up'] += getUncertainty(ppMap,ppVaryMap['up'][f])**2
        bgErr2['up'] += getUncertainty(bgMap,bgVaryMap['up'][f])**2
        apErr2['down'] += getUncertainty(apMap,apVaryMap['down'][f])**2
        ppErr2['down'] += getUncertainty(ppMap,ppVaryMap['down'][f])**2
        bgErr2['down'] += getUncertainty(bgMap,bgVaryMap['down'][f])**2
    return ap, pp, bg, ap*apErr2['up']**0.5, pp*ppErr2['up']**0.5, bg*bgErr2['up']**0.5, ap*apErr2['down']**0.5, pp*ppErr2['down']**0.5, bg*bgErr2['down']**0.5

def getCardValues(analysis,mode,mass,channels=[]):
    filename = 'working/{0}/{1}/higgsCombineTest.Asymptotic.mH{2}.root'.format(analysis,mode,mass)
    tfile = ROOT.TFile(filename)
    
    ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
    
    workspace = tfile.Get("w")
    
    allVars = getArgsetMap(workspace,'allVars')
    allPdfs = getArgsetMap(workspace,'allPdfs')
    allFuncs = getArgsetMap(workspace,'allFunctions')
    
    #print 'vars'
    #printDict(allVars)
    #print 'pdfs'
    #printDict(allPdfs)
    #print 'functions'
    #printDict(allFuncs)
    #print 'vars', len(allVars), 'pdfs', len(allPdfs), 'functions', len(allFuncs)

    apValSR, ppValSR, bgValSR, apErrUpSR, ppErrUpSR, bgErrUpSR, apErrDownSR, ppErrDownSR, bgErrDownSR = varyNuisances(allVars, allFuncs, doSB=False,channels=channels)
    apValSB, ppValSB, bgValSB, apErrUpSB, ppErrUpSB, bgErrUpSB, apErrDownSB, ppErrDownSB, bgErrDownSB = varyNuisances(allVars, allFuncs, doSB=True,channels=channels)

    return {
        'apSR': {'val': apValSR, 'errUp': apErrUpSR, 'errDown': apErrDownSR,},
        'ppSR': {'val': ppValSR, 'errUp': ppErrUpSR, 'errDown': ppErrDownSR,},
        'bgSR': {'val': bgValSR, 'errUp': bgErrUpSR, 'errDown': bgErrDownSR,},
        'apSB': {'val': apValSB, 'errUp': apErrUpSB, 'errDown': apErrDownSB,},
        'ppSB': {'val': ppValSB, 'errUp': ppErrUpSB, 'errDown': ppErrDownSB,},
        'bgSB': {'val': bgValSB, 'errUp': bgErrUpSB, 'errDown': bgErrDownSB,},
    }
    

analyses = ['Hpp3lAP','Hpp3lPP','Hpp4l','HppAP','HppPP','HppComb']
analyses = ['HppComb']
modes = ['ee100','em100','et100','mm100','mt100','tt100','BP1','BP2','BP3','BP4']
modes = ['ee100','em100','et100','mm100','mt100','tt100']
masses = [200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]

channels = {
    'ee100' : {'lll': ['eel'], 'llt': ['eet'], 'llll': ['eeee']},
    'em100' : {'lll': ['eml'], 'llt': ['emt'], 'llll': ['emem']},
    'mm100' : {'lll': ['mml'], 'llt': ['mmt'], 'llll': ['mmmm']},
    'et100' : {'lll': ['ell'], 'llt': ['elt'], 'ltl': ['etl'], 'ltt': ['ett'], 'llll': ['elel'], 'lllt': ['elet'], 'ltlt': ['etet']},
    'mt100' : {'lll': ['mll'], 'llt': ['mlt'], 'ltl': ['mtl'], 'ltt': ['mtt'], 'llll': ['mlml'], 'lllt': ['mlmt'], 'ltlt': ['mtmt']},
    'tt100' : {'lll': ['lll'], 'llt': ['llt'], 'ltl': ['ltl'], 'ltt': ['ltt'], 'ttl': ['ttl'], 'ttt': ['ttt'],  'llll': ['llll'], 'lllt': ['lllt'], 'lltt': ['lltt'], 'ltlt': ['ltlt'], 'lttt': ['lttt'], 'tttt': ['tttt']},
}

analyses = ['HppComb']
#modes = ['ee100']
#masses = [200]
data = {}
for analysis in analyses:
    data[analysis] = {}
    for mode in modes:
        data[analysis][mode] = {}
        for mass in masses:
            allVals = {}
            for chan in sorted(channels[mode]):
                print analysis,mode,mass,chan
                allVals[chan] = getCardValues(analysis,mode,mass,channels=channels[mode][chan])
            data[analysis][mode][mass] = allVals
            dumpResults(data,'limit_uncertainties')
