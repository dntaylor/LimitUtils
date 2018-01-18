#!/usr/bin/env python

import ROOT
import math
import json

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

def getVals(allFuncs):
    expMap = {}
    expMapB = {}
    for f,v in allFuncs.iteritems():
        #if 'proc' in f: continue
        if 'SB' in f: continue
        if 'bonly' in f:
            expMapB[f] = v.getVal()
        else:
            expMap[f] = v.getVal()
    return expMap, expMapB

def getUncertainty(expMap,shiftMap):
    total = 0.
    totalShift = 0.
    for f in expMap:
        total += expMap[f]
        totalShift += shiftMap[f]
    unc = abs(total-totalShift)/total if total else 0.
    return unc

def varyNuisances(allVars, allFuncs, *nuis):
    # get unvaried expected for each channel
    expMap, expMapB = getVals(allFuncs)
    # vary each nuisance independently and get change in expected
    varyMap = {'up':{}, 'down':{}}
    varyMapB = {'up':{}, 'down':{}}
    for n in nuis:
        if n not in allVars:
            print 'Unrecognized nuisance {0}'.format(n)
            continue
        allVars[n].setVal(1.)
        varyMap['up'][n], varyMapB['up'][n] = getVals(allFuncs)
        allVars[n].setVal(-1.)
        varyMap['down'][n], varyMapB['down'][n] = getVals(allFuncs)
        allVars[n].setVal(0.)
    # sum of squares the changes
    err2 = {'up':0.,'down':0.}
    for n in nuis:
        if n not in allVars: continue
        err2['up'] += getUncertainty(expMap,varyMap['up'][n])**2
        err2['down'] += getUncertainty(expMap,varyMap['down'][n])**2
    return (err2['up']**0.5 + err2['down']**0.5)/2.

def getCardUncertainties(analysis,mode,mass):
    filename = 'working/{0}/{1}/higgsCombineTest.Asymptotic.mH{2}.root'.format(analysis,mode,mass)
    tfile = ROOT.TFile(filename)
    
    ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
    
    workspace = tfile.Get("w")
    
    allVars = getArgsetMap(workspace,'allVars')
    allPdfs = getArgsetMap(workspace,'allPdfs')
    allFunctions = getArgsetMap(workspace,'allFunctions')
    
    #print 'vars'
    #printDict(allVars)
    #print 'pdfs'
    #printDict(allPdfs)
    #print 'functions'
    #printDict(allFunctions)
    
    uncertainties = {}
    uncertainties['lumi'] =      varyNuisances(allVars,allFunctions,*[x for x in allVars if x.startswith('lumi') and not x.endswith('In')])
    uncertainties['sigAP'] =       varyNuisances(allVars,allFunctions,*['sig_unc_AP'])
    uncertainties['sigPP'] =       varyNuisances(allVars,allFunctions,*['sig_unc_PP'])
    #uncertainties['charge'] =    varyNuisances(allVars,allFunctions,*[x for x in allVars if 'charge' in x and not x.endswith('In')])
    uncertainties['elec_id'] =   varyNuisances(allVars,allFunctions,*['elec_id'])
    uncertainties['muon_id'] =   varyNuisances(allVars,allFunctions,*['muon_id'])
    uncertainties['tau_id'] =    varyNuisances(allVars,allFunctions,*['tau_id'])
    uncertainties['stat'] =      varyNuisances(allVars,allFunctions,*[x for x in allVars if x.startswith('stat') and not x.endswith('In')])
    uncertainties['alpha_unc'] = varyNuisances(allVars,allFunctions,*[x for x in allVars if x.startswith('alpha_unc') and not x.endswith('In')])
    return uncertainties

analyses = ['Hpp3lAP','Hpp3lPP','Hpp3lPPR','Hpp4l','Hpp4lR','HppAP','HppPP','HppPPR','HppComb']
modes = ['ee100','em100','et100','mm100','mt100','tt100','BP1','BP2','BP3','BP4']
masses = [200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]
modes = ['mm100']
#masses = [200]
unc = {}
for mode in modes:
    unc[mode] = {}
    for mass in masses:
        unc[mode][mass] = getCardUncertainties('HppComb',mode,mass)

    keys = sorted(unc[mode][200].keys())

    print ' '.join(['{0:10}'.format(k) for k in [mode]+keys])
    for mass in masses:
        print ' '.join(['{0:10}'.format(x) for x in [mass]+['{0:10.4}'.format(unc[mode][mass][k]*100.) for k in keys]])
    print ''
