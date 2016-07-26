#!/usr/bin/env python
'''
A script to retrieve the limits

Author: Devin N. Taylor, UW-Madison
'''

import os
import sys
import glob
import pwd
import argparse
import errno
import socket
import signal
import logging
import math
import ROOT
import subprocess
from multiprocessing import Pool

masses = [200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500]

def python_mkdir(dir):
    '''A function to make a unix directory as well as subdirectories'''
    try:
        os.makedirs(dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(dir):
            pass
        else: raise

def limitsWrapper(args):
    getLimits(*args)

def runCommand(command):
    return subprocess.Popen(command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0]

def getLimits(analysis,mode,mass,outDir,prod='',doImpacts=False):
    '''
    Submit a job using farmoutAnalysisJobs --fwklite
    '''
    datacard = 'datacards/{0}/{1}/{2}{3}.txt'.format(analysis,mode,mass,prod)
    workspace = 'datacards/{0}/{1}/{2}{3}.root'.format(analysis,mode,mass,prod)
    impacts = 'impacts/{0}/{1}/{2}{3}.json'.format(analysis,mode,mass,prod)
    outimpacts = 'impacts/{0}/{1}/{2}{3}'.format(analysis,mode,mass,prod)

    # mkdirs
    python_mkdir('datacards/{0}/{1}'.format(analysis,mode))
    python_mkdir('impacts/{0}/{1}'.format(analysis,mode))

    # combine cards
    if analysis=='HppAP':
        # just cp
        runCommand('cp datacards/Hpp3l/{1}/{2}AP.txt {0}'.format(datacard,mode,mass))
    if analysis=='HppPP':
        # combine 3l PP and 4l PP
        runCommand('combineCards.py datacards/Hpp3l/{1}/{2}PP.txt datacards/Hpp4l/{1}/{2}.txt > {0}'.format(datacard,mode,mass))
    if analysis=='HppComb':
        # combein 3l AP, 3l PP, and 4l PP
        runCommand('combineCards.py datacards/Hpp3l/{1}/{2}.txt datacards/Hpp4l/{1}/{2}.txt > {0}'.format(datacard,mode,mass))

    # get datacard path relative to $CMSSW_BASE
    dfull = os.path.abspath(os.path.join(os.path.dirname(__file__),datacard))
    ddir = os.path.dirname(dfull)
    dname = os.path.basename(dfull)
    cmssw_base = os.environ['CMSSW_BASE']
    dsplit = dfull.split('/')
    srcpos = dsplit.index('src')
    drel = '/'.join(dsplit[srcpos:])
    dreldir = '/'.join(dsplit[srcpos:-1])

    # first, get the approximate bounds from asymptotic
    work = 'working/{0}{2}/{1}'.format(analysis,mode,prod)
    python_mkdir(work)
    workfull = os.path.join(os.path.dirname(__file__),work)
    combineCommand = 'combine -M Asymptotic {0} -m {1} --saveWorkspace'.format(dfull,mass)
    command = 'pushd {0}; nice {1};'.format(workfull,combineCommand) 
    logging.info('{0}:{1}:{2}: Finding Asymptotic limit: {3}'.format(analysis,mode,mass,datacard))
    logging.debug('{0}:{1}:{2}: {3}'.format(analysis,mode,mass,combineCommand))
    out = runCommand(command)

    fname = os.path.join(workfull, "higgsCombineTest.Asymptotic.mH{0}.root".format(mass))
    file = ROOT.TFile(fname,"READ")
    tree = file.Get("limit")
    if not tree: 
        logging.warning('{0}:{1}:{2}: Asymptotic presearch failed'.format(analysis,mode,mass))
        quartiles = [0., 0., 0., 0., 0., 0.]
    else:
        quartiles = []
        for i, row in enumerate(tree):
            quartiles += [row.limit]


    # now do the higgs combineharvester stuff
    if doImpacts:
        wfull = os.path.abspath(os.path.join(os.path.dirname(__file__),workspace))
        ifull = os.path.abspath(os.path.join(os.path.dirname(__file__),impacts))
        logging.info('{0}:{1}:{2}: text2workspace'.format(analysis,mode,mass))
        command = 'text2workspace.py {0} -m {1}'.format(datacard,mass)
        runCommand(command)
        logging.info('{0}:{1}:{2}: Impacts: initial fit'.format(analysis,mode,mass))
        command = 'pushd {0}; nice combineTool.py -M Impacts -d {1} -m {2} --doInitialFit --robustFit 1'.format(work,wfull,mass)
        runCommand(command)
        logging.info('{0}:{1}:{2}: Impacts: nuissance fits'.format(analysis,mode,mass))
        command = 'pushd {0}; nice combineTool.py -M Impacts -d {1} -m {2} --robustFit 1 --doFits'.format(work,wfull,mass)
        runCommand(command)
        logging.info('{0}:{1}:{2}: Impacts: saving/plotting'.format(analysis,mode,mass))
        command = 'pushd {0}; combineTool.py -M Impacts -d {1} -m {2} -o {3}'.format(work,wfull,mass,ifull)
        runCommand(command)
        command = 'plotImpacts.py -i {0} -o {1}'.format(ifull,outimpacts)
        runCommand(command)

    # now get the fullCLs
    args = [
        ['Expected 0.025', '--expectedFromGrid 0.025', 'higgsCombineTest.HybridNew.mH{0}.quant0.025.root'],
        ['Expected 0.160', '--expectedFromGrid 0.160', 'higgsCombineTest.HybridNew.mH{0}.quant0.160.root'],
        ['Expected 0.500', '--expectedFromGrid 0.500', 'higgsCombineTest.HybridNew.mH{0}.quant0.500.root'],
        ['Expected 0.840', '--expectedFromGrid 0.840', 'higgsCombineTest.HybridNew.mH{0}.quant0.840.root'],
        ['Expected 0.975', '--expectedFromGrid 0.975', 'higgsCombineTest.HybridNew.mH{0}.quant0.975.root'],
        ['Observed',       '',                         'higgsCombineTest.HybridNew.mH{0}.root'],
    ]
    
    ## merge the output
    #gridfile = 'grid_{0}.root'.format(mass)
    #sourceDir = '/hdfs/store/user/dntaylor/2016-01-21_allLimits_10KToys_100Points_v1/{0}/{1}/{2}'.format(analysis,mode,mass)
    #logging.info('{0}:{1}:{2}: Merging: {3}'.format(analysis,mode,mass,sourceDir))
    #haddCommand = 'hadd {0} {1}/*.root'.format(gridfile,sourceDir)
    #command = 'pushd {0}; nice {1};'.format(workfull,haddCommand)
    #logging.debug('{0}:{1}:{2}: {3}'.format(analysis,mode,mass,haddCommand))
    #out = runCommand(command)

    ## get CL
    #fullQuartiles = []
    #for i in range(len(args)):
    #    logging.info('{0}:{1}:{2}: Calculating: {3}'.format(analysis,mode,mass,args[i][0]))
    #    combineCommand = 'combine {0} -M HybridNew --freq --grid={1} -m {2} --rAbsAcc 0.001 --rRelAcc 0.001 {3}'.format(dfull, gridfile, mass, args[i][1])
    #    command = 'pushd {0}; nice {1};'.format(workfull, combineCommand)
    #    logging.debug('{0}:{1}:{2}: {3}'.format(analysis,mode,mass,combineCommand))
    #    outfile = '{0}/{1}'.format(workfull,args[i][2].format(mass))
    #    out = runCommand(command)

    #    # read the limit
    #    file = ROOT.TFile(outfile,"READ")
    #    tree = file.Get("limit")
    #    if not tree:
    #        logging.warning('HybridNew failed')
    #        val = 0.
    #    else:
    #        val = 0.
    #        for i, row in enumerate(tree):
    #            val = row.limit
    #    fullQuartiles += [val]

    # save the values
    quartileMap = {
        'asymptotic' : quartiles,
        #'fullCLs' : fullQuartiles,
    }
    for name in ['asymptotic']:
        fileDir = '{0}/{1}/{2}/{3}'.format(name,analysis,mode,mass)
        if outDir: fileDir = outDir + '/' + fileDir
        python_mkdir(fileDir)
        fileName = '{0}/limits{1}.txt'.format(fileDir,prod)
        with open(fileName,'w') as f:
            outline = ' '.join([str(x) for x in quartileMap[name]])
            logging.info('{0}:{1}:{2}: Limits: {3} - {4}'.format(analysis,mode,mass,name, outline))
            f.write(outline)



def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Process limits')

    parser.add_argument('-a','--analysis', nargs='?',type=str,const='',choices=['Hpp4l','Hpp3l','HppAP','HppPP','HppComb'],help='Analysis to process')
    parser.add_argument('-d','--directory', nargs='?',type=str,const='',help='Custom out directory')
    parser.add_argument('-m','--mass', nargs='?',type=str,const='500',help='Mass for Higgs combine')
    parser.add_argument('-bp','--branchingPoint',nargs='?',type=str,const='BP4',default='BP4',choices=['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4'],help='Choose branching point')
    parser.add_argument('-ab','--allBranchingPoints',action='store_true',help='Run over all branching points')
    parser.add_argument('-am','--allMasses',action='store_true',help='Run over all masses')
    parser.add_argument('-aa','--allAnalyses',action='store_true',help='Run over all anlayses')
    parser.add_argument('-i','--impacts',action='store_true',help='Do impacts (slower)')
    parser.add_argument('-l','--log',nargs='?',type=str,const='INFO',default='INFO',choices=['INFO','DEBUG','WARNING','ERROR','CRITICAL'],help='Log level for logger')

    args = parser.parse_args(argv)

    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')

    allowedAnalyses = ['Hpp4l','Hpp3l','HppAP','HppPP','HppComb'] if args.allAnalyses else [args.analysis]
    allowedBranchingPoints = ['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4'] if args.allBranchingPoints else [args.branchingPoint]
    allowedMasses = masses if args.allMasses else [args.mass]

    for an in allowedAnalyses:
        for bp in allowedBranchingPoints:
            if len(allowedMasses)==1:
                if an=='Hpp3l':
                    getLimits(an,bp,allowedMasses[0],args.directory,'AP',args.impacts)
                    getLimits(an,bp,allowedMasses[0],args.directory,'PP',args.impacts)
                else:
                    getLimits(an,bp,allowedMasses[0],args.directory,'',args.impacts)
            else:
                allArgs = []
                for m in allowedMasses:
                    if an=='Hpp3l':
                        newArgs = [an,bp,m,args.directory,'AP',args.impacts]
                        allArgs += [newArgs]
                        newArgs = [an,bp,m,args.directory,'PP',args.impacts]
                        allArgs += [newArgs]
                    else:
                        newArgs = [an,bp,m,args.directory,'',args.impacts]
                        allArgs += [newArgs]
                p = Pool(14)
                try:
                    p.map_async(limitsWrapper, allArgs).get(999999)
                except KeyboardInterrupt:
                    p.terminate()
                    print 'limits cancelled'
                    sys.exit(1)

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)

