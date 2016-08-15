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

def getLimits(analysis,mode,mass,outDir,prod='',doImpacts=False,retrieve=False,submit=False,dryrun=False,jobName='',skipAsymptotic=False):
    '''
    Submit a job using farmoutAnalysisJobs --fwklite
    '''
    srcdir = os.path.join(os.environ['CMSSW_BASE'],'src')
    datacard = 'datacards/{0}/{1}/{2}{3}.txt'.format(analysis,mode,mass,prod)
    workspace = 'datacards/{0}/{1}/{2}{3}.root'.format(analysis,mode,mass,prod)
    impacts = 'impacts/{0}/{1}/{2}{3}.json'.format(analysis,mode,mass,prod)
    outimpacts = 'impacts/{0}/{1}/{2}{3}'.format(analysis,mode,mass,prod)

    # mkdirs
    python_mkdir('{2}/datacards/{0}/{1}'.format(analysis,mode,srcdir))
    python_mkdir('{2}/impacts/{0}/{1}'.format(analysis,mode,srcdir))

    # combine cards
    if analysis=='HppAP':
        # just cp
        runCommand('pushd {0}; cp datacards/Hpp3l/{1}/{2}AP.txt {3}'.format(srcdir,mode,mass,datacard))
    if analysis=='HppPP':
        # combine 3l PP and 4l PP
        runCommand('pushd {0}; combineCards.py datacards/Hpp3l/{1}/{2}PP.txt datacards/Hpp4l/{1}/{2}.txt > {3}'.format(srcdir,mode,mass,datacard))
    if analysis=='HppComb':
        # combein 3l AP, 3l PP, and 4l PP
        runCommand('pushd {0}; combineCards.py datacards/Hpp3l/{1}/{2}.txt datacards/Hpp4l/{1}/{2}.txt > {3}'.format(srcdir,mode,mass,datacard))

    # get datacard path relative to $CMSSW_BASE
    dfull = os.path.abspath(os.path.join(os.environ['CMSSW_BASE'],'src',datacard))
    ddir = os.path.dirname(dfull)
    dname = os.path.basename(dfull)
    cmssw_base = os.environ['CMSSW_BASE']
    dsplit = dfull.split('/')
    srcpos = dsplit.index('src')
    drel = '/'.join(dsplit[srcpos:])
    dreldir = '/'.join(dsplit[srcpos:-1])

    # first, get the approximate bounds from asymptotic
    work = 'working/{0}{2}/{1}'.format(analysis,mode,prod)
    workfull = os.path.join(srcdir,work)
    python_mkdir(workfull)
    combineCommand = 'combine -M Asymptotic {0} -m {1} --saveWorkspace'.format(dfull,mass)
    command = 'pushd {0}; nice {1};'.format(workfull,combineCommand) 
    logging.info('{0}:{1}:{2}: Finding Asymptotic limit: {3}'.format(analysis,mode,mass,datacard))
    logging.debug('{0}:{1}:{2}: {3}'.format(analysis,mode,mass,combineCommand))


    name = 'asymptotic'
    fileDir = '{4}/{0}/{1}/{2}/{3}'.format(name,analysis,mode,mass,srcdir)
    python_mkdir(fileDir)
    fileName = '{0}/limits{1}.txt'.format(fileDir,prod)
    if skipAsymptotic and os.path.isfile(fileName):
        with open(fileName,'r') as f:
            quartiles = [float(x) for x in f.readlines()[0].split()]
    else:
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
            outline = ' '.join([str(x) for x in quartiles])
            logging.info('{0}:{1}:{2}: Limits: {3}'.format(analysis,mode,mass,outline))

        with open(fileName,'w') as f:
            outline = ' '.join([str(x) for x in quartiles])
            f.write(outline)

    if submit:
        sample_dir = '/nfs_scratch/{0}/{1}/{2}/{3}/{4}{5}'.format(pwd.getpwuid(os.getuid())[0], jobName, analysis, mode, mass, prod)

        # create submit dir
        submit_dir = '{0}/submit'.format(sample_dir)
        if os.path.exists(submit_dir):
            logging.warning('Submission directory exists for {0}.'.format(jobName))
            return
        # setup the job parameters
        rmin = 0.8*min(quartiles[:5])
        rmax = 1.2*max(quartiles[:5])
        num_points = 100
        points_per_job = 5
        toys = 10000

        # create dag dir
        dag_dir = '{0}/dags/dag'.format(sample_dir)
        os.system('mkdir -p {0}'.format(os.path.dirname(dag_dir)))
        os.system('mkdir -p {0}'.format(dag_dir+'inputs'))

        # output dir
        output_dir = 'srm://cmssrm.hep.wisc.edu:8443/srm/v2/server?SFN=/hdfs/store/user/{0}/{1}/{2}/{3}/{4}{5}'.format(pwd.getpwuid(os.getuid())[0], jobName, analysis, mode, mass, prod)

        # create file list
        rlist = [r*(rmax-rmin)/num_points + rmin for r in range(int(num_points/points_per_job))]
        input_name = '{0}/rvalues.txt'.format(dag_dir+'inputs')
        with open(input_name,'w') as file:
            for r in rlist:
                file.write('{0}\n'.format(r))

        # create bash script
        bash_name = '{0}/{1}.sh'.format(dag_dir+'inputs', jobName)
        bashScript = '#!/bin/bash\n'
        bashScript += 'printenv\n'
        bashScript += 'read -r RVAL < $INPUT\n'
        for i in range(points_per_job):
            dr = i*(rmax-rmin)/points_per_job
            bashScript += 'combine $CMSSW_BASE/{0} -M HybridNew --freq -s -1 --singlePoint $(bc -l <<< "$RVAL+{1}") --saveToys --fullBToys --clsAcc 0 --saveHybridResult -m {2} -n Tag -T {3} -i 2\n'.format(drel,dr,mass,toys)
        bashScript += 'hadd $OUTPUT higgsCombineTag.HybridNew.mH{0}.*.root\n'.format(mass)
        bashScript += 'rm higgsCombineTag.HybridNew.mH{0}.*.root\n'.format(mass)
        with open(bash_name,'w') as file:
            file.write(bashScript)
        os.system('chmod +x {0}'.format(bash_name))

        # create farmout command
        farmoutString = 'farmoutAnalysisJobs --infer-cmssw-path --fwklite --input-file-list={0} --assume-input-files-exist'.format(input_name)
        farmoutString += ' --submit-dir={0} --output-dag-file={1} --output-dir={2}'.format(submit_dir, dag_dir, output_dir)
        farmoutString += ' --extra-usercode-files="{0}" {1} {2}'.format(dreldir, jobName, bash_name)

        if not dryrun:
            logging.info('Submitting {0}/{1}/{2}/{3}{4}'.format(jobName,analysis,mode,mass,prod))
            os.system(farmoutString)
        else:
            print farmoutString


    # now do the higgs combineharvester stuff
    if doImpacts:
        wfull = os.path.abspath(os.path.join(os.environ['CMSSW_BASE'],'src',workspace))
        ifull = os.path.abspath(os.path.join(os.environ['CMSSW_BASE'],'src',impacts))
        logging.info('{0}:{1}:{2}: text2workspace'.format(analysis,mode,mass))
        command = 'text2workspace.py {0} -m {1}'.format(dfull,mass)
        runCommand(command)
        logging.info('{0}:{1}:{2}: Impacts: initial fit'.format(analysis,mode,mass))
        command = 'pushd {0}; nice combineTool.py -M Impacts -d {1} -m {2} --doInitialFit --robustFit 1'.format(workfull,wfull,mass)
        runCommand(command)
        logging.info('{0}:{1}:{2}: Impacts: nuissance fits'.format(analysis,mode,mass))
        command = 'pushd {0}; nice combineTool.py -M Impacts -d {1} -m {2} --robustFit 1 --doFits'.format(workfull,wfull,mass)
        runCommand(command)
        logging.info('{0}:{1}:{2}: Impacts: saving/plotting'.format(analysis,mode,mass))
        command = 'pushd {0}; combineTool.py -M Impacts -d {1} -m {2} -o {3}'.format(workfull,wfull,mass,ifull)
        runCommand(command)
        command = 'pushd {0}; plotImpacts.py -i {1} -o {2}'.format(srcdir,ifull,outimpacts)
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
    parser.add_argument('-j','--jobName', nargs='?',type=str,default='',help='Jobname for submission')
    parser.add_argument('-s','--submit',action='store_true',help='Submit Full CLs')
    parser.add_argument('-sa','--skipAsymptotic',action='store_true',help='Skip calculating asymptotic (read from file)')
    parser.add_argument('-r','--retrieve',action='store_true',help='Retrieve Full CLs')
    parser.add_argument('-dr','--dryrun',action='store_true',help='Dryrun for submission')
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
            if len(allowedMasses)==1 or args.submit:
                for m in allowedMasses:
                    if an=='Hpp3l':
                        getLimits(an,bp,m,args.directory,'AP',args.impacts,args.retrieve,args.submit,args.dryrun,args.jobName,args.skipAsymptotic)
                        getLimits(an,bp,m,args.directory,'PP',args.impacts,args.retrieve,args.submit,args.dryrun,args.jobName,args.skipAsymptotic)
                    else:
                        getLimits(an,bp,m,args.directory,'',args.impacts,args.retrieve,args.submit,args.dryrun,args.jobName,args.skipAsymptotic)
            else:
                allArgs = []
                for m in allowedMasses:
                    if an=='Hpp3l':
                        newArgs = [an,bp,m,args.directory,'AP',args.impacts,args.retrieve,args.submit,args.dryrun,args.jobName,args.skipAsymptotic]
                        allArgs += [newArgs]
                        newArgs = [an,bp,m,args.directory,'PP',args.impacts,args.retrieve,args.submit,args.dryrun,args.jobName,args.skipAsymptotic]
                        allArgs += [newArgs]
                    else:
                        newArgs = [an,bp,m,args.directory,'',args.impacts,args.retrieve,args.submit,args.dryrun,args.jobName,args.skipAsymptotic]
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

