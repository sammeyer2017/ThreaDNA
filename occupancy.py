#!/usr/bin/python
# coding: utf8
import sys,os
import numpy as np
import argparse

def occ(en,efft):
    return np.exp(-(en-np.mean(en))/efft)

def main(args):
    if args.efftemp==None:
        args.efftemp=.3
    else:
        args.efftemp=float(args.efftemp)
    if args.output==None:
        outputf=args.inf.split(".")[0]+"_occ.bed"    
    else:
        outputf=args.output

    
    # --- compute occupancy profile
    en=np.loadtxt(args.inf,skiprows=1,usecols=(1,2,3))
    occprof=occ(en[:,2],args.efftemp)
    tab=np.vstack((en[:,(0,1)].T,[occprof])).T
    
    # --- read name of chromosome
    a=open(args.inf,"r")
    b=a.readline()
    b=a.readline()
    se=b.split()[0]
    a.close()

    # --- write output file
    if float(np.__version__[:3])>=1.7:
        np.savetxt(outputf, tab, fmt=se+'\t%.1f\t%.1f\t%.3f', header='track type=bedGraph name=\"relative occupancy (a. u.)\"',comments="")
    else:
        np.savetxt(outputf, tab, fmt=se+'\t%.1f\t%.1f\t%.3f')#,header='track type=bedGraph name=\"raw binding free energy (k_B T)\"')
        os.system(r"sed -i -e '1itrack type=bedGraph name=\"relative occupancy (a. u.)\"\' %s"%outputf)


    # ---- optional: coverage
    if args.coverage!=None:
	print args.coverage
        covprof=np.convolve(occprof,np.ones(int(args.coverage))/int(args.coverage),"same")
        tab=np.vstack((en[:,(0,1)].T,[covprof])).T
        # cov output file: 
        if args.output==None:
            covoutput=args.inf.split(".")[0]+"_cov.bed"
        else:
            covoutput=args.output.split(".")[0]+"_cov.bed"
        # export
        if float(np.__version__[:3])>=1.7:
            np.savetxt(covoutput, tab, fmt=se+'\t%.1f\t%.1f\t%.3f', header='track type=bedGraph name=\"relative coverage (a. u.)\"',comments="")
        else:
            np.savetxt(covoutput, tab, fmt=se+'\t%.1f\t%.1f\t%.3f')#,header='track type=bedGraph name=\"raw binding free energy (k_B T)\"')
            os.system(r"sed -i -e '1itrack type=bedGraph name=\"relative coverage (a. u.)\"\' %s"%covoutput)
    # ----------------------

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='ThreaDNA helper program, that transforms an energy profile (output from the main program) into a protein occupancy profile along a given DNA sequence. Should preferably be applied along a single sequence.') 
    parser.add_argument('inf', metavar='input_file', type=str,help='Input energy profile in bedgraph format, typically an output from ThreaDNA.')
    parser.add_argument("-out","--output",type=str,help="Optional filename for the output occupancy profile. Default: same as input with _occ suffix.")
    parser.add_argument("-t","--efftemp",type=float,help="Optional value for the effective temperature (reduced units): scales how variations in deformation energy impact the occupancy profile. Default .3")
    parser.add_argument("-cov","--coverage",type=float,help="Optional: compute the coverage profile in addition to occupancy, i.e. the probability that a given position is covered by a protein (rather than the protein center). You must then provide the protein size in basepairs (e.g. 147 for a nucleosome). Output file has the same name as the main output, with the '_cov' suffix")
    args=parser.parse_args()
    main(args)
