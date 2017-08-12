#!/usr/bin/env python
# coding: utf8


import argparse
from Bio import motifs

def plot_pwm(filename,output=None):
    """
    Plots a mononuc probability matrix from a JASPAR file
    Caution: only ONE motif in the file !!!
    The output must be given as the name of a PDF file, otherwise automatic from the input name
    """
    if output is None:
        output=filename.split("/")[-1].split(".")[0]+".pdf"
    # Uses WebLogo from BioPython to plot the matrix in pdf format
    # Requires an internet connection
    fh = open(filename)
    for m in motifs.parse(fh, "jaspar"):
        try:
            m.weblogo(output,format="PDF")
        except:
            print "ERROR trying to plot the sequence logo using the online website WebLogo (http://weblogo.berkeley.edu/). A possible cause is the absence of a working internet connection, otherwise the pwm file %s may be corrupted. If you wish to plot the sequence logo, consider using the alternate website STAMP (http://www.benoslab.pitt.edu/stamp/)"%filename
    fh.close()
    return 0



# -------------- EXECUTION
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots a WebLogo sequence motif from a position weight-matrix file in JASPAR format, using the online wrapper of BioPython. Output in pdf format.') 
    parser.add_argument('input', type=str,help='Input JASPAR file')
    parser.add_argument("-o","--output",type=str,action="store",help="Output pdf file")
    args=parser.parse_args()
    plot_pwm(unicode(args.input),unicode(args.output)) #executes program
