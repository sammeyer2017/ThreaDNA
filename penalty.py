#!/usr/bin/python
# coding: utf8
import sys
import argparse

def seqpen(seqn,seq,char,pos,l,f,ref):
    i=0
    s=""
    while i<len(seq)-l+1:
        f.write(seqn+"\t"+str(ref+i)+"\t"+str(int(seq[i+pos-1]==char))+"\n")
        i+=1

def main(args):
    if args.reference==None:
        args.reference=(args.l+1)/2
    f=open(args.fi,"r")
    f2=open(args.fi.split("/")[-1].split(".")[0]+"_"+args.char+str(args.pos)+".bed","w")
    f2.write("Energy penalties for "+args.char+" in position "+str(args.pos)+"/"+str(args.l)+"\n")
    i=0
    li=f.read().split("\n")
    f.close()
    seq,seqn,si=[],[],[]
    while i<len(li):
        if li[i] and li[i][0]=='>':
            seqn.append(li[i][1:].replace(" ","_")) #replace whitespaces by underscores for further data utilisation (plot)
            si.append(i)
        i+=1
    si.append(i)
    i=0
    while i<len(si)-1:
        seq.append("".join(li[si[i]+1:si[i+1]]).upper())
        seqpen(seqn[i],seq[i],args.char,args.pos,args.l,f2,args.reference)
        i+=1
    f2.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='ThreaDNA helper program, to identify the protein positions matching given sequence patterns. Penalty is given the chosen nucleotide and its position on the protein, and returns the positions of the protein where this pattern is matched. ') #program parser and help
    parser.add_argument('fi', metavar='fasta_file', type=str,help='FASTA file containing the sequence')
    parser.add_argument('char',metavar='nucleotide',type=str,help='Nucleotide query (A, C, G, or T)')
    parser.add_argument("pos",metavar="position",type=int, help="Position of the nucleotide in the protein-DNA model")
    parser.add_argument("l",metavar="length",type=int,help="Length of the protein-DNA complex (in nucleotides)")
    parser.add_argument("-r","--reference",type=int,help="Reference nucleotide for this protein-DNA model (default (length+1)/2)")
    args=parser.parse_args()
    main(args)
