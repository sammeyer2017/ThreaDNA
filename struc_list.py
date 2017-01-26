#!/usr/bin/python
# coding: utf8
import argparse
import os
import sys

basedir=os.path.abspath(os.path.dirname(sys.argv[0])) #installation folder absolute path

def update_list(): #updates the current protein structres list
    os.chdir(basedir+"/structures/")
    prot,data,NP,ABCi,ABCs,ABColdi,ABColds=[],[],[],[],[],[],[]
    for p,d,f in os.walk(os.getcwd()): # explore struct directories
        if not not d :
            prot.append(d)
    for i,n in enumerate(prot[0]): # select protein name
        sub=[]
        for s in prot[i+1]: #select structure
            f=open(n+"/"+s+"/"+s+".info")
            res=["NA"] #default
            for line in f:
                if line[:11]=="RESOLUTION.": #extract resolution
                    res=float(line[14:18])
                if line[:4]=="SIZE": #extract size
                    size=int(line[7:-1])
                if line[:3]=="REF": #extract reference position
                    ref=int(line[6:-1])
                if line[:2]=="NP": #extract distribution for NP parameter
                    NPt=eval(line[5:-1])
                if line[:5]=="ABC_i": #extract distribution for NP parameter
                    ABCit=eval(line[8:-1])
                if line[:5]=="ABC_s": #extract distribution for NP parameter
                    ABCst=eval(line[8:-1])
                if line[:9]=="ABC_old_i": #extract distribution for NP parameter
                    ABCito=eval(line[12:-1])
                if line[:9]=="ABC_old_s": #extract distribution for NP parameter
                    ABCsto=eval(line[12:-1])
            f.close()
            sub.append([n,0,s,res,size,ref,NPt,ABCit,ABCst,ABCito,ABCsto]) #list by structures for each protein
        sub=sorted(sub,key = lambda struct : struct[3]) #sort structures by minimum resolution for each protein name
        for i,x in enumerate(sub):
            x[1]=i+1 #append ordered indexes to sub
        data+=[x[:6] for x in sub] #append protein data to final data
        NP+=[x[:3]+x[6] for x in sub]
        ABCi+=[x[:3]+x[7] for x in sub]
        ABCs+=[x[:3]+x[8] for x in sub]
        ABColdi+=[x[:3]+x[9] for x in sub]
        ABColds+=[x[:3]+x[10] for x in sub]
    l=open("list.dat","w") # write list in list.dat file (update if existing)
    n2=open("NP_distrib.dat","w")
    n3=open("ABC_i_distrib.dat","w")
    n4=open("ABC_s_distrib.dat","w")
    n5=open("ABC_old_i_distrib.dat","w")
    n6=open("ABC_old_s_distrib.dat","w")
    l.write("Protein\tRank\tStructure\tResolution\tSize\tPosition\n") #header
    n2.write("Protein\tRank\tStructure\tNP mean\tNP std\n") #header
    n3.write("Protein\tRank\tStructure\tABC_i mean\tABC_i std\n") #header
    n4.write("Protein\tRank\tStructure\tABC_s mean\tABC_s std\n") #header
    n5.write("Protein\tRank\tStructure\tABC_old_i mean\tABC_old_i std\n") #header
    n6.write("Protein\tRank\tStructure\tABC_old_s mean\tABC_old_s std\n") #header
    i=0
    while i<len(data):
        l.write("\t".join([str(x) for x in data[i]])+"\n") #write final data
        n2.write("\t".join([str(x) for x in NP[i]])+"\n") #write NP distribution        
        n3.write("\t".join([str(x) for x in ABCi[i]])+"\n") #write ABCi distribution 
        n4.write("\t".join([str(x) for x in ABCs[i]])+"\n") #write ABCs distribution 
        n5.write("\t".join([str(x) for x in ABColdi[i]])+"\n") #write ABColdi distribution 
        n6.write("\t".join([str(x) for x in ABColds[i]])+"\n") #write ABColds distribution 
        i+=1
    l.close()
    n2.close()
    n3.close()
    n4.close()
    n5.close()
    n6.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Updates the list of structures included in the database. This is useful only if you have added new structures manually.') #program parser and help
    args=parser.parse_args()
    update_list()

