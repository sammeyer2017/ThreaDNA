#!/usr/bin/python
# coding: utf8
import argparse
import numpy as np
import os
import sys
import urllib2
import struc_list as sl
from util import *

basedir=os.path.abspath(os.path.dirname(sys.argv[0])).decode('utf8') #installation folder absolute path

"""
def autoanalyze(id,out):     #analyse sequence on webserver    
    url=urllib2.urlopen('http://ndbserver.rutgers.edu/service/ndb/atlas/stfeatures?searchTarget='+id.lower()+'&ftrType=bpmp&type=csv')
    #open url in ndb for query intra parameters
    ip = url.read()
    intratab=[map(float,x.split(',')[3:]) for x in ip.split('\n')[1:-1]] #convert results to float numbers
    if not intratab:
        sys.exit("ID not found in NDB")
    url=urllib2.urlopen('http://ndbserver.rutgers.edu/service/ndb/atlas/stfeatures?searchTarget='+id.lower()+'&ftrType=bpmsp&type=csv')
    #open url in ndb for query step parameters
    sp = url.read()
    intertab=[map(float,x.split(',')[3:9]) for x in sp.split('\n')[1:-1]]  #convert results to float numbers             
    it=np.array(intratab)
    st=np.array(intertab)
    np.savetxt(out+"_i.dat",it,fmt='%.3f') #save results in files in the output directory
    np.savetxt(out+"_s.dat",st,fmt='%.3f')
    return it,st

def info(pdb,out): #pdb file analysis in tmp directory
    f1=open(pdb,"r") #open pdb file
    f2=open(out+".info","w") #open info file
    for line in f1:
        if any(line[:6]==x for x in ["HEADER","TITLE ","EXPDTA","AUTHOR"]):
            f2.write(line) #write informations
        if "RESOLUTION." in line:
            f2.write(line[11:]) #write resolution
    f1.close()
    f2.close()
    return
"""

def outf(infil,prot): #output prefix function
        out=unicode(infil.split("/")[-1]) #get structure prefix
        outd=basedir+u"/structures/"+prot+u"/"+out+u"/" #output absolute path
        os.system("mkdir -p "+outd.encode(encoding="UTF-8")) #create output directories
        #os.system("cp "+infil+".pdb "+outd) #copy pdb file (optional)
        return outd+out

def out_parameters(infil,out): #extract parameters from .out analysis
    # find both intra and step parameters
    if infil[-4:]!=".out":
        sys.exit("problem not a 3DNA .out file")
    try:
        f=open(infil,"r").readlines()
    except:
        sys.exit("problem with input file: not found")
    for i,li in enumerate(f): # reads data
        if "Local base-pair parameters" in li:
            ni=i
        if "Local base-pair step parameters" in li:
            ns=i
    # find nb of bp
    l=ns-ni-6
    if l<1:
        sys.exit("error: less than 1 bp!")
    intradat=f[ni+2:ni+l+2] #remove headers
    interdat=f[ns+2:ns+l+1]
    for x in intradat: 
        if "----" in x:
            sys.exit("Can't operate with missing base-pair values") #exit if missing values
    for x in interdat:
        if "----" in x:
            sys.exit("Can't operate with missing step values") #exit if missing values
    intratab=[map(float,x.split()[2:]) for x in intradat] #convert results to float numbers
    intertab=[map(float,x.split()[2:]) for x in interdat]
    it=np.array(intratab)
    st=np.array(intertab)
    np.savetxt(out+"_i.dat",it,fmt='%.3f') #save results in fils in the output directory
    np.savetxt(out+"_s.dat",st,fmt='%.3f')
    return it,st

def lis_parameters(infil,out): #extract parameters from .lis file (curve+)
    if infil[-4:]!=".lis":
        sys.exit("problem not a curve+ .lis file")
    try:
        f=open(infil,"r").readlines()
    except:
        sys.exit("problem with input file: not existent")
    for i,li in enumerate(f): # reads data    
        if "Intra-BP" in li:
            ni=i
        if "Inter-BP" in li:
            ns=i
    # find nb of bp
    l=ns-ni-7
    if l<1:
        sys.exit("error: less than 1 bp!")
    intradat=f[ni+4:ni+l+4] #remove headers
    interdat=f[ns+2:ns+l+1]
    intratab=[map(float,x.split()[4:10]) for x in intradat] #convert results to float numbers
    intertab=[map(float,x.split()[4:10]) for x in interdat]
    it=np.array(intratab)
    st=np.array(intertab)
    np.savetxt(out+"_i.dat",it,fmt='%.3f') #save results in fils in the output directory
    np.savetxt(out+"_s.dat",st,fmt='%.3f')
    return it,st

def std_params(f,Qi,Qs):
    print "Testing struture profile on a standard sequence"
    s=read_fasta(basedir+"/structures/standard.fa")[0][0]
    for ca in ["NP","ABC_i","ABC_s","ABC_old_i","ABC_old_s"]:
        q0,K,ind,n=load_dnamodel(ca,basedir+"/")
        if n%2:
            q=np.reshape(np.divide(Qi[:,(3,4,5,0,1,2)],adim(ca)),(-1,1,6))
            q2=np.reshape(np.divide(Qi[(len(Qi)-79)/2:(len(Qi)+79)/2,(3,4,5,0,1,2)],adim(ca)),(-1,1,6))              
        else:
            q=np.reshape(np.divide(Qs[:,(3,4,5,0,1,2)],adim(ca)),(-1,1,6))
            q2=np.reshape(np.divide(Qs[(len(Qs)-78)/2:(len(Qs)+78)/2,(3,4,5,0,1,2)],adim(ca)),(-1,1,6)) 
        E=calc_E(q,q0,K)
        E2=calc_E(q2,q0,K)
        P=np.shape(E)[0]
        P2=np.shape(E2)[0]
        Et=np.zeros(len(s)-n-P+2,dtype=np.float16) #initializes energies at 0
        Et2=np.zeros(len(s)-n-P2+2,dtype=np.float16)
        seq_tmp=np.array([ind[s[j:j+n]] for j in range(len(s)-n+1)],dtype=np.uint8) #list of subsequences indices in calculated E array
        j=0
        while j<P: #loop on protein positions
            Et+=E[j,seq_tmp[j:len(s)-n-P+j+2]] #add partial energy for position j of the protein
            j+=1
        j=0
        while j<P2: #loop on protein positions
            Et2+=E2[j,seq_tmp[j:len(s)-n-P2+j+2]]
            j+=1
        #print ca,np.mean(Et2,dtype=np.float128),np.std(Et2,dtype=np.float128)
        f.write(ca+"    "+str([np.mean(Et,dtype=np.float128),np.std(Et,dtype=np.float128)])+"\n")

def getNDBid(id,out): #Get the NDB ID to PDB ID
    print "Requesting 'http://www.rcsb.org/pdb/explore/explore.do?structureId="+id.lower()
    url=urllib2.urlopen('http://www.rcsb.org/pdb/explore/explore.do?structureId='+id.lower()) # PDB url of the protein
    pattern='href="http://ndbserver.rutgers.edu/service/ndb/atlas/summary?searchTarget='
    ip = url.read()
    ndbURL=None
    
    for line in ip.split('\n'):
      if pattern in line:
        ndbURL=line
    
    if not ndbURL: 
      sys.exit("There is not a NDB ID associated to the PDB ID")
    ndbURL=ndbURL.split('=')
    ndb=ndbURL[2:]  
    ndb=(" ".join(ndb))
    ndb=ndb[:-1]
    print "Found the NDB ID associated, that is : "+ndb
    return getCSVfile(id,out)
  
def getCSVfile(id,out):                           
        if id:
    url=urllib2.urlopen('http://ndbserver.rutgers.edu/service/ndb/atlas/stfeatures?searchTarget='+id.lower()+'&ftrType=bpmp&type=csv')
    #open url in ndb for query intra parameters
    ip = url.read()
    print ip
    intratab=[map(float,x.split(',')[3:]) for x in ip.split('\n')[1:-1]] #convert results to float numbers
    if intratab==[]:
        print "ERROR: could not load the intra-bp parameter file. It is probably absent from the NDB database. Try to analyze your pdb file through the Web3DNA webserver, and provide the output .out file to ThreaDNA."
        raise NameError("NDB coordinate file absent")
    url=urllib2.urlopen('http://ndbserver.rutgers.edu/service/ndb/atlas/stfeatures?searchTarget='+id.lower()+'&ftrType=bpmsp&type=csv')
    #open url in ndb for query step parameters
    sp = url.read()
    intertab=[map(float,x.split(',')[3:9]) for x in sp.split('\n')[1:-1]] #convert results to float numbers             
    it=np.array(intratab)
    st=np.array(intertab)
    np.savetxt(out+"_i.dat",it,fmt='%.3f') #save results in files in the output directory
    np.savetxt(out+"_s.dat",st,fmt='%.3f')
    return it,st

def outDirectory(id,name):
    print "Creation of the structure directory "+id
    out=basedir+u"/structures/"+name+u"/"+id+u"/" # Output absolute path
    os.system("mkdir -p "+out.encode(encoding="UTF-8")) #Create output directories
    return out

def getPDBfile(id,name):
    url='https://files.rcsb.org/download/'+id.lower()+'.pdb'
    if not url:
      sys.exit("ID not found in PDB")
    out=outDirectory(id,name)+id
    print "Downloading the PDB file" + url
    PDBfile = urllib2.urlopen(url)
    data = PDBfile.read()
    with open(out+".pdb", "wb") as code:
      code.write(data)
    PDBfile.close()
    return getInfoPDBfile(out,name)
      
def getInfoPDBfile(out,name): 
    if ".pdb" in out:
	directory=os.path.abspath(out)[:-4]
	print "Creation of the structure directory"+name
    else:
	directory=out
    PDBfile=open(directory+".pdb", "r")
    directory=outf(directory,name)
    info=open(directory+".info", "w") #Create info file
    for line in PDBfile:
      if any(line[:6]==x for x in ["HEADER","TITLE ","EXPDTA","AUTHOR"]):
	info.write(line) #write informations
      if "RESOLUTION." in line:
        info.write(line[11:]) #write resolution
    PDBfile.close()
    info.close()
    return directory

def main(name,f,pdb=None,ref=None):
    ext=f.split(".")

    if len(ext)==1: # ! no file info
	out=getPDBfile(f,name)
	Qi,Qs=getNDBid(f,out)

    elif ext[-1]=="pdb":
	out=getInfoPDBfile(f,name)
	Qi,Qs=getNDBid(f[:-4],out)

    elif ext[-1]=="out": # ! no file info
	f=os.path.abspath(f) #converts to absolute path
	out=outf(f[:-4],name) #output absolute path and prefix
	Qi,Qs=out_parameters(f,out)#extraction of parameters from .out file to .dat files
	
    elif ext[-1]=="lis": # ! no file info
	f=os.path.abspath(f) #converts to absolute path
	out=outf(f[:-4],name) #output absolute path and prefix

	Qi,Qs=lis_parameters(f,out)
    else:
	sys.exit("Unknown extension")
    f2=open(out+".info","a+") #reopen info file
    f2.write("SIZE    "+str(Qi.shape[0])+"\n") # append base pairs count if not present
    print "Got parameters of size "+str(Qi.shape[0])+" for this structure"
    if ref==None:
	print "Warning : no reference nucleotide loaded, will take the central nucleotide ("+str((Qi.shape[0]+1)/2)+") as reference"
	f2.write("REF    "+str((Qi.shape[0]+1)/2)+"\n")
    else:
	f2.write("REF    "+str(ref)+"\n")
    std_params(f2,Qi,Qs)
    f2.close()
    sl.update_list() #update struct list
    print("List updated with new structure !")
    return "List updated with new structure !"


    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Imports the DNA deformation state from a new protein-DNA complex structural model, and adds it to the local database for use in threaDNA.') #program parser and help
    parser.add_argument('name', metavar='prot_name', type=str,help='Name of the protein. Caution: if several models will be analyzed for this protein, always use the same name/case.')
    parser.add_argument('f',metavar='input',type=str,help='Input of the DNA base coordinates: (i) either a PDB ID (http://www.rcsb.org/pdb/home/home.do) or (ii) a file containing the coarse-grained coordinates of the DNA within the complex (.out from 3DNA or the webserver http://w3dna.rutgers.edu, or .lis from Curves+) or (iii) the .pdb atomic coordinates of the protein-DNA complex if the software x3dna is installed on the computer and accessible in the path from the current directory.')
    parser.add_argument("-p","--pdb",type=str,action='store', help="PDB file corresponding to the coordinates given in input. This allows threaDNA to extract additional information on the deformation model, in particular the structure resolution.")
    parser.add_argument("-r","--ref",type=int,action="store",help="Position of the reference nucleotide in this structure (default = central nucleotide). This is useful when combining/comparing models where the protein binds different numbers of DNA bases: you should then set a reference to align the binding profiles, for instance the index of a basepair that contacts a given protein residue. Typically, for the nucleosome it is the dyad basepair.")
    args=parser.parse_args()

    main(args.name,args.f,args.pdb,args.ref)
