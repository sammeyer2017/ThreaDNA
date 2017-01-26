#!/usr/bin/env python
# coding: utf8

import argparse
#from subprocess import call
import os
import sys
import numpy as np
from itertools import product
import time
from util import *
import re

def getref(prot,struct,ca):
    al,mi,si=[],[],[]
    f1=open(loc+"structures/list.dat","r")
    l1=np.array([li[:-1].split("\t") for li in f1.readlines()[1:]])
    l1=l1[l1[:,0]==prot,:] #filter by prot name
    f1.close()
    f2=open(loc+"structures/"+ca+"_distrib.dat","r")
    l2=np.array([li[:-1].split("\t") for li in f2.readlines()[1:]])
    l2=l2[l2[:,0]==prot,:] #filter by prot name
    f2.close()
    if len(l1)==0:
        print "ERROR: This protein name is not in database. Mind the case!"
        print "You can find the list of available structures in INSTALL_DIR/structures/list.dat"
        quit()
    i=0
    while i<len(struct):
        if struct[i].isdigit():
            try:
                al.append(l1[int(struct[i])-1,5].astype(int)) #take step before
                struct[i]=l1[int(struct[i])-1,2] #get structure name by order
            except:
                sys.exit("No structure numbered "+str(struct[i]))
        else:
            try:
                al.append(l1[l1[:,2]==struct[i],5].astype(int)) 
            except:
                sys.exit("No structure named "+str(struct[i]))
        mi.append(l2[l2[:,2]==struct[i],3].astype(float))
        si.append(l2[l2[:,2]==struct[i],4].astype(float))
        i+=1
    return struct,np.reshape(al,-1),np.array(mi),np.array(si)

def seq_ind(seq,ind,n,i): #returns list of indexes from subsequence
    return ind[seq[i:i+n]]

def E_seq(E,seq,seqn,ind,n,st,tmp,pat,P,off): #calculates E from short sequences list (lenght of the protein), or matching patterns
    Et,pos,k=[],[],[]
    i=0
    while i<len(seqn): #for each sequence
        pos_t,Es=[],[]
        s=seq[i]
        for m in re.finditer(pat,s[(n-1)/2:len(s)-(n-1)/2]): #find all matches along sequence (! very important -(n/2) and NOT -n/2)
            Es.append(np.trace(E[:,[ind[s[j:j+n]] for j in range(m.start(),m.start()+P)]],dtype=np.float32))
            pos_t.append(m.start()+(n-1)/2)
        if Es: #match found
            Et.append(Es)
            pos.append(pos_t)
            k.append(seqn[i])
        i+=1
    if tmp and Et: #write temp files by structures in current directory
        writeprofile(st,Et,P,k,off,n,pos=pos)
        print "Writing : "+str(time.time()-t)
    elif tmp:
        print "Warning : couldn't write tmp file with no matching data for "+st  
    return Et,pos,k #energy and position of sequences

'''
def E_tot(E,seq,seqn,ind,n,st,tmp,al,P): #calculate total energy by position along all sequence
    #seq_vec=np.vectorize(seq_ind,excluded=['seq','ind','n'])
    E_s=[]
    i=0
    seqm=[]
    while i<len(seqn):
        si=seq[i]
        if len(si)-n-P+2<0:
            i+=1
            continue
        print "Analyzing sequence "+seqn[i]+"\n"#+s+"\n"
        seqm.append(seqn[i])
        E_tmp=np.zeros(len(si)-n-P+2,dtype=np.float32) #initializes energies at 0
        seq_tmp=np.array([ind[si[j:j+n]] for j in range(len(si)-n+1)],dtype=np.uint8) #list of subsequences indices in calculated E array
        #print "Preprocessing : "+str(time.time()-t) #display preprocessing time for this structure (cumulated)
        j=0
        #f=open("test"+st+i+".txt","w")
        while j<P: #loop on protein positions
            #print "pos "+str(j+1)+"/"+str(P)
            E_tmp=np.add(E_tmp,E[j,seq_tmp[j:len(si)-n-P+j+2]],dtype=np.float32) #add partial energy for position j of the protein
            #print E_tmp
            #f.write(str(E[j,seq_tmp[j:len(seq[i])-n-P+j+2]][0])+"\t"+str(j)+"\t"+str(P-j)+"\n")
            j+=1
        #f.close()
        #print "Energy calculation : "+str(time.time()-t) #display energy calculation time for this structure (cumulated)
        #print np.mean(E),np.mean(E_tmp),np.mean(E_tmp)/np.mean(E)
        E_s.append(E_tmp)
        i+=1
    #print E_s
    if tmp: #write temp files by structures in current directory 
        writetmp(st,E_s,P,seqm,al,n) #alignment parameter used as offset
        #print "Writing : "+str(time.time()-t) #display temp writing time for this structure (cumulated)
    return E_s,seqm
'''



def E_tot(E,seqs_tmp,seqn,ind,n,st,tmp,al,P): #calculate total energy by position along all sequence
    E_s=[]
    i=0
    seqm=[]
    while i<len(seqn):
        seq_tmp=seqs_tmp[i]
        if len(seq_tmp)-P+1<0:
            i+=1
            continue
        print "Analyzing sequence "+seqn[i]+"\n"
        seqm.append(seqn[i])
        E_tmp=np.zeros(len(seq_tmp)-P+1,dtype=np.float32) #initializes energies at 0
        #seq_tmp=np.array([ind[si[j:j+n]] for j in range(len(si)-n+1)],dtype=np.uint8) #list of subsequences indices in calculated E array
        #print "Preprocessing : "+str(time.time()-t) #display preprocessing time for this structure (cumulated)
        j=0
        #f=open("test"+st+i+".txt","w")
        while j<P: #loop on protein positions
            #print "pos "+str(j+1)+"/"+str(P)
            E_tmp=np.add(E_tmp,E[j,seq_tmp[j:len(seq_tmp)-P+j+1]],dtype=np.float32) #add partial energy for position j of the protein
            #print E_tmp
            #f.write(str(E[j,seq_tmp[j:len(seq[i])-n-P+j+2]][0])+"\t"+str(j)+"\t"+str(P-j)+"\n")
            j+=1
        #f.close()
        #print "Energy calculation : "+str(time.time()-t) #display energy calculation time for this structure (cumulated)
        #print np.mean(E),np.mean(E_tmp),np.mean(E_tmp)/np.mean(E)
        E_s.append(E_tmp)
        i+=1
    #print E_s
    if tmp: #write temp files by structures in current directory 
        writeprofile(st,E_s,P,seqm,al,n) #alignment parameter used as offset
        #print "Writing : "+str(time.time()-t) #display temp writing time for this structure (cumulated)
    return E_s,seqm

def writeprofile(name,E,P,seqn,off,n,pos=None,start=None): #write temp files by structures in current directory
    if isinstance(start,int):
        startoff=start
    else:
        startoff=[off for se in seqn]
    if type(P)==list:
        P=max(P)
    if len(seqn)==1:   # only one sequence in the fasta file
        k=0
        # simple case: use numpy.savetxt
        se=seqn[0][-24:]   # chromosome sequence: keep only 24 last chars to make a small file
        ee=E[k]
        if pos:
            posraw=np.array(pos[k],dtype=int)
            posref=posraw+off
            pos1=posraw+1
            posP=posraw+P+1-n%2
        else:
            posraw=np.arange(len(ee))
            posref=posraw+startoff[k]+(n-1)/2
            pos1=posraw+(n+1)/2
            posP=posraw+P+n/2
        if float(np.__version__[:3])>=1.7:
            #print os.getcwd().decode('utf8')
            #print name
            #print curdir
            np.savetxt(name,np.transpose((posref,posref,E[0],pos1,posP)), fmt=se+'\t%d\t%d\t%.2f\t%d\t%d',header='track type=bedGraph name=\"binding free energy (a. u.)\"',comments="")
        else:
            np.savetxt(name, np.transpose((posref,posref,E[0],pos1,posP)), fmt=se+'\t%d\t%d\t%.2f\t%d\t%d')#,header='track type=bedGraph name=\"raw binding free energy (k_B T)\"')
            os.system(r"sed -i -e '1itrack type=bedGraph name=\"binding free energy (a. u.)\"\' %s"%name)
    else:
        # several sequences
        f=open(name,"w")
        f.write('track type=bedGraph name=\"raw binding free energy (k_B T)\"\n') #bed header
        k=0
        while k<len(seqn):
            j=0
            se=seqn[k][-24:]
            ee=E[k]
            if pos:
                posraw=np.array(pos[k],dtype=int)
                posref=posraw+off
                pos1=posraw+1
                posP=posraw+P+1-n%2
            else:
                posraw=np.arange(len(ee))
                posref=posraw+startoff[k]+(n-1)/2
                pos1=posraw+(n+1)/2
                posP=posraw+P+n/2
            while j<len(ee): #writes sequence name, start position, end, position and energies for each energy calculated
                f.write("%s\t%d\t%d\t%0.2f\t%d\t%d\n"%(se,posref[j],posref[j],ee[j],pos1[j],posP[j]))
                j+=1
            k+=1                        
        f.close()

def writetmp_old(struct,E,P,seqn,off,n,pos=None): #write temp files by structures in current directory
    f=open(struct+".bed","w")
    f.write('track type=bedGraph name=\"raw binding free energy (k_B T)\"\n') #bed header
    k=0
    while k<len(seqn):
        j=0
        while j<len(E[k]): #writes sequence name, start position, end, position and energies for each energy calculated
            if pos:
                f.write(seqn[k]+"\t"+str(pos[k][j]+off)+"\t"+str(pos[k][j]+off)+"\t"+str(E[k][j])+"\t"+str(pos[k][j]+1)+"\t"+str(pos[k][j]+P+1-n%2)+"\n")
            else:
                f.write(seqn[k]+"\t"+str(j+off+(n-1)/2)+"\t"+str(j+off+(n-1)/2)+"\t"+str(E[k][j])+"\t"+str(j+(n+1)/2)+"\t"+str(j+P+n/2)+"\n")
            j+=1
        k+=1
    f.close()



def align(E,al,full,seqn): #align several sequences from different lengths
    i=0
    #if np.shape(E)[1]>1:
    start,end=[],[]
    #if np.shape(E)[1]==1:
    #    stat=np.reshape(stat,np.shape(E)[0])
    while i<len(E): #protein structures
        j=0
        start_tmp,end_tmp=[],[]
        while j<len(E[i]): #DNA sequences
            start_tmp.append(al[i])
            end_tmp.append(al[i]+len(E[i][j]))
            j+=1
        start.append(start_tmp)
        end.append(end_tmp)
        i+=1
    start,inds=np.amax(start,0),np.argmax(start,0) #max starting pos
    end,inde=np.amin(end,0),np.argmin(end,0) #min ending pos
    #else:
    #    print E[:,0]
    #    stat=min(len(E[:,0]))
    E2=[]
    i=0
    while i<len(E[0]): #DNA sequences
        j=0
        tmp=[]
        while j<len(E): #protein structures /!\ can take a lot of memory to copy tabs
            tmp.append(E[j][i][start[i]-al[j]:end[i]-al[j]])
            #if len(E[j][x][start[i]-al[j]:end[i]-al[j]])!=1:
                #print "error",len(E[j][x][start[i]-al[j]:end[i]-al[j]]),start,end #1 seq/struct but != length for sequences -> error in shape
            j+=1
        E2.append(np.array(tmp))
        i+=1  
    return E2,start,end #aligned inverted array, matching sequences list and ref positions

def writeweights(w,name,m,st,sym): #write the weights used for each structure
    f= open(name.split(".")[0]+".w","w")
    #f.write("Weight percentage by DNA sequence for each protein structure\n") #removed for easier data reading
    f.write("Sequence")
    for j in st:
            f.write("\t"+j) #structure names as header
            if sym:
                f.write("\t"+j+"_rev")
    f.write("\n")
    i=0
    while i<len(w):
        f.write(m[i]) #matching sequence name
        for j in w[i]: #structures weights for this sequence
            f.write("\t"+str(j)) 
        f.write("\n")
        i+=1
    f.close()

def writenorm_old(y,m,n,name,P,pos,off,struct,sym,al,start,end):
    i=0
    while i<len(P)/(1+sym): #protein structure
        j=0
        f=open(struct[i]+"_norm.bed","w")
        f.write('track type=bedGraph name=\"normalized binding free energy (a.u.)\"\n') #bed header
        if sym:
            f2=open(struct[i]+"_rev_norm.bed","w")
            f2.write('track type=bedGraph name=\"normalized binding free energy (a.u.)\"\n') #bed header
        while j<np.shape(y)[0]:  #DNA sequence
            k=0
            while k<end[j]-start[j]:  #position on sequence
                if pos:
                    f.write(m[j]+"\t"+str(pos[j][k]+off)+"\t"+str(pos[j][k]+off)+"\t"+str(y[j][i*(sym+1),k])+"\t"+str(pos[j][k]+1)+"\t"+str(pos[j][k]+P[i*(sym+1)]+1)+"\n")
                else:
                    f.write(m[j]+"\t"+str(start[j]+k+(n-1)/2)+"\t"+str(start[j]+k+(n-1)/2)+"\t"+str(y[j][i*(sym+1),k])+"\t"+str(start[j]-al[i*(sym+1)]+k+(n+1)/2)+"\t"+str(start[j]-al[i*(sym+1)]+k+P[i]+(n+1)/2)+"\n")
                if sym:
                    if pos:
                        f2.write(m[j]+"\t"+str(pos[j][k]+off)+"\t"+str(pos[j][k]+off)+"\t"+str(y[j][2*i+1,k])+"\t"+str(pos[j][k]+1)+"\t"+str(pos[j][k]+P[i]+1)+"\n")
                    else:
                        f2.write(m[j]+"\t"+str(start[j]+k+(n-1)/2)+"\t"+str(start[j]+k+(n-1)/2)+"\t"+str(y[j][2*i+1,k])+"\t"+str(start[j]-al[2*i+1]+k+(n+1)/2)+"\t"+str(start[j]-al[2*i+1]+k+P[i]+(n+1)/2)+"\n")
                k+=1
            j+=1
        f.close()
        if sym:
            f2.close()
        i+=1

def writeprofile_old(E,m,n,name,P,pos,off,start): #writes the final profile regrouping all structures for the protein
    if len(E)!=len(m): #different data sequences and matching sequences number
        sys.exit("Error while constructing global profile for several FASTA sequences")
    print "Writing profile in "+name+".bed"
    f=open(curdir+"/"+name+".bed","w")
    f.write('track type=bedGraph name=\"combined binding free energy (AU)\"\n') #bed header
    k=0
    if type(P)==list:
        P=max(P)
    while k<len(m): #DNA sequences
        j=0
        while j<len(E[k]): #positions on DNA, writes sequence name, start and stop positions, energy
            if pos:
                f.write(m[k]+"\t"+str(pos[k][j]+off)+"\t"+str(pos[k][j]+off)+"\t"+str(E[k][j])+"\t"+str(pos[k][j]+1)+"\t"+str(pos[k][j]+P+1-n%2)+"\n")
            else:
                f.write(m[k]+"\t"+str(j+start[k]+(n-1)/2)+"\t"+str(j+start[k]+(n-1)/2)+"\t"+str(E[k][j])+"\t"+str(j+(n+1)/2)+"\t"+str(j+P+n/2)+"\n")
            j+=1
        k+=1
    f.close()

def writeseqE(E,seq,name):
    f=open(name.split(".")[0]+"_seq.bed","w")
    f.write('track type=bedGraph name=\"sequence binding free energy (AU)\"\n') #bed header
    i=0
    while i<len(seq):
        f.write(seq[i]+"\t"+str(E[i])+"\n")
        i+=1
    f.close()

def profile(E,p,m,n,name,P,pos,mi,s,al,start,end): #calculates global profile from structures profiles
    sym=(p["Symmetrization"]=="Yes")
    try:
        b=1/float(p["Temperature factor"])
    except:
        b=1.
    E2,w,z=[],[],[]
    i=0
    while i<len(E):   #DNA sequences
        if len(mi)==(1+sym):
            si=np.array(s,dtype=np.float64) #unique structure correction
        else:
            si=np.array(s*np.std(mi/s,ddof=1,dtype=np.float64)) #unique structure correction
        if len(si.shape)==1:
            si=np.reshape(si,(-1,1))
        y=np.divide(E[i],si,dtype=np.float64) #corrected profiles by structure
        w.append(np.exp(-b/si)) #weights computation
        z.append(y)
        y=np.log(np.sum(np.exp(-b*y),0))/-b #final profile computation from structure profiles
        if np.isnan(y).any():
            sys.exit("nan found while calculating profile")
        E2.append(y)
        i+=1
    if p.get("Keep Weights")=="Yes" and w: #if weights are available
        writeweights(w,name,m,p["Structures"],p["Symmetrization"]=="Yes")  #save weights
    if type(P)==int:
         P=len(mi)*[P]
    #if p["Keep tmp"]=="Yes":
    #    writenorm(z,m,n,name,P,pos,p.get("Offset"),p["Structures"],p["Symmetrization"]=="Yes",al,start,end)
    if p.get("Sequence Energy")=="Yes":
        z=[np.log(sum(np.exp(-b*E)))/-b for E in E2]
        writeseqE(z,m,name)
    writeprofile(name,E2,P,m,p.get("Offset"),n,pos=pos,start=start) #write final profile
    # ---------------- old version
    # writeprofile(E2,m,n,name,P,pos,p.get("Offset"),start) #write final profile
    #print "Profile writing : "+str(time.time()-t) #displays the profile writing time (cumulated)

def main(fasta,params,output):
    t=time.time() #get program starting time
    print "Initializing..."
    # ---- Parameters
    p=read_params(params) #load input parameters (keys:values)
    p["Pattern"],P=convert(p.get("Pattern"),p["DNA parameter"]) #converts input pattern to regular expression
    p["Structures"],al,mi,si=getref(p["Protein"],p["Structures"].split(","),p["DNA parameter"])
    #if P:
        #al=len(p["Structures"])*[(P+2-(p["DNA parameter"]=="ABCi"))/2]
    q0,K,ind,n=load_dnamodel(p["DNA parameter"],loc,p.get("Temperature")) #load dna model parameters and index
    try:
        supercoiling=float(p["Superhelical density"])
    except:
        supercoiling=0.
    # --- sequence
    seq,seqn=read_fasta2(fasta) #load fasta squences from file
    if output!="None":
        name=output
    else:
        name=fasta.split("/")[-1].split(".")[0]+"_"+params.split("/")[-1].split(".")[0]+".bed" #output file
    seqs_tmp=np.array([[ind[s[j:j+n]] for j in range(len(s)-n+1)] for s in seq],dtype=np.uint8)
    # --------- HANDLING OF SUPERCOILING ----------------------- #
    if supercoiling != 0.:
        # sum up the sequence elements to get the global twist stiffness
        invKtorque=1/K[:,2,2]
        twists=q0[0,:,2]
        indexoccurences=np.bincount(np.ravel(seqs_tmp),minlength=len(invKtorque))
        ktot=1/np.sum(indexoccurences*invKtorque)
        thetatot=sum(indexoccurences*twists)
        # compute the suitable torque
        gamma=ktot*supercoiling*thetatot
        # for each sequence type, compute the displacement of equilibrium position
        Kinv=np.linalg.inv(K)
        dq0=gamma*Kinv[:,:,2]
        # TEST : the displacement should affect mostly the twist parameter ! 
        # print dq0/q0
        # update the value of q0
        q0+=dq0
    # ---------------------------------------------------------- # 
    #print "Initialization : "+str(time.time()-t) #displays initialization time
    #if not os.path.exists("/tmp/"+name): #make working directory in /tmp
    #    os.makedirs("/tmp/"+name)
    #os.chdir("/tmp/"+name) #going to tmp newly created dir
    if p["Keep tmp"]=="Yes" and not os.path.exists(name+"_tmp"):
        os.mkdir(name+"_tmp")
    Et,a_rev=[],None
    pos,pos_r=0,0
    for i,s in enumerate(p["Structures"]): #calculates E for each structure
        q=read_q(p["Protein"],s,p["DNA parameter"],loc,P,p.get("Offset"),al[i]) #unit conversion and data completion for q
        E=calc_E(q,q0,K) #energy for each set of nucleotides
        try:
            p["Offset"]=int(p["Offset"])
        except:
            p["Offset"]=(np.shape(q)[0]+2-(p["DNA parameter"]=="ABCi"))/2
        if p["Pattern"]: #uses pattern finding and short sequences energy calculating function
            print "Calculating Energy on patterns for structure "+s
            Etmp,pt,sm=E_seq(E,seq,seqn,ind,n,name+"_tmp/"+s+".bed",p["Keep tmp"]=="Yes",p["Pattern"],P,p["Offset"])
            if Etmp:
                Et.append(Etmp) #energies for each subsequence
                seqm=sm
                pos=pt #positions of subsequences in global sequences
        else:
            print "Calculating Energy for structure "+s
            if type(P)==list:
                P.append(np.shape(q)[0]) #get size of the protein from q
                Etmp,sm=E_tot(E,seqs_tmp,seqn,ind,n,name+"_tmp/"+s+".bed",p["Keep tmp"]=="Yes",al[i],P[i]) #calculate total energy along the whole sequences
                Et.append(Etmp) #calculate total energy along the whole sequences
                seqm=sm
            else:
                Etmp,sm=E_tot(E,seqs_tmp,seqn,ind,n,name+"_tmp/"+s+".bed",p["Keep tmp"]=="Yes",p["Offset"],P)
                Et.append(Etmp) #calculate total energy along the whole sequences
                seqm=sm
        if p["Symmetrization"]=="Yes": #calculates total energy for symmetric sructure
            #np.savetxt(s+"_q.dat",np.reshape(q[al[i]-5:al[i]+5],(10,6)),delimiter="\t")
            q=rev_q(q) #reverses q
            #np.savetxt(s+"_revq.dat",np.reshape(q[P[i]-al[i]+2-n%2-5:P[i]-al[i]+2-n%2+5],(10,6)),delimiter="\t")
            E=calc_E(q,q0,K) #recalculates energy
            if p["Pattern"]: #uses pattern finding and short sequences energy calculating function
                print "Calculating Energy on patterns for symmetric structure "+s+"_rev"
                Etmp,pt,sm=E_seq(E,seq,seqn,ind,n,name+"_tmp/"+s+"_rev.bed",p["Keep tmp"]=="Yes",p["Pattern"],P,p["Offset"])
                if Etmp: 
                    Et.append(Etmp)
            else:
                print "Calculating Energy for symmetric structure "+s+"_rev"
                if type(P)==list:
                    Etmp,sm=E_tot(E,seqs_tmp,seqn,ind,n,name+"_tmp/"+s+"_rev.bed",p["Keep tmp"]=="Yes",P[i]-al[i]+2-n%2,P[i])
                    Et.append(Etmp) #calculate total energy along the whole sequences
                else:
                    Etmp,sm=E_tot(E,seqs_tmp,seqn,ind,n,name+"_tmp/"+s+"_rev.bed",p["Keep tmp"]=="Yes",P-p["Offset"]+2-n%2,P)
                    Et.append(Etmp)
    print "Calculation : "+str(time.time()-t) #displays the total calculation time
    #if len(Et)>1: #global profile calculation for several structures
    if type(P)==int:
        al=np.array([p["Offset"]]*len(p["Structures"]))
    if p["Symmetrization"]=="Yes":
        mi=np.reshape([x for pair in zip(mi,mi) for x in pair],-1)
        si=np.reshape([x for pair in zip(si,si) for x in pair],-1)
        al=np.reshape([x for pair in zip(al,P-al+2-n%2) for x in pair],-1)
    Et,start,end=align(Et,al,not p["Pattern"],seqn) #structures alignment
    profile(Et,p,seqm,n,name,P,pos,mi,si,al,start,end)
    #p.get("Keep Weights")=="Yes",p["Structures"],p["Offset"],p["Symmetrization"]=="Yes",p["Keep tmp"]=="Yes") #profile computation
    print "Global profile built !"
    #else:
    #    print "Warning : not enough data to build global profile, make sure you save the partial data"
    #os.chdir("/tmp") #going back to /tmp
    #a="mv "+name+" "+curdir.replace(" ","\ ")
    #if p["Keep tmp"]=="Yes":
    #    if not os.path.exists(curdir+"/"+name): 
    #        os.system(("mv "+name+" "+curdir.replace(" ","\ ")).encode(encoding="utf-8")) #moving temp dir in user cwd
    #    else:
    #        os.system(("mv "+name+"/* "+curdir.replace(" ","\ ")+"/"+name).encode(encoding="utf-8")) #moving only files in existing dir
    #else: 
    #    os.system("rm -r "+name) #cleaning tmp dir, warning for only one structure will remove the only partial data
    print "Execution time :"+str(time.time()-t) #prints total exeution time
    sys.exit(0)


loc=os.path.abspath(os.path.dirname(sys.argv[0])).decode('utf8')+u"/" #saving exec. dir
#curdir=os.getcwd().decode('utf8')#saving cwd
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Computes the DNA deformation energy profile for a DNA-binding protein model along a DNA sequence, taking into account the base-pair or base-pair-step deformations.') #argument parser and help
    parser.add_argument('fasta',type=str,help='FASTA file containing DNA sequence(s)')
    parser.add_argument('params', type=str,help='Parameter file describing the protein model and computation parameters')
    parser.add_argument("-o","--output",type=str,action="store",help="Name of .bed output file")
    #parser.add_argument('data',metavar='DNA_params', type=str,help='DNA parameters q0 and K in npz format')
    args=parser.parse_args()
    main(unicode(args.fasta),unicode(args.params),unicode(args.output)) #executes program

#os.system("make") #compiling c code
#t1=os.system("./main "+str(args.fasta)) #calling c program and storing return value
#t2=call(["./main", "test", "to", "main"])
#print "python got return value "+str(t1/256)
