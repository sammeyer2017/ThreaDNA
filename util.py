# coding: utf8

from globvar import *
import numpy as np
import sys
import string

def read_fasta(fa): #reads input fasta file to dictionary (seq_name:sequence)
    seq,seqn=[],[]
    try:
        f=open(fa,"r")
    except:
        sys.exit("File not found : "+fa)
    i=-1
    for l in f:
        if l[0]=='>':
            i+=1
            seqn.append(l[1:-1].replace(" ","_")) #replace whitespaces by underscores for further data utilisation (plot)
            seq.append("")
        else:
            if i<0:
               sys.exit("Incorrect FASTA file "+fa)
            seq[i]+=l[:-1].upper()
    f.close()
    print("sequences loaded : "+",".join(seqn))
    return np.array(seq),np.array(seqn)

def read_fasta2(fa): #reads input fasta file to dictionary (seq_name:sequence)
    seq,seqn,si=[],[],[]
    try:
        f=open(fa,"r")
    except:
        sys.exit("File not found : "+fa)
    i=0
    l=f.read().split("\n")
    f.close()
    while i<len(l):
        if l[i] and l[i][0]=='>':
            seqn.append(l[i][1:].replace(" ","_")) #replace whitespaces by underscores for further data utilisation (plot)
            si.append(i)
        i+=1
    si.append(i)
    i=0
    while i<len(si)-1:
        seq.append("".join(l[si[i]+1:si[i+1]]).upper())
        i+=1
    print("sequences loaded : "+",".join(seqn))
    return np.array(seq),np.array(seqn)


def read_params(p): #reads parameters file (.cnf)
    try:
        f=open(p,"r")
        p=dict(l.split(':') for l in filter(None,f.read().split('\n'))) #construct a dictionary from parameters
        f.close()
    except:
        sys.exit("Incorrect parameters file : "+p)
    if len(p)<6:
        sys.exit("Missing parameters in "+p) #not enough basic parameters
    #try:
    #    prot=l[0].split(':')[1]
    #    struct=l[1].split(':')[1].split(',')
    #    ca=l[2].split(':')[1]
    #    sym=(l[3].split(':')[1]=="Yes")
    #    tmp=(l[4].split(':')[1]=="Yes")
    #    Tf=float(l[5].split(':')[1])
    #    if len(l)>6:
    #        T=float(l[6].split(':')[1])
    #    else:
    #        T=300
    #    if len(l)>7:
    #        pat=l[7].split(':')[1]
    #    else:
    #        pat=None
    #except:
    #    sys.exit("Incorrect parameters file")
    return p

def convert(seq,ca): #convert seq expression to regular expression
    if not seq: #no input pattern
        return None,[]
    try:
        return None,int(seq)-1+(ca=="ABC_i" or ca=="ABC_old_i")
    except:
        pass
    j=0
    pat=r'(?=(' #pattern to find overlapping sequences
    for i in seq:
        if i=='N' or i=='X': #any nucleotide
            j+=1
        else:
            if j!=0:
                pat+='[ATGCN]{'+str(j)+'}' #computes the total of "any nucleotide" length
                j=0
            if i=='R': #converts IUPAC code to regular expression
 	        pat+='[AG]'
            elif i=='Y':
 	        pat+='[CT]'
            elif i=='S':
 	        pat+='[GC]'
            elif i=='W':
 	        pat+='[AT]'
            elif i=='K':
 	        pat+='[GT]'
            elif i=='M':
 	        pat+='[AC]'
            elif i=='B':
 	        pat+='[CGT]'
            elif i=='D':
 	        pat+='[AGT]'
            elif i=='H':
 	        pat+='[ACT]'
            elif i=='V':
                pat+='[ACG]'
            elif i in 'ATGC':
                pat+=i
            else: #unidentified nucleotide in pattern -> error
                sys.exit("Unidentified nucleotide in pattern : "+i)
    if j!=0:
        pat+='[ATGCN]{'+str(j)+'}' #computes the total of "any nucleotide" length at the end of the seq
    return pat+'))',len(seq)-1+(ca=="ABC_i" or ca=="ABC_old_i")

def comp(pat): #complementary pattern (regular expression)
    return pat.translate(string.maketrans("ATCG","TAGC"))

def load_dnamodel(ca,out,T=None): #load DNA model data (q0 ad K)
    if ca=="ABC_i": #load appropriate data
        filename="param/stiff_ABC_1mus_intra_tri.npz"
        n=3
    elif ca=="ABC_s":
        filename="param/stiff_ABC_1mus_step_tetra.npz"
        n=4
    elif ca=="ABC_old_i": #load appropriate data
        filename="param/stiff_ABC_50ns_intra_tri.npz"
        n=3
    elif ca=="ABC_old_s":
        filename="param/stiff_ABC_50ns_step_tetra.npz"
        n=4
    elif ca=="NP":
        filename="param/stiff_NP_di.npz"
        n=2
    else: #unknown parameter
        sys.exit("Please choose DNA parameters between NP, ABC_i, ABC_s, ABC_old_i, ABC_old_s")
    if T: #load temperature-related complementary data
        T=int(T)
        if ca=="ABC_i" or ca=="ABC_old_i":
            tempname="param/stiff_temp_intra_mono.npz"
        else:
            tempname="param/stiff_temp_step_di.npz"
        dt=np.load(out+tempname)
    d=np.load(out+filename)
    seq=sequences(ca) #generate sequences list
    tmp={}
    ind=[] #index to be constructed
    m=[] #q0
    s=[] #K
    for i,se in enumerate(seq):
        if selfcompseq(se): #self-complementary-sequence
            mi=np.divide(d["m"][i]*[0,1,1,0,1,1],adim(ca)) #q0
            si=.5*(d["s"][i]+np.multiply(np.outer(reversemat, reversemat),d["s"][i])) #K
            if T: #variation of parameters if temperature !=300K
                ti=shortseqind(ca,se)
                mt=np.divide(dt["m"][ti]*[0,1,1,0,1,1],adim(ca))
                st=.5*(dt["s"][ti]+np.multiply(np.outer(reversemat, reversemat),d["s"][ti]))
                mi-=(T/300.-1.)*mt
                si-=(T/300.-1.)*st
            m.append(mi)
            s.append(si)
            ind.append(se)
        else:               
            mi=np.divide(d["m"][i],adim(ca)) #q0
            si=d["s"][i] #K
            if T: #variation of parameters if temperature !=300K
                ti=shortseqind(ca,se)
                if ti>0:
                    mt=np.divide(dt["m"][ti],adim(ca))
                    st=dt["s"][ti]
                    mi-=(T/300.-1.)*mt
                    si-=(T/300.-1.)*st
                else:
                    it=-1*ti
                    mt=np.divide(np.multiply(dt["m"][it],reversemat),adim(ca))
                    st=np.multiply(np.outer(reversemat, reversemat),dt["s"][it])
                    mi-=(T/300.-1.)*mt
                    si-=(T/300.-1.)*st
            m.append(mi)
            s.append(si)
            ind.append(se)
            # conjugate sequence
            mi=np.divide(np.multiply(d["m"][i],reversemat),adim(ca))
            si=np.multiply(np.outer(reversemat, reversemat),d["s"][i])
            if T: #variation of parameters if temperature !=300K
                ti=shortseqind(ca,se)
                if ti>0:
                    mt=np.divide(np.multiply(dt["m"][ti],reversemat),adim(ca))
                    st=np.multiply(np.outer(reversemat, reversemat),dt["s"][ti])
                    mi-=(T/300.-1.)*mt
                    si-=(T/300.-1.)*st
                    #print mi
                else:
                    it=-1*ti
                    mt=np.divide(dt["m"][it],adim(ca))
                    st=dt["s"][it]
                    #print "tot",(T/300.-1.)*mt
                    mi-=(T/300.-1.)*mt
                    si-=(T/300.-1.)*st
                    #print mi
            m.append(mi)
            s.append(si)
            ind.append(compl_string(se))
    ind=dict(zip(ind,range(len(ind)))) #index dictionary construction
    print("Loaded "+ca+" parameters")
    return np.reshape(m,(1,-1,6)),np.array(s,dtype=np.float32),ind,n #reshaping for further utilisation (q0,K,ind,n)

def read_q(p,s,ca,out,l,off,al): #reads protein-bound DNA parameters
    end="_s.dat" #step parameters
    if ca=="ABC_i" or ca=="ABC_old_i": #intra pamaeters
        end="_i.dat"
    try:
        f=open(out+"/".join(["structures",p,s,s])+end,"r") #search in added structures data
    except:
        sys.exit("File not found : "+"/".join(["structures",p,s,s])+end)
    q=np.reshape(np.array([map(float,li[:-1].split(" ")) for li in f]),(-1,6)) #q extraction
    if type(l)==int: #pattern length modifying q length #center by default if P-l is odd
        if off:
            off=int(off)
        else:
            off=(l+2-(ca=="ABC_i" or ca=="ABC_old_i"))/2
        if al-off+l>np.shape(q)[0] or off>al:
            sys.exit("Pattern of length "+str(l)+" beginning at position "+str(al-off)+" exceeds protein size")
        #print al-off,al-off+l
        q=q[al-off:al-off+l,:]
    q=np.divide(q[:,(3,4,5,0,1,2)],adim(ca)) #q conversion
    f.close()
    return np.reshape(q,(-1,1,6)) #reshaping q for further utilisation

def rev_q(q): #reverses q for symetrization
    return np.multiply(q[::-1],reversemat,dtype=np.float64)

def calc_E(q,q0,K): #energy calculation by nucleotide combination
    qi=q-q0 #modification induced by protein
    E=np.einsum('ijk,jkl,ijl->ij',qi,K,qi)/2 #optimal computation of qi.T*K*qi
    #print np.shape(qi),np.shape(K),np.shape(E)
    return np.array(E)



            

def writepwm(Emat,ind,filename):
    """
    Saves a PWM matrix into a textfile of PSSM format. 
    Caution: in the main ThreaDNA program, the saved matrix contains the ELASTIC ENERGY associated to each sequence, not the frequency. The unit is arbitrary. An energy scale must then be used to get an absolute frequency. 
    Params: Matrix of energy (nb_elements_in_prot x all_possible_nucl); dictonary of sequences_indexes; filename
    The lowest value of the matrix is arbitrarily set to 0. 
    """
    initseqs=ind.keys()
    inds=ind.values()
    # change order
    nis=sorted(initseqs)
    Emat-=np.min(np.ravel(Emat))
    f=open(filename,"w")
    f.write(">%s\n"%(filename.split("/")[-1].split(".")[0]))
    for i,s in enumerate(nis):
        f.write("%s\t["%s)
        for x in Emat.T[ind[s]]:
            f.write("\t%.4f"%x)
        f.write("\t]\n")
    f.close()
    return 0
