# coding: utf8

from globvar import *
#from util import *
import numpy as np
import os
from Bio import motifs


# --------- JASPAR delimiter
jdel="\t"


def get_aligned_sequences(aligned_sequence_file):
    """
    reads sequence file in fasta format or list of sequences, and returns the list of sequences with names
    """
    try: 
        seqs=SeqIO.parse(aligned_sequence_file,"fasta")
        li=[(s.Seq,s.name) for s in seqs]
    except:
        l=open(aligned_sequence_file,"r")
        li=[(x.split("\n")[0],'') for x in l.readlines()]
        l.close()
        return li

def writepwm(mat,ind,filename):
    """
    Saves a PWM matrix into a textfile of PSSM format.
    It may be an ENERGY matrix (as provided by ThreaDNA) or a PROBABILITY matrix (or whatever other quantity, log-odds etc but we don't use it for the moment). In the latter case, the sum of each column is 1. 
    Params: Matrix (nb_elements_in_prot x all_possible_nucl) or proba; dictonary of sequences_indexes; filename
    """
    initseqs=ind.keys()
    inds=ind.values()
    # change order
    nis=sorted(initseqs)
    #Emat-=np.min(np.ravel(Emat))
    f=open(filename,"w")
    f.write(">%s\n"%(filename.split("/")[-1].split(".")[0]))
    for i,s in enumerate(nis):
        f.write("%s\t["%s)
        for x in mat.T[ind[s]]:
            f.write("%s%.3f"%(jdel,x))
        f.write("%s]\n"%jdel)
    f.close()
    return 0


def load_pwm(filename):
    """
    Load generalized PWM from file in JASPAR format
    Caution: limited to a single-PWM file... cannot handle JASPAR files with multiple PWMs
    """
    a=open(filename,"r")
    ls=a.readlines()
    seqs=[]
    mat=[]
    seqn=len(ls)-1
    for l in ls:
        if l[0]!=">":
            q=l.split(jdel)
            le=len(q)-3
            seqs.append(q[0])
            if q[1]!="[":
                print "problem: wrong file format: line %s"%l
                return 1
            mat.append(q[2:le-1])
            if q[-1][0]!="]":
                print "problem: wrong file format: line %s"%l
                return 1    
    a.close()
    ind=dict([(b,i) for i,b in enumerate(seqs)])
    matrix=np.array(mat, dtype=float).T
    return ind, matrix


def energ_pwm_to_proba_PWM(Emat,ind,b=1.):
    """
    Takes the matrix of energies pwm, and computes a matrix of probabilities
    Emat is a matrix of energies, with a priori unknown energy scale. 
    b is the multiplicative constant that gives 1 kBT in the chosen energy scale (typically larger than 1)
    ind is the dictionary of correspondence index-sequence
    """
    pmat=np.exp(-b*Emat)
    pmat=np.divide(pmat,np.sum(pmat,axis=1,keepdims=True))
    return pmat


def dinuc_PWM_to_mononuc(pmat, ind):
    """
    takes in an indirect readout probability matrix for dinucs, and returns the mononuc matrix 
    the sequence dictionary for mononucs is fixed by convention in the order A, C, G, T
    """
    # simplify to mononuc
    seqs=ind.keys()
    inds=ind.values()
    nucid=(len(seqs[0])-1)/2 # index of the mononuc in sequence
    mmat=np.zeros((len(pmat)+1,4))
    for si,l in enumerate((pmat.T)/2.):
        # find corresponding sequence
        s=seqs[inds.index(si)]
        mmat[:(-1),bases.index(s[nucid])]+=l
        mmat[1:,bases.index(s[nucid+1])]+=l
    # add a half value to the first and last dinucs, since these were undercounted
    mmat[0]*=2
    mmat[-1]*=2
    # test if all mononucprobabilities have a sum equal to 1
    #if np.not_equal(np.sum(mmat,axis=1),1.*np.ones(len(mmat))).any():
    if np.max(np.sum(mmat,axis=1))>1.000001:
        print "problem: the mononucleotide probability matrix is not correct. sum of element-wise values: %s"%str(np.sum(mmat,axis=1))
    return mmat


def plot_PWM(filename,output=None):
    """
    Plots a mononuc probability matrix from a JASPAR file
    Caution: only ONE motif in the file !!!
    The output must be given as the name of a PDF file, otherwise automatic from the input name
    """
    if output is None:
        output=filename.split(".")[0]+".pdf"
    # Uses WebLogo from BioPython to plot the matrix in pdf format
    # Requires an internet connection
    fh = open(filename)
    for m in motifs.parse(fh, "jaspar"):
        try:
            m.weblogo(output,format="PDF")
        except:
            print "ERROR trying to plot the sequence logo using the online website WebLogo (http://weblogo.berkeley.edu/). A possible cause is the absence of a working internet connection, otherwise the PWM file %s may be corrupted. If you wish to plot the sequence logo, consider using the alternate website STAMP (http://www.benoslab.pitt.edu/stamp/)"%filename
    fh.close()
    return 0


def make_and_plot_proba_matrices_from_energ_PWM(filename, b=1., plot=True):
    """
    - reads an energy motif file coming from ThreaDNA, in JASPAR format
    - outputs the probability matrix for dinucs in JASPAR format as well as mononuc, and plots the latter
    """
    ind,Emat=load_pwm(filename)
    pmat=energ_pwm_to_proba_PWM(Emat,ind)
    mmat=dinuc_PWM_to_mononuc(pmat, ind)
    # export pmat and mmat in JASPAR format
    writepwm(pmat,ind,filename.split(".")[0]+"_proba.pwm")
    mind=dict([(b,i) for i,b in enumerate(bases)])
    mmf=filename.split(".")[0]+"_mono.pwm"
    writepwm(mmat,mind,mmf)
    # plot the matrix
    if plot:
        mmf2=filename.split(".")[0]+"_mono2.pwm"
        writepwm(100*mmat,mind,mmf2)
        plot_PWM(mmf2)
        #os.system("rm %s"%mmf2)
    return 0


def seq_indexes_from_ind(seqs, ind):
    """
    takes a list of sequences and the dictionary of indexes of all dinuc/tetranuc types, and returns a numpy array with occurences of each index in sequences
    """
    dinucl=len(ind.keys()[0])
    return np.array([[ind[s[j:j+dinucl]] for j in range(len(s)-dinucl+1)] for s in seqs],dtype=np.uint8)

def compute_energies_from_energ_mat(seq_inds, pmat, ind):
    """
    makes the computation of energy profile from a matrix: intermediate function
    """
    e=[]
    P=len(pmat) # size of matrix, = size of prot -1 for dinuc, and size of prot for mononuc
    for i,s in enumerate(seqs):
        E_tmp=np.zeros(len(seq_tmp)-P+1,dtype=np.float32)
        for j in range(P): #loop on protein positions
            E_tmp=np.add(E_tmp,pmat[j,seq_inds[j:len(s)-P+j+1]],dtype=np.float32)
        e.append(E_tmp)
    return e

def compute_energy_profiles_from_PWM(fastafile, pwmfile):
    ind, mat = load_pwm(pwmfile)
    seqs=get_aligned_sequences(fastafile)
    seq_inds=seq_indexes_from_ind(seqs, ind)
    compute_energies_from_energ_mat(seq_inds, pmat, ind)


def compute_direct_PWM(aligned_sequence_file, indirect_proba_PWM=None, indirect_energy_PWM=None, b=1., protsize=None, filename=None):
    """
    Takes in a fasta file containing the aligned sequences for PWM, and the indirect (dinuc) PWM. 
    The latter can be given either as an energy matrix such as given by the main ThreaDNA program, or as a probability matrix. In the former case, an energy scale (beta factor) can be given. 
    Computes the indirect (dinuc and mononuc) and direct (dinuc and mononuc) matrices, as well as an information file with e.g. the weight of indirect vs direct readout
    Optional: restrict protein size to avoid side effects. Normally, the indirect output file is the same as input, with restricted size
    CAUTION: the alignments MUST MATCH, i.e. the centers of the sequences must also be the center of the provided PWM. i.e. without protsize, all sequences must have the length of indirect_proba_PWM + 1. 
    """
    # get sequences
    seqs=get_aligned_sequences(aligned_sequence_file)
    print seqs
    # get probability PWM, or compute it from energy PWM
    if indirect_proba_PWM != None:
        ind, mat = load_pwm(indirect_proba_PWM)
    elif indirect_energy_PWM != None:
        ind, emat = load_pwm(indirect_energy_PWM)
        mat=energ_pwm_to_proba_PWM(emat,ind,b)
    else:
        print "You must provide either a probability matrix or an energy matrix"
        return 1
    totlen=len(mat)
    # output name
    if filename is None:
        filename=indirect_proba_PWM.split(".")[0]+"_with_"+aligned_sequence_file.split("/")[-1].split(".")[0]
    # -----------------------------------------
    # compute the matrix part and the sequence length that we must incorporate
    if protsize is None:
        protsize=totlen+1
    dinucl=len(ind.keys()[0])     # size of sequence dictionary
    matlen=protsize-1
    seqlen=protsize+dinucl-2
    # correction of sizes
    seqlist=[]
    for s in seqs:
        l=len(s)
        # check size
        if l<seqlen:
            print("ERROR: sequence %s is too short for required size given by protein size or by matrix size"%s)
            return 1
        mid=l/2.
        seqlist.append(s[(int(mid-seqlen/2.)):(int(mid+seqlen/2.))])
    mid=(protsize-1)/2.
    pmat_indir=mat[(int(mid-(protsize-1)/2.)):(int(mid+(protsize-1)/2.))]
    # -------------------------------------------
    # compute dinuc sequences
    inds_in_fasta=seq_indexes_from_ind(seqlist, ind) # inds_in_fasta has EXACTLY the same size as pmat
    freqmat=np.zeros((matlen,len(ind.keys())))
    for li in inds_in_fasta:
        for ix,x in enumerate(li):
            freqmat[ix,x]+=1
    freqmat/=float(len(seqlist))
    if np.max(np.sum(freqmat,axis=1))>1.000001:
        print "ERROR: probabilities do not count to 1 in matrix: %s"%str(freqmat)
    # freqmat is the observed PWM, with the sum of energies
    # probabilities are multiplied :
    pmat_dir=freqmat/pmat_indir
    pmat_dir/=np.sum(pmat_dir,axis=1,keepdims=True)
    pwm_dir=dinuc_PWM_to_mononuc(pmat_dir, ind)
    mind=dict([(b,i) for i,b in enumerate(bases)])
    # --------------------------------------------
    writepwm(pmat_indir, ind, filename+"_ind.pwm")
    writepwm(pmat_dir, ind, filename+"_dir_dinuc.pwm")
    writepwm(pwm_dir, mind, filename+"_dir.pwm")
    return 0


compute_direct_PWM("CRP_seqs.fasta", indirect_proba_PWM=None, indirect_energy_PWM="test_CRP_en.pwm", b=14., protsize=22, filename="bla")
    
    





#ind, mmat = load_pwm("test_CRP.pwm")
#make_and_plot_proba_matrices_from_energ_PWM("test_CRP.pwm", b=14.3, plot=True)


# see STAMP motifs for checking !!!
