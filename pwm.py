# coding: utf8

from globvar import *
#from util import *
import numpy as np
import os
from Bio import motifs,SeqIO


# --------- JASPAR delimiter
jdel="\t"


def get_aligned_sequences(aligned_sequence_file):
    """
    reads sequence file in fasta format or list of sequences, and returns the list of sequences with names
    """
    ext=aligned_sequence_file.split('.')[-1]
    if ext=="fa" or ext=="fasta" or ext=="fsa":
        seqs=list(SeqIO.parse(aligned_sequence_file,"fasta"))
        li=[str(s.seq) for s in seqs]
        names=[s.id for s in seqs]
    else:
        print "opened sequence file as a normal text file"
        l=open(aligned_sequence_file,"r")
        li=[x.split("\n")[0] for x in l.readlines()]
        names=["" for s in li]
        l.close()
    return li,names

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
            le=len(q)
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


def energ_pwm_to_proba_pwm(Emat,ind,b=1.):
    """
    Takes the matrix of energies pwm, and computes a matrix of probabilities
    Emat is a matrix of energies, with a priori unknown energy scale. 
    b is the multiplicative constant that gives 1 kBT in the chosen energy scale (typically larger than 1)
    ind is the dictionary of correspondence index-sequence
    """
    pmat=np.exp(-b*Emat)
    pmat=np.divide(pmat,np.sum(pmat,axis=1,keepdims=True))
    return pmat





def dinuc_pwm_to_mononuc(pmat, ind):
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


def make_and_plot_proba_matrices_from_energ_pwm(filename, b=1., plot=True):
    """
    - reads an energy motif file coming from ThreaDNA, in JASPAR format
    - b gives the energy scale to go to kT
    - outputs the probability matrix for dinucs in JASPAR format as well as mononuc, and plots the latter
    """
    ind,Emat=load_pwm(filename)
    pmat=energ_pwm_to_proba_pwm(Emat,ind,b)
    mmat=dinuc_pwm_to_mononuc(pmat, ind)
    # export pmat and mmat in JASPAR format
    writepwm(pmat,ind,filename.split(".")[0]+"_proba.pwm")
    mind=dict([(ba,i) for i,ba in enumerate(bases)])
    mmf=filename.split(".")[0]+"_mono.pwm"
    writepwm(mmat,mind,mmf)
    # plot the matrix
    if plot:
        mmf2=filename.split(".")[0]+"_mono2.pwm"
        writepwm(100*mmat,mind,mmf2)
        plot_pwm(mmf2)
        #os.system("rm %s"%mmf2)
    return 0


def seq_indexes_from_ind(seqs, ind):
    """
    takes a list of sequences and the dictionary of indexes of all dinuc/tetranuc types, and returns a numpy array with occurences of each index in sequences
    """
    dinucl=len(ind.keys()[0])
    return np.array([[ind[s[j:j+dinucl]] for j in range(len(s)-dinucl+1)] for s in seqs],dtype=np.uint8)


def compute_pwm_from_sequences(aligned_sequence_file, outfile=None):
    if outfile==None:
        outfile=aligned_sequence_file.split(".")[0]+".pwm"
    seqs,names=get_aligned_sequences(aligned_sequence_file)
    mind=dict([(ba,i) for i,ba in enumerate(bases)])
    inds_in_fasta=seq_indexes_from_ind(seqs, mind)
    freqmat=np.zeros((len(inds_in_fasta[0]), 4))
    for li in inds_in_fasta:
        for ix,x in enumerate(li):
            freqmat[ix,x]+=1
    writepwm(freqmat/float(len(inds_in_fasta)),mind,outfile)
    print freqmat/float(len(inds_in_fasta))
    return freqmat/float(len(inds_in_fasta))


def compute_energies_from_energ_mat(seq_inds, emat):
    """
    makes the computation of energy profile from a matrix: intermediate function
    """
    e=[]
    P=len(emat) # size of matrix, = size of prot -1 for dinuc, and size of prot for mononuc
    for i,s in enumerate(seq_inds):
        E_tmp=np.zeros(len(s)-P+1,dtype=np.float32)
        for j in range(P): #loop on protein positions
            E_tmp=np.add(E_tmp,emat[j,s[j:len(s)-P+j+1]],dtype=np.float32)
        e.append(E_tmp)
    return e


def firstind_float(pmat, ind):
    """
    gives the first index of the energy value computed along the sequence... in 1-index !!!
    caution: returns a float
    """
    return len(pmat)/2.+len(ind.keys()[0])/2


def writeprofile_simple(filename,E,seqnames,firstind): #write temp files by structures in current directory
    """
    Write energy profile as a bed file
    - output filename
    - list of arrays of energies for different sequences
    - seqnames 
    - first index of each energy value: a floating number or an integer
    """
    if len(seqnames)==1:   # only one sequence in the fasta file: use Numpy for quicker export if the sequence is very long
        k=0
        # simple case: use numpy.savetxt
        se=seqnames[0]   # chromosome sequence: keep only 24 last chars to make a small file
        ee=E[k]
        posref=np.arange(len(ee))+firstind
        if float(np.__version__[:3])>=1.7:
            # directly export file with header
            np.savetxt(filename,np.transpose((posref-0.5,posref+0.5,E[0])), fmt=se+'\t%.1f\t%.1f\t%.2f',header='track type=bedGraph name="binding free energy (a. u.)"',comments="")
        else:
            # export file and add header afterwards
            np.savetxt(filename, np.transpose((posref-0.5,posref+0.5,E[0])), fmt=se+'\t%.1f\t%.1f\t%.2f')
            os.system(r"sed -i -e '1itrack type=bedGraph name=\"binding free energy (a. u.)\"' %s"%filename)
    else:
        # several sequences
        f=open(filename,"w")
        f.write('track type=bedGraph name="raw binding free energy (a. u.)"\n') #bed header
        k=0
        while k<len(seqnames):
            j=0
            se=seqnames[k]
            ee=E[k]
            posref=np.arange(len(ee))+firstind
            while j<len(ee): #writes sequence name, start position, end, position and energies for each energy calculated
                f.write("%s\t%.1f\t%.1f\t%0.2f\n"%(se,posref[j]-0.5,posref[j]+0.5,ee[j]))
                j+=1
            k+=1                        
        f.close()
    return 0


def compute_energy_profiles_from_pwm(fastafile, energ_pwmfile=None, proba_pwmfile=None, bedname=None):
    """ 
    Main function to compute an energy profile from a pwm file... generally an energy PWM file coming from ThreaDNA, but a probability PWM is also possible. In that case, the energy profile is given in kT
    """
    if energ_pwmfile != None:
        ind, mat = load_pwm(energ_pwmfile)
    elif proba_pwmfile != None:
        ind, pmat = load_pwm(proba_pwmfile)
        mat=-np.log(pmat+10**-8)
        print "Energy profile provided in kT unit"
    else:
        print "You must provide either a probability matrix or an energy matrix for the indirect contribution"
        return 1
    if bedname is None:
        bedname=fastafile.split(".")[0]+'_'+fastafile.split("/")[-1].split(".")[0]+".bed"
    #ind, mat = load_pwm(energ_pwmfile)
    #print np.shape(mat)
    seqs,names=get_aligned_sequences(fastafile)
    #print seqs,names
    seq_inds=seq_indexes_from_ind(seqs, ind)
    e=compute_energies_from_energ_mat(seq_inds, mat)
    fi=firstind_float(mat, ind)
    writeprofile_simple(bedname,e,names,fi)
    return 0




def compute_direct_pwm_from_dinuc_distributions_by_counting_dinuc_biases(aligned_sequence_file, indirect_proba_pwm=None, indirect_energy_pwm=None, b=1., protsize=None, filename=None):
    """
    Takes in a fasta file containing the aligned sequences for PWM, and the indirect (dinuc) PWM. 
    The latter can be given either as an energy matrix such as given by the main ThreaDNA program, or as a probability matrix. In the former case, an energy scale (beta factor) can be given. 
    Computes the indirect (dinuc and mononuc) and direct (dinuc and mononuc) matrices, as well as an information file with e.g. the weight of indirect vs direct readout
    THIS VERSION of the computation is based on a strong assumption: the dinucleotides are present in the list of sequences in proportions representative of their free energy values, and these proportions are used to estimate the relative energies, and thus of the direct contributions. Effectively, this requires a considerable number of known sites, larger than usually available. Also, we give the same weight to different sequences which have different affinities. In practice, the approach seems too demanding to be effective. 
    Optional: restrict protein size to avoid side effects. Normally, the indirect output file is the same as input, with restricted size
    CAUTION: the alignments MUST MATCH, i.e. the centers of the sequences must also be the center of the provided PWM. i.e. without protsize, all sequences must have the length of indirect_proba_PWM + 1. 
    """
    # get sequences
    seqs,names=get_aligned_sequences(aligned_sequence_file)
    print seqs
    # get probability PWM, or compute it from energy PWM
    if indirect_proba_pwm != None:
        ind, mat = load_pwm(indirect_proba_pwm)
    elif indirect_energy_pwm != None:
        ind, emat = load_pwm(indirect_energy_pwm)
        mat=energ_pwm_to_proba_pwm(emat,ind,b)
    else:
        print "You must provide either a probability matrix or an energy matrix for the indirect contribution"
        return 1
    totlen=len(mat)
    # output name
    if filename is None:
        filename=indirect_proba_pwm.split(".")[0]+"_with_"+aligned_sequence_file.split("/")[-1].split(".")[0]
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
    pwm_dir=dinuc_pwm_to_mononuc(pmat_dir, ind)
    mind=dict([(b,i) for i,b in enumerate(bases)])
    # --------------------------------------------
    writepwm(pmat_indir, ind, filename+"_ind.pwm")
    writepwm(pmat_dir, ind, filename+"_dir_dinuc.pwm")
    writepwm(pwm_dir, mind, filename+"_dir.pwm")
    return 0


def compute_direct_pwm_from_dinuc_distributions_by_reweighting_whole_sequences(aligned_sequence_file, indirect_proba_pwm=None, indirect_energy_pwm=None, b=1., indirect_seqsize=None, filename=None, percentile=90.):
    """
    Takes in a fasta file containing the aligned sequences for PWM, and the indirect (dinuc) PWM. 
    The latter can be given either as an energy matrix such as given by the main ThreaDNA program, or as a probability matrix. In the former case, an energy scale (beta factor) can be given. 
    Computes the indirect (dinuc and mononuc) and direct (dinuc and mononuc) matrices, as well as an information file with e.g. the weight of indirect vs direct readout
    IN THIS VERSION, we simply give a global weight to each sequence, computed from its indirect readout energy. This is not very precise, since we forget where indirect readout acts in the sequence, but the more precise version seemed to demanding. 
    Optional: restrict protein size for indirect readout to location where no direct contacts occur. Normally, the indirect output file is the same as input, with restricted size
    CAUTION: the alignments MUST MATCH, i.e. the centers of the sequences must also be the center of the provided PWM. i.e. without protsize, all sequences must have the length of indirect_proba_PWM + 1. 
    Percentile is the fraction of input sequences that we keep in the computation of the PWM. Indeed, false positives tend to get a strong deformation energy, and thus to be over-counted in the subsequent PWM... better throw them out by taking only e.g. 90% of the sequences
    """
    # get sequences
    seqs,names=get_aligned_sequences(aligned_sequence_file)
    #print seqs
    # get the indirect energy PWM in kt units
    if indirect_proba_pwm != None:
        ind, pmat = load_pwm(indirect_proba_pwm)
        emat=-np.log(pmat_indir)
    elif indirect_energy_pwm != None:
        ind, mat = load_pwm(indirect_energy_pwm)
        emat=mat*b
    else:
        print "You must provide either a probability matrix or an energy matrix for the indirect contribution"
        return 1
    totlen_indir=len(emat) 
    # output name
    if filename is None:
        filename=indirect_proba_pwm.split(".")[0]+"_with_"+aligned_sequence_file.split("/")[-1].split(".")[0]
    # -----------------------------------------
    # compute the matrix part and the sequence length that we must incorporate
    if indirect_seqsize is None:
        indirect_seqsize=totlen_indir+1
    dinucl=len(ind.keys()[0])     # size of sequence dictionary
    indir_matlen=indirect_seqsize-1
    seqlen=indirect_seqsize+dinucl-2
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
    mid=(indirect_seqsize-1)/2.
    emat_indir=emat[(int(mid-indir_matlen/2.)):(int(mid+indir_matlen/2.))]
    # -------------------------------------------
    # compute mononuc sequences corrected for whole energy
    mind=dict([(b,i) for i,b in enumerate(bases)])
    dinucinds_in_fasta=seq_indexes_from_ind(seqlist,ind) # inds_in_fasta has EXACTLY the same size as pmat
    es=np.array(compute_energies_from_energ_mat(dinucinds_in_fasta, emat_indir))[:,0] # normally each sequence has exactly ONE energy value because the sizes match exactly
    weights=np.exp(es-np.mean(es))
    limval=np.percentile(weights, percentile)
    correct_seqs=np.where(weights<limval)[0]
    print correct_seqs                        
    inds_in_fasta=seq_indexes_from_ind(seqs,mind)
    freqmat=np.zeros((len(seqs[0]),len(ind.keys())))
    for il,li in enumerate(inds_in_fasta):
        if il in correct_seqs:
            print li, weights[il]
            for ix,x in enumerate(li):
                freqmat[ix,x]+=weights[il]
    freqmat/=np.sum(weights[correct_seqs])
    if np.max(np.sum(freqmat,axis=1))>1.000001:
        print "ERROR: probabilities do not count to 1 in matrix: %s"%str(freqmat)
    # --------------------------------------------
    writepwm(emat_indir, ind, filename+"_ind_en.pwm")
    writepwm(freqmat, mind, filename+"_dir_proba.pwm")
    return 0









# compute_direct_pwm("CRP_seqs.fasta", indirect_proba_pwm=None, indirect_energy_pwm="test_CRP_en.pwm", b=14., protsize=22, filename="bla")
#compute_energy_profiles_from_pwm("crp_lindemose.fa", "crp_lindemose_CRP_en.pwm", bedname="lindemose_from_threadna.bed")
#compute_energy_profiles_from_pwm("crp_lindemose.fa", energ_pwmfile=None, proba_pwmfile="CRP_seqs.pwm", bedname="lindemose_from_rdbpwm.bed")





# -------------- EXECUTION
   
#loc=os.path.abspath(os.path.dirname(sys.argv[0])).decode('utf8')+u"/" #saving exec. dir
# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description='Computes the DNA deformation energy profile along a DNA sequence from a position weight matrix file in the JASPAR format, either classical (mononucleotides) or indirect (dinucleotides).') 
#     parser.add_argument('', type=str,help='Parameter file describing the protein model and computation parameters')
#     parser.add_argument("-s","--sequence",type=str,help='FASTA file containing DNA sequence(s)')
#     parser.add_argument("-o","--output",type=str,action="store",help="Name of .bed/.dpwm output file")
#     args=parser.parse_args()
#     main(unicode(args.params),unicode(args.sequence),unicode(args.output)) #executes program


#ind, mmat = load_pwm("test_CRP.pwm")
#make_and_plot_proba_matrices_from_energ_pwm("crp_lindemose_CRP_en.pwm", b=20., plot=True)


# see STAMP motifs for checking !!!
compute_direct_pwm_from_dinuc_distributions_by_reweighting_whole_sequences("CRP_seqs.txt", indirect_proba_pwm=None, indirect_energy_pwm="crp_lindemose_CRP_en.pwm", b=10., indirect_seqsize=8, filename="bla")
