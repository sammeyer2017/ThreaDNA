# coding: utf8
npar=41
parameters=['shear', 'stretch', 'stagger', 'buckle', 'propel', 'opening', 'xdisp', 'ydisp', 'inclin', 'tip', 'ax-bend', 'shift', 'slide', 'rise', 'tilt', 'roll', 'twist', 'h-ris', 'h-twi', 'alphaW', 'betaW', 'gammaW', 'deltaW', 'epsilW', 'zetaW', 'chiW', 'phaseW', 'ampW', 'alphaC', 'betaC', 'gammaC', 'deltaC', 'epsilC', 'zetaC', 'chiC', 'phaseC', 'ampC', 'minw', 'mind', 'majw', 'majd']

#par=parameters[11:17]
par=parameters[14:17]+parameters[11:14]   #angles first
intrapar=parameters[3:6]+parameters[:3]

global bases, alldinucs, adinucs, cdinucs, nextdinucs, nndinucs, monos
                  
bases=['A','C','G','T']

monos=["A","C"]

alldinucs=["AA", "AC", "AG", "AT", "CA", "CG", "GA", "GC", "GG", "TA"] #dinucleotides
adinucs=["AA", "AC", "AG", "AT","TA"]
cdinucs=["CA", "CG", "GA", "GC", "GG"]

# dinucs=["AA", "AC", "AG", "AT", "CA", "CG", "GA", "GC", "GG", "TA"]

nextdinucs=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "ATA", "ATC", "ATG", "CAA", "CAC", "CAG", "CCA", "CCC", "CCG", "CGA", "CGC", "CTA", "CTC", "GAA", "GAC", "GCA", "GCC", "GGA", "GTA", "TAA", "TCA"] #trinucleotides

#tetranucleotides
nndinucs=["AAAA", "AAAC", "AAAG", "AAAT", "AACA", "AACC", "AACG", "AACT", \
"AAGA", "AAGC", "AAGG", "AAGT", "AATA", "AATC", "AATG", "AATT", \
"ACAA", "ACAC", "ACAG", "ACAT", "ACCA", "ACCC", "ACCG", "ACCT", \
"ACGA", "ACGC", "ACGG", "ACGT", "ACTA", "ACTC", "ACTG", "AGAA", \
"AGAC", "AGAG", "AGAT", "AGCA", "AGCC", "AGCG", "AGCT", "AGGA", \
"AGGC", "AGGG", "AGTA", "AGTC", "AGTG", "ATAA", "ATAC", "ATAG", \
"ATAT", "ATCA", "ATCC", "ATCG", "ATGA", "ATGC", "ATGG", "ATTA", \
"ATTC", "ATTG", "CAAA", "CAAC", "CAAG", "CACA", "CACC", "CACG", \
"CAGA", "CAGC", "CAGG", "CATA", "CATC", "CATG", "CCAA", "CCAC", \
"CCAG", "CCCA", "CCCC", "CCCG", "CCGA", "CCGC", "CCGG", "CCTA", \
"CCTC", "CGAA", "CGAC", "CGAG", "CGCA", "CGCC", "CGCG", "CGGA", \
"CGGC", "CGTA", "CGTC", "CTAA", "CTAC", "CTAG", "CTCA", "CTCC", \
"CTGA", "CTGC", "CTTA", "CTTC", "GAAA", "GAAC", "GACA", "GACC", \
"GAGA", "GAGC", "GATA", "GATC", "GCAA", "GCAC", "GCCA", "GCCC", \
"GCGA", "GCGC", "GCTA", "GGAA", "GGAC", "GGCA", "GGCC", "GGGA", \
"GGTA", "GTAA", "GTAC", "GTCA", "GTGA", "GTTA", "TAAA", "TACA", \
"TAGA", "TATA", "TCAA", "TCCA", "TCGA", "TGAA", "TGCA", "TTAA"]


global xc,xg,seqs,datadir,savedatadir,plotdir


def longc(ca): #unused function
    if ca is "i":
        return "intra"
    else:
        return "step"

def shortc(ca): #unused function
    if ca is "intra":
        return "i"
    else:
        return "s"


#######################################################

# Some useful functions
def oligoseq(se): return "GC"+se[2:]+se+se+se+"GC"

def compl_base(b): #bases complement
    if b=="A":
        return "T"
    if b=="T":
        return "A"
    if b=="C":
        return "G"
    if b=="G":
        return "C"

def compl_string(seq): #sequence reverse complement
    se=[compl_base(x) for x in seq]
    se.reverse()
    return ''.join(se)

def selfcompseq(seq): #test if sequence is self-complementary
    return compl_string(seq)==seq


# list of substring positions in string
def str_occur(str,substr):
    l=len(substr)
    res=[]
    for i in range(len(str)-l+1):
        if substr==str[i:i+l]:
            res.append(i)
    return res

# list of tetramers corresponding to a dinuc
def tetramers(dinuc):
    return filter(lambda x: (x[1]==dinuc[0] and x[2]==dinuc[1],nndinucs))

def adim(case): #convert parameters units
    if case=="intra" or case=="ABC_i" or case=="ABC_old_i":
        return [11.62, 9.37, 4.53, 0.302, 0.117, 0.425]
    elif case=="step" or case=="ABC_old_s" or case=="ABC_s" or case=="NP":
        return [4.46, 7.10, 6.76, 0.712, 0.712, 0.345]


'''
########################################################
# BUILD SEVERAL DIRECTORIES CORRESPONDING TO DIFFERENT SITUATIONS
# IN ALL CASE, WE KEEP DIFFERENT VALUES OF KG300, BUT WE TAKE A VARIABLE NUMBER OF TM
##  1. INTRA
#  intra_complete : compute separately all dinuc: this is justified because the dinuc correspond to different neighboring sequences for the considered bp
# intra_mononuc : compute separately the AT and the CG (2 dinucleotides) : this is justified because the values are close
# intra_neutral : compute ALL together : a single value of Tm for all sequences

## 2. BPSTEP PARAMS
# step_complete : separately all tetramers
# step_neutral : all together


# Note : to switch from the complete set to other sets
# nuc i : divide the K_i by Ko_i : scales ALL dinucs on the same level
# linreg on this data : yields a Ko_glob and Tm_glob
# then the Tm_glob is for all, and Ko_i^new is Ko_i*Ko_glob


global dirs
dirs=["a"]#,"e"]

global allcases,allcasestot
allcases=["ABC_s","ABC_i","NP","ABC_old_s","ABC_old_i"]
allcasestot=["ABCs","ABCi","t","NP","NI"]


def pars(case): #unused function
    if case is "intra": 
        return intrapar
    else:
        return par


def unit(ipar): #unused function
    if ipar<=2:
        return r'$\deg$'
    else:
        return r'$\AA$'


def cov_unit(ipar1,ipar2): #unused function
    return r'$[\sigma_{%d} \sigma_{%d}]$' % (ipar1+1,ipar2+1)

def leg(ipar): #unused function
    return par[ipar]+" ("+ unit(ipar)+")"

def stleg(ipar1,ipar2): #unused function
    return r'$[k_B T_0]/[\sigma_{%d} \sigma_{%d}]$' % (ipar1+1,ipar2+1)
    
    ## if ipar1<=2 and ipar2<=2:
    ##     return r'$k_B T_0 /\deg^2$'
    ## elif ipar1>=3 and ipar2>=3:
    ##     return r'$k_B T_0 /\AA^2$'
    ## else:
    ##     return r'$k_B T_0 / (\AA.\deg)$'

'''


def sequences(case): #find approprioate nucleotides table
    if case is "intra" or case=="ABC_i" or case=="ABC_old_i":
        return nextdinucs
    if case=="step" or case=="ABC_s" or case=="ABC_old_s":
        return nndinucs
    if case=="NP":
        return alldinucs

def allsequences(case): #returns all nucleotides combinations
    return sequences(case)+filter(not selfcomplseq,map(compl_seq,sequences(case)))
    
def othersequences(case):
    return filter(allsequences(case),lambda se: se not in sequences(case))

# give the index of a longseq in the list of short seqs
# here longseq is trinuc for intra, tetra OR di for step
# shortseq is mono,di
# CAUTION: in case where the sequence is complementary, put a minus sign!!
def shortseqind(ca,se): #get nucleotide index in data
    # compute the short seq
    if ca=="ABC_i" or case=="ABC_old_i" or ca=="intra":
        if len(se) is 3:
            sse=se[1]
            sl=monos
        else:
            print "ERROR the seq must be a trinuc"
            return 1
    elif ca=="ABC_s" or case=="ABC_old_s" or ca=="step":
        if len(se) is 4:
            sse=se[1:3]
            sl=alldinucs
        else: 
            print "ERROR the seq must be a tetra"
            return 1
    elif ca=="NP":
        sse=se
        sl=alldinucs
    # test if in the shortseq list 
    try:
        return sl.index(sse)
    except:
        return -1*sl.index(compl_string(sse))

global reversemat #matrix for reversed sequences
reversemat=[-1,1,1,-1,1,1]

def seqlength(ca): #number of nucleotides
    if ca=="ABC_i" or case=="ABC_old_i" or ca=="intra":
        return 3
    if ca=="ABC_s" or case=="ABC_old_s" or ca=="step":
        return 4
    if ca=="NP":
        return 2

###### diverse functions

def nonetofloat(str):
    if str=="None":
        return float('nan')
    else:
        try:
            return float(str)
        except:
            return str
