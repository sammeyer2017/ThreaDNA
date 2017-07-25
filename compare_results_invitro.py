import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

plt.rcParams.update({'pdf.fonttype': 42})
plt.rcParams.update({'ps.fonttype': 42})
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'legend.fontsize': 10})
plt.rcParams.update({'font.family': "Arial"})

#import statsmodels.api as sm

kT_to_kcalmol=1/0.593

def extract_indiv_values_from_bed(bedfile, outfile=None, position=None):
    """
    Takes a bedfile with profiles of several sequences, and computes the energy for a single (actual) position. The code checks that the provided position is the one with minimal energy. 
    """
    # assumes all sequences have same length!!!
    if outfile==None:
        outfile=bedfile.split(".")[0]+".dat"
    bedres=pd.read_table(bedfile,header=None,skiprows=1)
    names=bedres[0].unique()
    nl=len(names)
    rl=len(bedres)/nl
    mid=rl/2
    try: 
        v=np.reshape(np.array([bedres[3]]), (nl, rl))
    except:
        print "it seems that all sequences do not have the same size !! here we need this"
        return 1
    mins=np.argmin(v, axis=1)
    mi=mins[0]
    if not np.all(np.equal(mins, mi*np.ones(len(mins)))):
        print "it seems that some sequences have not their minimal energy at the same position; you should check this !"
        print "indexes of minimal energy: %s"%str(mins)
    if position==None:
        vs=v[:,mi]
    elif isinstance(position, int):
        if position==mi:
            print "the given position is indeed the minimum of energy"
        else:
            print "CAUTION: the given position is NOT the minimum of energy of the first sequence, this is embarrassing: please check your position, otherwise you may get absurd results"
            print "indexes of minimal energy: %s"%str(min)
        vs=v[:,position]
    d=pd.DataFrame({"name": names, "energy": vs})
    d.to_csv(outfile, sep="\t", columns=("name", "energy"), header=False, index=False)



def compare_predicted_to_measured_affinities(exp_file, outfile, indirect_energy_file=None, direct_energy_file=None, indir_factor=1., dir_factor=1., plot=False):
    """
    takes in one or two energy files, and computes the comparison of energies with exp
    outputs the results in outfile, and optionally as a pdf graph
    input exp file: name and foldchange
    # CAUTION: the provided values are FOLDCHANGES, i.e. proportional to 1/affinity in concentration !!! the larger, the most affine here.  
    the exp and computed files must have the same sequence names in the first column!
    """
    expfc=pd.read_table(exp_file,header=None,sep='\t')
    # compute energies in kT from expfc
    e=-np.log(expfc[1])
    # print e
    t=np.zeros(len(e))
    if indirect_energy_file!=None:
        indirect=pd.read_table(indirect_energy_file,header=None,sep='\t')
        indir=indir_factor*np.array(indirect[1]-indirect[1][0])
        t+=indir
    if direct_energy_file!=None:
        direct=pd.read_table(direct_energy_file,header=None,sep='\t')
        dire=dir_factor*np.array(direct[1]-direct[1][0])
        t+=dire
    # print t
    if max(t)==0.:
        print "Must provide at least one predicted file !!"
        return 1
    # ---------------------
    if direct_energy_file!=None:
        ad, _, _, _ = np.linalg.lstsq(np.array([e]).T, dire)
        r2d=pearsonr(e, dire)
    if indirect_energy_file!=None:
        ai, _, _, _ = np.linalg.lstsq(np.array([e]).T, indir)
        r2i=pearsonr(e, indir)
    ac, _, _, _ = np.linalg.lstsq(np.array([e]).T, t)
    #rlm_model = sm.RLM(e, t, M=sm.robust.norms.HuberT()).fit()
    #a=rlm_model.params[0]
    r2c=pearsonr(e, t)
    print "Slope of the regression: %.2f"%ac
    print "Pearson's r^2 and p-val: %.3f, %.2e"%(r2c[0]**2, r2c[1])
    o=open(outfile, 'w')
    o.write("Experimental file: %s\n"%exp_file)
    o.write("Direct energy file: %s\n"%direct_energy_file)
    o.write("Indirect energy file: %s\n"%indirect_energy_file)
    o.write("Energy factor for direct energy: %.1f\n"%dir_factor)
    o.write("Energy factor for indirect energy: %.1f\n"%indir_factor)
    o.write("---------------------\n")
    o.write("Pearson's corr r^2 and p-val: %.3f, %.2e\n"%(r2c[0]**2, r2c[1]))
    o.write("Slope of the regression: %.2f\n"%ac)
    o.close()
    # ---------------------
    
    if plot:
        plt.figure(figsize=(4.5,3))
        plt.plot([-2,2],[-2,2],color="black",ls="-", lw=.8)
        if direct_energy_file!=None:
            plt.plot(e/kT_to_kcalmol,dire/kT_to_kcalmol/ad,ls="",marker="<",ms=3, color="gray",label="PWM: $r^2=%.2f$"%r2d[0]**2)
        if indirect_energy_file!=None:
            plt.plot(e/kT_to_kcalmol,indir/kT_to_kcalmol/ai,ls="",marker=">",ms=3, color="black",label="ThreaDNA: $%.2f$"%r2i[0]**2)
        plt.plot(e/kT_to_kcalmol,t/kT_to_kcalmol/ac,ls="",marker="o",ms=4, color="blue",label="combination: $%.2f$"%(r2c[0]**2))
        plt.xlabel(r"$\Delta G_{exp}$ (kcal/mol)")
        plt.ylabel(r"$\Delta G_{pred}$ (kcal/mol)")
        plt.plot([-2,5],[-2,5],color="black",ls="-", lw=.8)
        #plt.ylim(-3,3)
        #plt.xlim(-2,2)
        plt.axhline(0,color="black")
        plt.axvline(0,color="black")
        plt.legend(loc=2,)
        plt.tight_layout()
        plt.savefig(outfile.split(".")[0]+".pdf")
    return 0
    

def plot_correlations(values, files, outfile, minmax):
    corr=[]
    for iv, v in enumerate(values):
        with open(files[iv],"r") as f:
            corr.append(float(f.readlines()[-2].split(":")[1].split(",")[0]))
    plt.figure(figsize=(4,2.5))
    plt.subplots_adjust(left=0.15, right=0.86, top=0.95, bottom=0.18)
    ax1=plt.subplot()
    plt.plot(values, corr, color="blue")
    maxcorr=max(corr)
    mci=corr.index(maxcorr)
    plt.axvline(values[mci], ls="--", color="blue", lw=1)
    plt.axhline(maxcorr, ls="--", color="blue", lw=1)
    plt.text(.67, .83, "combination", color="blue")
    plt.xlabel("Fraction of direct contribution")
    plt.ylabel("Pearson's correlation $r^2$")
    plt.axhline(corr[0], color="black", ls="--", lw=.5)
    plt.axhline(corr[-1], color="gray", ls="--", lw=.5)
    plt.yticks(np.arange(0,1,.1))
    plt.ylim(minmax)
    plt.xlim(0,1)
    plt.text(0.01, .48, "ThreaDNA", color="black")
    plt.text(.85, .25, "PWM", color="gray")
    ax2 = ax1.twinx()
    mi,ma=minmax
    fact=1/(ma-mi)
    pvalslog=-np.arange(2,12)
    corrvals=[.496, .607, .689, .752, .7985, .837, .8671, .89145, .9112, .92726]#, .95425]
    pvalsloglab=[str(x) for x in pvalslog]
    pvalsloglab[-2]=""
    pvalsloglab[-4]=""
    ax2.set_yticks([(x*x-mi)*fact for x in corrvals])
    ax2.set_yticklabels(pvalsloglab,fontsize=11)
    ax2.set_ylabel("correlation p-value ($\log_{10}$)")
    #plt.tight_layout()
    plt.savefig(outfile)
    plt.close()
    return values[mci]

 


"""
# simple pwm
#extract_indiv_values_from_bed("lindemose_from_rdbpwm.bed", position=8)
compare_predicted_to_measured_affinities("crp_lindemose_edited.csv", outfile="comparison_simple_pwm.dat", indirect_energy_file=None, direct_energy_file="lindemose_from_rdbpwm.dat", dir_factor=2.55, plot=True)

# pure threadna
#extract_indiv_values_from_bed("crp_lindemose_CRP.bed", position=6)
compare_predicted_to_measured_affinities("crp_lindemose_edited.csv", outfile="comparison_pure_threadna.dat", indirect_energy_file="crp_lindemose_CRP.dat", indir_factor=50., plot=True)


# combination of threadna and original pwm
for dirfrac in np.arange(0, 1.01, .05):
    indirfrac=1-dirfrac
    compare_predicted_to_measured_affinities("crp_lindemose_edited.csv", outfile="comparison_threadna_+_%d_original_pwm.dat"%(100*dirfrac), indirect_energy_file="crp_lindemose_CRP.dat", indir_factor=indirfrac*50., direct_energy_file="lindemose_from_rdbpwm.dat", dir_factor=dirfrac*2.55, plot=False)
"""

"""
bestfrac=plot_correlations(np.arange(0, 1.01, .05), ["comparison_threadna_+_%d_original_pwm.dat"%(100*dirfrac) for dirfrac in np.arange(0, 1.01, .05)], "comparisons_threadna_purePWM.pdf", (.23, .9))
print bestfrac
compare_predicted_to_measured_affinities("crp_lindemose_edited.csv", outfile="comparison_threadna_+_%d_original_pwm.dat"%(100*bestfrac), indirect_energy_file="crp_lindemose_CRP.dat", indir_factor=(1-bestfrac)*50., direct_energy_file="lindemose_from_rdbpwm.dat", dir_factor=bestfrac*2.55, plot=True)


# pure threadna transformed into monoPWM: show that it does not work !!
extract_indiv_values_from_bed("lindemose_from_threadna_monopwm.bed", position=6)
compare_predicted_to_measured_affinities("crp_lindemose_edited.csv", outfile="comparison_pure_threadna_monopwm.dat", indirect_energy_file="lindemose_from_threadna_monopwm.dat", indir_factor=50., plot=True)

"""
