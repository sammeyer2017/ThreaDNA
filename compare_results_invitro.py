import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

plt.rcParams.update({'pdf.fonttype': 42})
plt.rcParams.update({'ps.fonttype': 42})
plt.rcParams.update({'font.size': 11})
plt.rcParams.update({'legend.fontsize': 10})
plt.rcParams.update({'font.family': "Arial"})


kT_to_kcalmol=1/0.593



def extract_indiv_values_from_bed(bedfile, outfile=None, position=None):
    """
    Takes a bedfile with profiles of several sequences, and computes the energy for a single (actual) position. 
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



def compare_predicted_to_measured_affinities(exp_file, outfile, indirect_energy_file=None, direct_energy_file=None, indir_factor=1., plot=False):
    """
    takes in one or two energy files, and computes the comparison of energies with exp
    outputs the results in outfile, and optionally as a pdf graph
    input exp file: name and foldchange
    the exp and computed files must have the same sequence names in the first column!
    """
    expfc=pd.read_table(exp_file,header=None,sep='\t')
    # compute energies in kT from expfc
    e=-np.log(expfc[1])
    print e
    t=np.zeros(len(e))
    if indirect_energy_file!=None:
        indir=pd.read_table(indirect_energy_file,header=None,sep='\t')
        t+=indir_factor*np.array(indir[1])
    if direct_energy_file!=None:
        direct=pd.read_table(direct_energy_file,header=None,sep='\t')
        t+=np.array(direct[1]-direct[1][0])
        print t
    if max(t)==0.:
        print "Must provide at least one predicted file !!"
        return 1
    # ---------------------
    a, _, _, _ = np.linalg.lstsq(np.array([e]).T, t)
    r2=pearsonr(e, t)
    print "Slope of the regression: %.2f"%a
    print "Pearson's r^2 and p-val: %.3f, %.2e"%r2
    o=open(outfile, 'w')
    o.write("Experimental file: %s\n"%exp_file)
    o.write("Direct energy file: %s\n"%direct_energy_file)
    o.write("Indirect energy file: %s\n"%indirect_energy_file)
    o.write("Energy factor for indirect energy: %.1f\n"%indir_factor)
    o.write("---------------------")
    o.write("Pearson's r^2 and p-val: %.3f, %.2e\n"%r2)
    o.write("Slope of the regression: %.2f\n"%a)
    o.close()
    # ---------------------
    
    if plot:
        plt.figure(figsize=(4,3))
        plt.plot(e/kT_to_kcalmol,t/kT_to_kcalmol/a,ls="",marker="o",color="black",label="$r^2=%.2f$"%r2[0])
        plt.xlabel(r"$\Delta G_{exp}$ (kcal/mol)")
        plt.ylabel(r"$\Delta G_{pred}$ (kcal/mol)")
        plt.axhline(0,color="black")
        plt.axvline(0,color="black")
        plt.plot([-2,2],[-2,2],color="black",ls="--")
        plt.legend(loc=2,)
        plt.tight_layout()
        plt.savefig(outfile.split(".")[0]+".pdf")
    return 0
    
    
    
# simple pwm
extract_indiv_values_from_bed("lindemose_from_rdbpwm.bed", position=8)
compare_predicted_to_measured_affinities("crp_lindemose_edited.csv", outfile="comparison_simple_pwm.dat", indirect_energy_file=None, direct_energy_file="lindemose_from_rdbpwm.dat", plot=True)
