#!/usr/bin/env python
# coding: utf8

from pwm import *
import matplotlib.pyplot as plt
import matplotlib.patches as pat

plt.rcParams.update({'pdf.fonttype': 42})
plt.rcParams.update({'ps.fonttype': 42})
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'legend.fontsize': 10})
#plt.rcParams.update({'mathtex.fontset': "cm"})
plt.rcParams.update({'font.family': "Arial"})


basecolors=dict([("A","green"), ("C", "blue"), ("T", "red"), ("G", "yellow")])


def get_colors(ind):
    revind=dict([(ind[b],b) for b in ind.keys()])
    n=len(ind.keys())
    if n==16:
        # dinuc
        return [(basecolors[revind[i][0]], basecolors[revind[i][1]]) for i in range(n)]
    elif n==4:
        # mononuc
        return [(basecolors[revind[i][0]], basecolors[revind[i][0]]) for i in range(n)]



def plot_pwm(proba_pwmfile, outfile=None, figsize=6, resol=300):
    """
    plots a PWM probability file using the standard representation
    """
    if outfile==None:
        outfile=proba_pwmfile.split(".")[0]+".pdf"
    ind,mat=load_pwm(proba_pwmfile)
    colors=get_colors(ind)
    fig=plt.figure(figsize=(figsize,2))
    ax=plt.subplot()
    plt.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.32)
    # height of each number
    n=len(ind.keys())
    le=len(mat)
    totheights=np.sum(np.log2((mat+10**-8)*n)*mat, axis=1, keepdims=True)
    heights=mat*totheights
    if n==16:
        # dinuc:
        mononucs=np.arange(1,le+2, 1)
        plt.xlim(.5,le+1.5)
        plt.xticks(range(1,le+2))
        xs=np.arange(1.5, le+.5, 1)
    elif n==4:
        # mononuc:
        mononucs=np.arange(1,le+1, 1)
        plt.xlim(0.5,le+.5)
        plt.xticks(range(1,le+1))
        xs=np.arange(1, le+1, 1)
    t=0
    for ix, x in enumerate(xs):
        hes=heights[ix]
        top=np.cumsum(hes)
        if max(top)>t:
            t=max(top)
        for ic,c in enumerate(colors):
            if hes[ic]!=0.:
                ax.add_patch(pat.Rectangle((x-.3,top[ic]-hes[ic]),.3,hes[ic],color=c[0]))
                ax.add_patch(pat.Rectangle((x,top[ic]-hes[ic]),.3,hes[ic],color=c[1]))
    plt.axvline(np.mean(mononucs), ls="--", color="black", lw=.8)
    if t<2:
        plt.ylim(0,2)
    else:
        plt.ylim(0,t*1.1)
    plt.ylabel("bits")
    plt.xlabel("basepair along protein")
    plt.savefig(outfile, dpi=resol)
    plt.close()

    
def plot_pwm_from_sequences(aligned_sequence_file, outfile=None, probafile=None):
    if outfile==None:
        outfile=aligned_sequence_file.split(".")[0]+".pdf"
    if probafile==None:
        probafile=aligned_sequence_file.split(".")[0]+".pwm"
    pwm=compute_pwm_from_sequences(aligned_sequence_file, outfile=probafile)
    plot_pwm(probafile, outfile=outfile)
    return 0




# -------------- EXECUTION
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot a mono/dinucleotide position weight-matrix in the style of ThreaDNA. For a dinucleotide PWM, it is recommended to plot the probability rather than the energy PWM. Output in pdf format. Execution requires the MatPlotLib library.')
    # list of operations
    # classical pwm from aligned sequences
    # proba/energy pwm to energy profile along sequence
    parser.add_argument('input', type=str,help='Input PWM file in JASPAR format')
    parser.add_argument("-o","--output",type=str,action="store",help="Output pdf file")
    args=parser.parse_args()
    plot_pwm(unicode(args.input),unicode(args.output)) #executes program
    
