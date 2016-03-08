#!/usr/bin/python
import numpy as np
import matplotlib as mpl
from scipy.optimize import bisect

mpl.use('pgf')

def figsize(scale):
    fig_width_pt = 550.0                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause scatters to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "text.fontsize": 10,
    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # scatters will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)

import matplotlib.pyplot as plt

# I make my own newfig and savefig functions
def newfig(width):
    plt.clf()
    fig = plt.figure(figsize=figsize(width))
    ax = fig.add_subplot(111)
    return fig, ax

def savefig(filename):
    plt.savefig('{}.pdf'.format(filename))


# Simple scatter
fig, ax  = newfig(0.6)

pure_nr_x19 = np.loadtxt("pure_nr_x19.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_x19 = np.loadtxt("moded_nr_x19.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_wls_x19 = np.loadtxt("moded_nr_wls_x19.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_bfgs_x19 = np.loadtxt("moded_nr_bfgs_x19.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_bfgs_wls_x19 = np.loadtxt("moded_nr_bfgs_wls_x19.txt",delimiter=",",skiprows=1,ndmin=2)

ax.scatter(pure_nr_x19[:,0], pure_nr_x19[:,1])
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Pure NR, x=.19')
savefig('pure_nr_x19')

ax.clear()
ax.scatter(moded_nr_x19[:,0], moded_nr_x19[:,1])
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified NR, x=.19')
savefig('moded_nr_x19')

ax.clear()
ax.scatter( moded_nr_wls_x19[:,0], moded_nr_wls_x19[:,1])
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified NR with Line Search, x=.19')
savefig('moded_nr_wls_x19')

ax.clear()
ax.scatter(moded_nr_bfgs_x19[:,0], moded_nr_bfgs_x19[:,1])
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified NR with BFGS, x=.19')
savefig('moded_nr_bfgs_x19')

ax.clear()
ax.scatter( moded_nr_bfgs_wls_x19[:,0], moded_nr_bfgs_wls_x19[:,1])
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified BR with BFGS and Line Search, x=.19')
savefig('moded_nr_bfgs_wls_x19')


pure_nr_x30 = np.loadtxt("pure_nr_x30.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_x30 = np.loadtxt("moded_nr_x30.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_wls_x30 = np.loadtxt("moded_nr_wls_x30.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_bfgs_x30 = np.loadtxt("moded_nr_bfgs_x30.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_bfgs_wls_x30 = np.loadtxt("moded_nr_bfgs_wls_x30.txt",delimiter=",",skiprows=1,ndmin=2)

ax.clear()
ax.scatter(pure_nr_x30[:,0], pure_nr_x30[:,1])
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Pure NR, x=.30')
savefig('pure_nr_x30')

ax.clear()
ax.scatter(moded_nr_x30[:,0], moded_nr_x30[:,1])
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified NR, x=.30')
savefig('moded_nr_x30')

ax.clear()
ax.scatter( moded_nr_wls_x30[:,0], moded_nr_wls_x30[:,1])
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified NR with Line Search, x=.30')
savefig('moded_nr_wls_x30')

ax.clear()
ax.scatter(moded_nr_bfgs_x30[:,0], moded_nr_bfgs_x30[:,1])
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified NR with BFGS, x=.30')
savefig('moded_nr_bfgs_x30')

ax.clear()
ax.scatter( moded_nr_bfgs_wls_x30[:,0], moded_nr_bfgs_wls_x30[:,1])
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified BR with BFGS and Line Search, x=.30')
savefig('moded_nr_bfgs_wls_x30')

ax.clear()
pure_nr_x19_conv = np.loadtxt("pure_nr_x19_conv.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_x19_conv = np.loadtxt("moded_nr_x19_conv.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_wls_x19_conv = np.loadtxt("moded_nr_wls_x19_conv.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_bfgs_x19_conv = np.loadtxt("moded_nr_bfgs_x19_conv.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_bfgs_wls_x19_conv = np.loadtxt("moded_nr_bfgs_wls_x19_conv.txt",delimiter=",",skiprows=1,ndmin=2)

ax.clear()
ax.scatter( range(1,len(pure_nr_x19_conv)+1), pure_nr_x19_conv)
ax.set_xlabel('Load Step')
ax.set_ylabel('Num Iterations')
ax.set_title('Pure NR, x=.19')
savefig('pure_nr_x19_conv')

ax.clear()
ax.scatter(range(1,len(moded_nr_x19_conv)+1) ,moded_nr_x19_conv)
ax.set_xlabel('Load Step')
ax.set_ylabel('Num Iterations')
ax.set_title('Modified NR, x=.19')
savefig('moded_nr_x19_conv')

ax.clear()
ax.scatter( range(1,len(moded_nr_wls_x19_conv)+1), moded_nr_wls_x19_conv)
ax.set_xlabel('Load Step')
ax.set_ylabel('Num Iterations')
ax.set_title('Modified NR with Line Search, x=.19')
savefig('moded_nr_wls_x19_conv')

ax.clear()
ax.scatter( range(1,len(moded_nr_bfgs_x19_conv)+1), moded_nr_bfgs_x19_conv)
ax.set_xlabel('Load Step')
ax.set_ylabel('Num Iterations')
ax.set_title('Modified NR with BFGS, x=.19')
savefig('moded_nr_bfgs_x19_conv')

ax.clear()
ax.scatter( range(1,len(moded_nr_bfgs_wls_x19_conv)+1), moded_nr_bfgs_wls_x19_conv)
ax.set_xlabel('Load Step')
ax.set_ylabel('Num Iterations')
ax.set_title('Modified BR with BFGS and Line Search, x=.19')
savefig('moded_nr_bfgs_wls_x19_conv')


pure_nr_x30_conv = np.loadtxt("pure_nr_x30_conv.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_x30_conv = np.loadtxt("moded_nr_x30_conv.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_wls_x30_conv = np.loadtxt("moded_nr_wls_x30_conv.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_bfgs_x30_conv = np.loadtxt("moded_nr_bfgs_x30_conv.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_bfgs_wls_x30_conv = np.loadtxt("moded_nr_bfgs_wls_x30_conv.txt",delimiter=",",skiprows=1,ndmin=2)

ax.clear()
ax.scatter( range(1,len(pure_nr_x30_conv)+1), pure_nr_x30_conv )
ax.set_xlabel('Load Step')
ax.set_ylabel('Num Iterations')
ax.set_title('Pure NR, x=.30')
savefig('pure_nr_x30_conv')

ax.clear()
ax.scatter( range(1,len(moded_nr_x30_conv)+1), moded_nr_x30_conv)
ax.set_xlabel('Load Step')
ax.set_ylabel('Num Iterations')
ax.set_title('Modified NR, x=.30')
savefig('moded_nr_x30_conv')

ax.clear()
ax.scatter( range(1,len(moded_nr_wls_x30_conv)+1), moded_nr_wls_x30_conv)
ax.set_xlabel('Load Step')
ax.set_ylabel('Num Iterations')
ax.set_title('Modified NR with Line Search, x=.30')
savefig('moded_nr_wls_x30_conv')

ax.clear()
ax.scatter( range(1,len(moded_nr_bfgs_x30_conv)+1), moded_nr_bfgs_x30_conv)
ax.set_xlabel('Load Step')
ax.set_ylabel('Num Iterations')
ax.set_title('Modified NR with BFGS, x=.30')
savefig('moded_nr_bfgs_x30_conv')

ax.clear()
ax.scatter( range(1,len(moded_nr_bfgs_wls_x30_conv)+1), moded_nr_bfgs_wls_x30_conv)
ax.set_xlabel('Load Step')
ax.set_ylabel('Num Iterations')
ax.set_title('Modified BR with BFGS and Line Search, x=.30')
savefig('moded_nr_bfgs_wls_x30_conv')

def N19_opt(myd1, myload):
    return myload - (0.19*myd1*myd1*myd1 -2.0*myd1*myd1 + 5.6*myd1)

def N30_opt(myd1, myload):
    return myload - (0.30*myd1*myd1*myd1 -2.0*myd1*myd1 + 5.6*myd1)

def N19_scatter(myd1):
    return 0.19*myd1*myd1*myd1 -2.0*myd1*myd1 +5.6*myd1

def N30_scatter(myd1):
    return 0.30*myd1*myd1*myd1 -2.0*myd1*myd1 +5.6*myd1

loads = np.linspace(.25, 40.0*.25, num=40)
n19_d = []
n30_d = []

for load in loads:
    n19_d.append( (bisect(N19_opt, -15.0,15.0, args=(load,) )) )
    n30_d.append( (bisect(N30_opt, -15.0,15.0, args=(load,) )) )

n19_d = np.array(n19_d)
n30_d = np.array(n30_d)

n19_n = N19_scatter(n19_d)
n30_n = N30_scatter(n30_d)

ax.clear()
ax.scatter(n19_d, n19_n)
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Exact N1 vs d1, x=.19')
savefig('exact_x19')

ax.clear()
ax.scatter(n30_d, n30_n)
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Exact N1 vs d1, x=.30')
savefig('exact_x30')
