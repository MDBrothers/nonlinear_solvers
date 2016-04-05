#!/usr/bin/python
import numpy as np
import matplotlib as mpl
from tabulate import tabulate


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

def N_explicit(myd1, x):
    return x*np.power(myd1,3.0) -2.0*np.power(myd1, 2.0) +5.6*myd1

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


header = ["pure d", "pure N","modified d", "modified N","modified wls d", "modified wls N","BFGS d", "BFGS N","BFGS wls d", "BFGS wls N"]
content = np.array([pure_nr_x19[:,0], pure_nr_x19[:,1], moded_nr_x19[:,0], moded_nr_x19[:,1], moded_nr_wls_x19[:,0], moded_nr_wls_x19[:,1], moded_nr_bfgs_x19[:,0], moded_nr_bfgs_x19[:,1], moded_nr_bfgs_wls_x19[:,0], moded_nr_bfgs_wls_x19[:,1]], ndmin=2)
print "\\begin{figure}"
print tabulate(content.transpose(), header, tablefmt="latex")
print "\\caption{Solution approximations, $x=.19$}"
print "\\end{figure}"

ax.scatter(pure_nr_x19[:,0], pure_nr_x19[:,1])
ax.plot(np.linspace(0,8.0, 100), N_explicit(np.linspace(0,8.0,100), .19))
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Pure NR, x=.19')
savefig('pure_nr_x19')

ax.clear()
ax.scatter(moded_nr_x19[:,0], moded_nr_x19[:,1])
ax.plot(np.linspace(0,8.0, 100), N_explicit(np.linspace(0,8.0,100), .19))
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified NR, x=.19')
savefig('moded_nr_x19')

ax.clear()
ax.scatter( moded_nr_wls_x19[:,0], moded_nr_wls_x19[:,1])
ax.plot(np.linspace(0,8.0, 100), N_explicit(np.linspace(0,8.0,100), .19))
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified NR with Line Search, x=.19')
savefig('moded_nr_wls_x19')

ax.clear()
ax.scatter(moded_nr_bfgs_x19[:,0], moded_nr_bfgs_x19[:,1])
ax.plot(np.linspace(0,8.0, 100), N_explicit(np.linspace(0,8.0,100), .19))
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified NR with BFGS, x=.19')
savefig('moded_nr_bfgs_x19')

ax.clear()
ax.scatter( moded_nr_bfgs_wls_x19[:,0], moded_nr_bfgs_wls_x19[:,1])
ax.plot(np.linspace(0,8.0, 100), N_explicit(np.linspace(0,8.0,100), .19))
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified BR with BFGS and Line Search, x=.19')
savefig('moded_nr_bfgs_wls_x19')


pure_nr_x30 = np.loadtxt("pure_nr_x30.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_x30 = np.loadtxt("moded_nr_x30.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_wls_x30 = np.loadtxt("moded_nr_wls_x30.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_bfgs_x30 = np.loadtxt("moded_nr_bfgs_x30.txt",delimiter=",",skiprows=1,ndmin=2)
moded_nr_bfgs_wls_x30 = np.loadtxt("moded_nr_bfgs_wls_x30.txt",delimiter=",",skiprows=1,ndmin=2)

header = ["pure d", "pure N","modified d", "modified N","modified wls d", "modified wls N","BFGS d", "BFGS N","BFGS wls d", "BFGS wls N"]
content = np.array([pure_nr_x30[:,0], pure_nr_x30[:,1], moded_nr_x30[:,0], moded_nr_x30[:,1], moded_nr_wls_x30[:,0], moded_nr_wls_x30[:,1], moded_nr_bfgs_x30[:,0], moded_nr_bfgs_x30[:,1], moded_nr_bfgs_wls_x30[:,0], moded_nr_bfgs_wls_x30[:,1]], ndmin=2)
print "\\begin{figure}"
print tabulate(content.transpose(), header, tablefmt="latex")
print "\\caption{Solution approximations, $x=.30$}"
print "\\end{figure}"

ax.clear()
ax.scatter(pure_nr_x30[:,0], pure_nr_x30[:,1])
ax.plot(np.linspace(0,4.5, 100), N_explicit(np.linspace(0,4.5,100), .30))
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Pure NR, x=.30')
savefig('pure_nr_x30')

ax.clear()
ax.scatter(moded_nr_x30[:,0], moded_nr_x30[:,1])
ax.plot(np.linspace(0,4.5, 100), N_explicit(np.linspace(0,4.5,100), .30))
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified NR, x=.30')
savefig('moded_nr_x30')

ax.clear()
ax.scatter( moded_nr_wls_x30[:,0], moded_nr_wls_x30[:,1])
ax.plot(np.linspace(0,4.5, 100), N_explicit(np.linspace(0,4.5,100), .30))
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified NR with Line Search, x=.30')
savefig('moded_nr_wls_x30')

ax.clear()
ax.scatter(moded_nr_bfgs_x30[:,0], moded_nr_bfgs_x30[:,1])
ax.plot(np.linspace(0,4.5, 100), N_explicit(np.linspace(0,4.5,100), .30))
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified NR with BFGS, x=.30')
savefig('moded_nr_bfgs_x30')

ax.clear()
ax.scatter( moded_nr_bfgs_wls_x30[:,0], moded_nr_bfgs_wls_x30[:,1])
ax.plot(np.linspace(0,4.5, 100), N_explicit(np.linspace(0,4.5,100), .30))
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Modified BR with BFGS and Line Search, x=.30')
savefig('moded_nr_bfgs_wls_x30')

ax.clear()
pure_nr_x19_conv = np.loadtxt("pure_nr_x19_conv.txt",delimiter=",",skiprows=1,ndmin=1)
moded_nr_x19_conv = np.loadtxt("moded_nr_x19_conv.txt",delimiter=",",skiprows=1,ndmin=1)
moded_nr_wls_x19_conv = np.loadtxt("moded_nr_wls_x19_conv.txt",delimiter=",",skiprows=1,ndmin=1)
moded_nr_bfgs_x19_conv = np.loadtxt("moded_nr_bfgs_x19_conv.txt",delimiter=",",skiprows=1,ndmin=1)
moded_nr_bfgs_wls_x19_conv = np.loadtxt("moded_nr_bfgs_wls_x19_conv.txt",delimiter=",",skiprows=1,ndmin=1)

header = ["pure","modified","modified wls","BFGS","BFGS wls"]
content = np.ndarray(shape =(40,5))
content[:,0] = pure_nr_x19_conv
content[:,1] = moded_nr_x19_conv
content[:,2] = moded_nr_wls_x19_conv
content[:,3] = moded_nr_bfgs_x19_conv
content[:,4] = moded_nr_bfgs_wls_x19_conv

print "\\begin{figure}"
print tabulate(content, header, tablefmt="latex")
print "\\caption{Number of iterations, $x=.19$}"
print "\\end{figure}"

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


pure_nr_x30_conv = np.loadtxt("pure_nr_x30_conv.txt",delimiter=",",skiprows=1,ndmin=1)
moded_nr_x30_conv = np.loadtxt("moded_nr_x30_conv.txt",delimiter=",",skiprows=1,ndmin=1)
moded_nr_wls_x30_conv = np.loadtxt("moded_nr_wls_x30_conv.txt",delimiter=",",skiprows=1,ndmin=1)
moded_nr_bfgs_x30_conv = np.loadtxt("moded_nr_bfgs_x30_conv.txt",delimiter=",",skiprows=1,ndmin=1)
moded_nr_bfgs_wls_x30_conv = np.loadtxt("moded_nr_bfgs_wls_x30_conv.txt",delimiter=",",skiprows=1,ndmin=1)

header = ["pure","modified","modified wls","BFGS","BFGS wls"]
content = np.ndarray(shape =(40,5))
content[:,0] = pure_nr_x30_conv
content[:,1] = moded_nr_x30_conv
content[:,2] = moded_nr_wls_x30_conv
content[:,3] = moded_nr_bfgs_x30_conv
content[:,4] = moded_nr_bfgs_wls_x30_conv

print "\\begin{figure}"
print tabulate(content, header, tablefmt="latex")
print "\\caption{Number of iterations, $x=.30$}"
print "\\end{figure}"

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

ax.clear()
ax.plot(np.linspace(0,8.0, 100), N_explicit(np.linspace(0,8.0,100), .19))
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Exact N1 vs d1, x=.19')
savefig('exact_x19')

ax.clear()
ax.plot(np.linspace(0,4.5, 100), N_explicit(np.linspace(0,4.5,100), .30))
ax.set_xlabel('d1')
ax.set_ylabel('N1')
ax.set_title('Exact N1 vs d1, x=.30')
savefig('exact_x30')
