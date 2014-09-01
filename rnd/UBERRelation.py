#!/usr/bin/env python
from ase.structure import bulk
from kimcalculator import KIMCalculator
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab
import scipy.optimize as opt
import operator

def setPlotOptions(labelsize=20,tickmajor=20,tickminor=10,markersize=10,legendsize=20,legendspacing=1.5,labelsizexy=16):
    pylab.rcdefaults()
    pylab.rcParams.update({'xtick.labelsize':labelsizexy,\
            'xtick.major.size':tickmajor,\
            'xtick.minor.size':tickminor,\
            'ytick.labelsize':labelsizexy,\
            'ytick.major.size':tickmajor,\
            'ytick.minor.size':tickminor,\
            'lines.markersize':markersize,\
            'axes.labelsize':labelsize,\
            'legend.fontsize':legendsize,\
            'legend.columnspacing':legendspacing,\
            })

setPlotOptions()

def energy(a, slab, cell, positions):
    """ Compute the energy of a lattice given the lattice constant.
     Multiply the positions for every lattice type since in this
     test only diamond has multiple sites, all other positions are zero.
    """
    slab.set_positions( positions * a )
    slab.set_cell( cell * a )
    energy = slab.get_potential_energy()
    return energy / len(slab)

def run_get_energies(args):
    from multiprocessing import Process, Queue
    q = Queue()
    p = Process(target=get_energies, args=args+(q,))
    p.start()
    p.join()
    try:
        return q.get(block=False)
    except Exception as e:
        raise AttributeError("Calc failed for %r" % args)

def get_energies(model, symbol, lattice, q):
    calc = KIMCalculator(str(model))
    
    slab = bulk(str(symbol), str(lattice), a=1)
    #slab = slab.repeat((6,6,6))
    cell = slab.get_cell()
    positions = slab.get_positions()
    slab.set_calculator(calc)
    lattice_constant = 6

    aopt_arr, eopt, iterations, funccalls, warnflag = \
            opt.fmin(energy, lattice_constant, args=(slab,cell,positions),
                    full_output=True, disp=False, xtol=1e-12, maxfun=100)
    lattice_constant = aopt_arr
    cohesive_energy = -eopt/len(slab)
    dx = 1e-6*lattice_constant
    d2v = (energy(lattice_constant - dx, slab, cell, positions) -
            2*energy(lattice_constant, slab, cell, positions) +
            energy(lattice_constant + dx, slab, cell, positions)) / (dx*dx)

    l = np.sqrt( cohesive_energy / d2v )

    list_a = np.linspace(lattice_constant/2.0, lattice_constant*2., 101)
    list_E = np.array([energy(ta, slab, cell, positions) for ta in list_a])

    astar = (list_a - lattice_constant) / l
    Estar = list_E / cohesive_energy

    #return (astar, Estar)
    q.put((astar,Estar))

from pipeline import kimobjects
bads = [
    "EDIP_BOP_Bazant_Kaxiras_Si__MO_958932894036_001",
    "EMT_Asap_MetalGlass_CuMgZr__MO_655725647552_001",
    "EMT_Asap_Standard_Jacobsen_Stoltze_Norskov_AlAgAuCuNiPdPt__MO_118428466217_001",
    "Glue_Ercolessi_Adams_Al__MO_324507536345_001",
    "IMD_EAM_Schopf_AlNiCo_A__MO_122703700223_001",
    "MEAM_2NN_Fe_to_Ga__MO_145522277939_001",
    "kcc_meam_LiSi__MO_596436139350_001",
    "model_ArCHHeXe_BOP_AIREBO__MO_154399806462_001",
    "model_He_P_AbIn__MO_126942667206_001",
]
lattices = ['fcc', 'bcc', 'sc']
combos = []
for mo in kimobjects.Model.all():
    for l in lattices:
        for species in mo.kimspec['species']:
            if (mo.kim_code_version == '001' and str(mo) not in bads
                    and not str(mo).startswith("EMT")
                    and not str(mo).startswith("MEAM")):
                combos.append((mo, species, l)) 
                #print mo, l, species

alist = []
Elist = []
for i in xrange(0,len(combos)):
    print i, combos[i][0]
    try:
        astar, Estar = run_get_energies(combos[i])
        alist.append(astar)
        Elist.append(Estar)
    except:
        print "failed."

fig = pylab.figure()
for i, (a, E) in enumerate(zip(alist, Elist)):
    if i == 0:
        pylab.plot(a, E, 'k-', linewidth=1.3, alpha=0.08, label='Model predictions')  # use 0.05 if using all lines
    else:
        pylab.plot(a, E, 'k-', linewidth=1.3, alpha=0.08)  # use 0.05 if using all lines

xs = np.linspace(-1.3, 5, 10000)
ys = -(1+xs)*np.exp(-xs)
pylab.plot(xs, ys, 'b-', linewidth=2.5, label='Rydberg function')

pylab.xlabel(r"$\rm{Scaled\,separation}\,\,a^{\star}$")
pylab.ylabel(r"$\rm{Scaled\,energy}\,\,E^{\star}$")
pylab.xlim(-1.3, 3)
pylab.ylim(-1, 1)
pylab.legend(loc='upper right')
pylab.show()
pylab.savefig("uber.png")
