from hoomd import *
from hoomd import md
import numpy as np


sigma = 3.741 #A
la =  6.66270
lr = 23.000
kBby10 = 0.83144622
eps = 353.55 * kBby10
ms = 1

N = 3500

Lx = 46.824
Ly = Lx
Lz = 207.999

Nequilibrium = 5000000
Nproduction = 12000000
Ntotal  = Nequilibrium + Nproduction

context.initialize("")

system = init.read_gsd(filename='init.gsd', restart='restart.gsd')


nl = md.nlist.cell()
rcut = 6*sigma # A

#Setting Mie interaction parameters
mie = md.pair.mie(r_cut= rcut, nlist=nl)
mie.pair_coeff.set('CCG', 'CCG', epsilon= eps, sigma= sigma,
                   n = lr, m = la)

all = group.all()
standard = md.integrate.mode_standard(dt = 0.003)

#write restart each 100k steps
restart = dump.gsd(filename="restart.gsd", group = all, truncate=True, 
		period = 100000, phase=0)

T = 240.0 #K
kT = T * kBby10
isothermal = md.integrate.nvt(group = all, kT = kT, tau = 1.0)
run_upto(Nequilibrium, quiet = False)
isothermal.disable()


period = 1000
#data to be logged
loglist = ['potential_energy', 'kinetic_energy',
                     'temperature', 'pressure_zz', 'pressure_xx', 'pressure_yy']
log = analyze.log(filename = 'log.dat', quantities = loglist, period = period, 
		overwrite=False, phase = -1)


nve = md.integrate.nve(group = all)
run_upto(Ntotal)