from hoomd import *
from hoomd import md
import numpy as np

sigma = 3.741 #A
la =  6.66270
lr = 23.000
kBby10 = 0.83144622
eps = 353.55 * kBby10
mass = 44.01

#3500 beads
N = 3500

Lx = 46.824
Ly = Lx
Lz = 207.999


context.initialize("")

uc = lattice.fcc(a = 2.)

snap = uc.get_snapshot()
snap.replicate(5, 5, 35)
snap.particles.N
box = data.boxdim(Lx =  Lx, Ly =  Ly, Lz = Lz)
snap.box = box

snap.particles.types = ['CCG']
snap.particles.mass[:] = mass
snap.particles.diameter[:] = sigma

system = init.read_snapshot(snap)

nl = md.nlist.cell()

rcut = 6*sigma # A

mie = md.pair.mie(r_cut= rcut, nlist=nl)

mie.pair_coeff.set('CCG', 'CCG', epsilon= eps, sigma= sigma,
                   n = lr, m = la)

all = group.all()
standard = md.integrate.mode_standard(dt=0.003)
nve = md.integrate.nve(group = all, limit = 0.1)
run_upto(5e3, quiet = False)

Tquench = 600 #K
kTquench = Tquench * kBby10
nve.disable()

isothermal = md.integrate.nvt(group = all, kT = kTquench, tau = 1.0)
run_upto(1e5, quiet = False)

dump.gsd(filename="init.gsd", group=group.all(), period=None, time_step=0
         , dynamic=['attribute', 'momentum', 'topology'])