from hoomd import *
from hoomd import md
import numpy as np

sigma = 4.629 #A
la =  6.00
lr = 19.61
kBby10 = 0.83144622
eps = 414.90 * kBby10
ms = 3
mass = 142.2817 / ms

#equals to 1668 molecules of 3 beads
N = 1668

Lx = 57.63
Ly = Lx
Lz = 345.78


context.initialize("")

uc = lattice.fcc(a = 2.)

snap = uc.get_snapshot()
snap.replicate(3, 3, 139)
snap.particles.N
box = data.boxdim(Lx =  Lx, Ly =  Ly, Lz = Lz)
snap.box = box

snap.particles.types = ['DCG']
snap.particles.mass[:] = mass
snap.particles.diameter[:] = sigma

nbonds = int(N * ms * 2 / 3)
snap.bonds.resize(nbonds) 
snap.bonds.group[:] = [[i, i+1] for j in range(0, N*ms, ms) for i in [j, j+1]]
snap.bonds.types = ['bndDCG']

nangles = int(N * ms * 1 / 3)
snap.angles.resize(nangles)
snap.angles.group[:] = [[i, i+1, i+2] for j in range(0, N*ms, ms) for i in [j]]
snap.angles.types = ['angDCG']

system = init.read_snapshot(snap)

nl = md.nlist.cell()

rcut = 6*sigma # A

mie = md.pair.mie(r_cut= rcut, nlist=nl)

mie.pair_coeff.set('DCG', 'DCG', epsilon= eps, sigma= sigma,
                   n = lr, m = la)

#Angles and Bond obtained from raaSAFT
kcal2joule=4184
AlkaneBondConstant = 2*7.540*kcal2joule*0.1
AlkaneAngleZero = 157.6*np.pi/180.0
AlkaneAngleConstant = 2*2.650*kcal2joule*0.1

bond = md.bond.harmonic() 
bond.bond_coeff.set('bndDCG', k=AlkaneBondConstant , r0=sigma)

angle = md.angle.harmonic()
angle.angle_coeff.set('angDCG', k = AlkaneAngleConstant, t0 = AlkaneAngleZero)


all = group.all()
standard = md.integrate.mode_standard(dt=0.002)
nve = md.integrate.nve(group = all, limit = 0.1)
run_upto(5e3, quiet = False)

Tquench = 1000 #K
kTquench = Tquench * kBby10
nve.disable()

isothermal = md.integrate.nvt(group = all, kT = kTquench, tau = 1.0)
run(1e5, quiet = False)

dump.gsd(filename="init.gsd", group=group.all(), period=None, time_step=0
         , dynamic=['attribute', 'momentum', 'topology'])