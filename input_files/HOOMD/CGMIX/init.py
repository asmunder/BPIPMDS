from hoomd import *
from hoomd import md
import numpy as np

kBby10 = 0.83144622

sigmaCO2 = 3.741 #A
laCO2 =  6.66270
lrCO2 = 23.000
epsCO2 = 353.55 * kBby10
msCO2 = 1
NCO2 = 2256

sigmaC10 = 4.629 #A
laC10 =  6.0
lrC10 = 19.61
epsC10 = 414.90 * kBby10
msC10 = 3
NC10 = 1248

Ntotal = NCO2 + msC10 * NC10

kij = 0.075
sigmaij = (sigmaCO2 + sigmaC10)/2
laij = ((laC10 - 3) * (laCO2 - 3))**0.5 + 3
lrij = ((lrC10 - 3) * (lrCO2 - 3))**0.5 + 3
epsij = (epsC10 * epsCO2)**0.5
epsij *= (sigmaC10**3 * sigmaCO2**3)**0.5 / sigmaij**3
epsij *= (1 - kij)

Lx = 76.5
Ly = Lx
Lz = 459.0
box = data.boxdim(Lx =  Lx, Ly =  Ly, Lz = Lz)

context.initialize("")

uc = lattice.fcc(a = 2.)

snap = uc.get_snapshot()
snap.replicate(10, 10, 15)
snap.particles.N
snap.box = box

snap.particles.types = ['CCG', 'DCG']

snap.particles.mass[:NCO2] = 44.0100000
snap.particles.diameter[:NCO2] = sigmaCO2
snap.particles.typeid[:NCO2] = 0

snap.particles.mass[NCO2:] = 142.281700 / msC10
snap.particles.diameter[NCO2:] = sigmaC10 
snap.particles.typeid[NCO2:] = 1

nbonds = int(NC10 * msC10 * 2 / 3)
snap.bonds.resize(nbonds) 
snap.bonds.group[:] = [[i, i+1] for j in range(NCO2, NCO2 + NC10*msC10, msC10) for i in [j, j+1]]
snap.bonds.types = ['bndDCG']

nangles = int(NC10 * msC10 * 1 / 3)
snap.angles.resize(nangles)
snap.angles.group[:] = [[i, i+1, i+2] for j in range(NCO2, NCO2 + NC10*msC10, msC10) for i in [j]]
snap.angles.types = ['angDCG']

system = init.read_snapshot(snap)

nl = md.nlist.cell()

rcut = 28.00 # A

mie = md.pair.mie(r_cut= rcut, nlist=nl)

mie.pair_coeff.set('CCG', 'CCG', epsilon= epsCO2, sigma= sigmaCO2,
                   n = lrCO2, m = laCO2)

mie.pair_coeff.set('DCG', 'DCG', epsilon= epsC10, sigma= sigmaC10,
                   n = lrC10, m = laC10)

mie.pair_coeff.set('CCG', 'DCG', epsilon= epsij, sigma= sigmaij,
                   n = lrij, m = laij)

kcal2joule=4184
AlkaneBondConstant = 2*7.540*kcal2joule*0.1
AlkaneAngleZero = 157.6*np.pi/180.0
AlkaneAngleConstant = 2*2.650*kcal2joule*0.1

bond = md.bond.harmonic() 
bond.bond_coeff.set('bndDCG', k=AlkaneBondConstant , r0=sigmaC10)

angle = md.angle.harmonic()
angle.angle_coeff.set('angDCG', k = AlkaneAngleConstant, t0 = AlkaneAngleZero)

all = group.all()
standard = md.integrate.mode_standard(dt=0.003)
nve = md.integrate.nve(group = all, limit = 0.1)
run_upto(5e3, quiet = False)

Tquench = 1000 #K
kTquench = Tquench * kBby10
nve.disable()

isothermal = md.integrate.nvt(group = all, kT = kTquench, tau = 1.0)
run_upto(1e5, quiet = False)

dump.gsd(filename="init.gsd", group=group.all(), period=None, time_step=0
         , dynamic=['attribute', 'momentum', 'topology'])
