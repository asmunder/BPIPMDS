from hoomd import *
from hoomd import md
import numpy as np

kBby10 = 0.83144622

#CO2 Parameters
sigmaCO2 = 3.741 #A
laCO2 =  6.66270
lrCO2 = 23.000
epsCO2 = 353.55 * kBby10
msCO2 = 1
NCO2 = 2256


#C10 Parameters
sigmaC10 = 4.629 #A
laC10 =  6.0
lrC10 = 19.61
epsC10 = 414.90 * kBby10
msC10 = 3
NC10 = 1248


#Total number of beads
Ntotal = NCO2 + msC10 * NC10

#Cross interaction parameters
kij = 0.075
sigmaij = (sigmaCO2 + sigmaC10)/2
laij = ((laC10 - 3) * (laCO2 - 3))**0.5 + 3
lrij = ((lrC10 - 3) * (lrCO2 - 3))**0.5 + 3
epsij = (epsC10 * epsCO2)**0.5
epsij *= (sigmaC10**3 * sigmaCO2**3)**0.5 / sigmaij**3
epsij *= (1 - kij)


Lx = 76.5
Ly = 76.5
Lz = 459.0

Nequilibrium = 25000000
Nproduction = 75000000
Ntotal  = Nequilibrium + Nproduction

context.initialize()

system = init.read_gsd(filename='init.gsd', restart='restart.gsd')

nl = md.nlist.cell()
#Set 1-4 exclusions
nl.reset_exclusions(exclusions = ['1-2', '1-3', '1-4'])

rcut = 28.00 # A

mie = md.pair.mie(r_cut= rcut, nlist = nl)

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

#write restart each 100k steps
restart = dump.gsd(filename="restart.gsd", group = all, truncate=True, period = 100000, phase=0)

update.balance()

T = 344.15 #K
kT = T * kBby10
isothermal = md.integrate.nvt(group = all, kT = kT, tau = 1.0)
run_upto(Nequilibrium)
isothermal.disable()

period = 1000
#data to be logged
loglist = loglist = ['potential_energy', 'kinetic_energy',
                     'temperature', 'pressure_zz', 'pressure_xx', 'pressure_yy']
log = analyze.log(filename = 'log.dat', quantities = loglist, period = period, overwrite=False, phase = -1)


nve = md.integrate.nve(group = all)
run_upto(Ntotal)

