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

Nequilibrium = 5000000
Nproduction = 12000000
Ntotal  = Nequilibrium + Nproduction

context.initialize("")

system = init.read_gsd(filename='init.gsd', restart='restart.gsd')


nl = md.nlist.cell()
#Set 1-4 exclusions
nl.reset_exclusions(exclusions = ['1-2', '1-3', '1-4'])

rcut = 6*sigma # A

#Setting Mie interaction parameters
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
standard = md.integrate.mode_standard(dt=0.003)

#write restart each 100k steps
restart = dump.gsd(filename="restart.gsd", group = all, truncate=True, 
		period = 100000, phase=0)

T = 500.0 #K
kT = T * kBby10
isothermal = md.integrate.nvt(group = all, kT = kT, tau = 1.0)
run_upto(Nequilibrium, quiet = False)
isothermal.disable()
nve = md.integrate.nve(group = all)

Na = 6.022e23
nbins = 250 #number of bins to compute density profile
period = 1000 #period for log and density profile

dz = Lz/(nbins -1)
vbin = Lx * Ly * dz #A3

def density_profile(timestep):
    snap = system.take_snapshot()
    if comm.get_rank() == 0:     
        f = open('zdensity','ab')
        posz = snap.particles.position[:, 2]
        h, _ = np.histogram(posz, bins = nbins - 1, range=(-Lz/2, Lz/2))
        den =  h / ms #divide by chain lenght -> molecule/bin
        den /= vbin #divide by bin volume  -> molecule/A3
        den *= 10**27 / Na #convert to mol/l
        np.savetxt(f, den, newline=' ', delimiter=',')
        f.write(b'\n')
        f.close()
#To read the file use np.loadtxt('zdensity')
    
callback = analyze.callback(callback = density_profile, period = period, phase = -1 )


zbins = np.linspace(-Lz/2, Lz/2, nbins)
zgroups = {}
zcompute = {}
stensorlist = []
for i in range(nbins-1):
    zname = 'z' + str(i)
    zgroups[zname] = group.cuboid(zname, zmin = zbins[i], zmax = zbins[i+1])
    zcompute[zname] = compute.thermo(zgroups[zname])
    stensorlist .append('pressure_xx_'+ zname)
    stensorlist .append('pressure_yy_'+ zname)
    stensorlist .append('pressure_zz_'+ zname)

def group_update(timestep):
    for group in zgroups.values():
        group.force_update()

#These callbacks are mean to update the group in the step period-1 and period
group_update1 = analyze.callback(callback = group_update, period = period, phase = -1)
run(1, quiet = True)
group_update2 = analyze.callback(callback = group_update, period = period, phase = -1)
pressurelog = analyze.log(filename = 'pressure.dat', quantities = stensorlist , period = period,
               overwrite=False, phase = -1)


#data to be logged
loglist = loglist = ['potential_energy', 'kinetic_energy',
                     'temperature', 'pressure_zz', 'pressure_xx', 'pressure_yy']
log = analyze.log(filename = 'log.dat', quantities = loglist, period = period, 
		overwrite=False, phase = -1)

run_upto(Ntotal)