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
nve = md.integrate.nve(group = all)



Na = 6.022e23
nbins = 250
period = 1000

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
loglist = ['potential_energy', 'kinetic_energy',
                     'temperature', 'pressure_zz', 'pressure_xx', 'pressure_yy']
log = analyze.log(filename = 'log.dat', quantities = loglist, period = period, 
		overwrite=False, phase = -1)

run_upto(Ntotal)