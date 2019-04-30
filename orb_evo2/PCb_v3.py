import rebound
import numpy as np
import time

m_p = 4.50524e-6
r_p = 1.0e-4
M = 0.123
semi = 0.049

############################################

#Power Functions to be used
def rand_powerlaw(slope, min_v, max_v):
    y = np.random.uniform()
    pow_max = pow(max_v, slope+1.)
    pow_min = pow(min_v, slope+1.)
    return pow((pow_max-pow_min)*y + pow_min, 1./(slope+1.))

def rand_uniform(minimum, maximum):
    return np.random.uniform()*(maximum-minimum)+minimum

def rand_rayleigh(sigma):
    return sigma*np.sqrt(-2*np.log(np.random.uniform()))

def rand_radii(radii):
    return radii*np.random.uniform(1., 2.)

def calc_hill(semi,m,M):
        r_hill = semi*((m/(3.*M))**(1./3.))
        return r_hill

def get_roch(R):
#       r_roche = a*(((2.*M)/m)**(1./3.))
        r_roche = 2.44*R
        return r_roche

r_pl = float(get_roch(r_p))-0.0001

print str(r_p)
print str(calc_hill(semi,m_p,M))

sim = rebound.Simulation()

#integrator options
sim.integrator = "mercurius" #hermes, ias15

sim.dt = 0.00001     # make at least 1/10 period of closest orbit
sim.testparticle_type = 1
sim.ri_ias15.min_dt = 1e-6 # ensure that close encounters do not stall the integration

#collision and boundary options
sim.collision = "direct"
sim.collision_resolve = "merge"
sim.collision_resolve_keep_sorted = 1
sim.boundary = "open"
boxsize = 0.02 #slightly bigger than hill radius
sim.configure_box(boxsize)
sim.track_energy_offset = 1

#mass in solar and radius in AU
sim.add(m=m_p,r=r_pl) #Proxima Centauri b (1.5 Me, r=4.6846e-5)

#Moonlet disk
Mtot_disk =(1e-2)*m_p   #Mass of disk, 7.97e-6 similar to Pluto-Charon, 2e-4 similar to Jupiter, 1e-2 similar to Earth
halfdisk = Mtot_disk / 2      #Mass of half disk for 50/50 ratio
N_e = 20     #Number of lunar embryos
r_e = 2e-6     # radius of lunar embryos, 300 km
m_e = halfdisk/N_e
N_ml = 100          #Number of moonlets
r_ml = 6.69e-9          #Radius of each moonlet, 1km
m_ml = halfdisk/N_ml

np.random.seed()
# Inserting embryos into disk
while sim.N < (21):
    a = rand_powerlaw(0,get_roch(r_p),calc_hill(semi,m_p,M))
    e = rand_rayleigh(0.01)
    inc = rand_rayleigh(0.005)
    f = rand_uniform(-np.pi,np.pi) #Moonlets near PCb, (-pi, pi) spawns moolets all around center body
    sim.add(m=m_e,r=r_e,a=a,e=e,inc=inc,Omega=0, omega=0,f=f)
sim.N_active = sim.N

# Inserting moonlets into disk
while sim.N < (121):
    a = rand_powerlaw(0,get_roch(r_p),calc_hill(semi,m_p,M))
    e = rand_rayleigh(0.01)
    inc = rand_rayleigh(0.005)
    f = rand_uniform(-np.pi,np.pi) #Moonlets near PCb, (-pi, pi) spawns moolets all around center body
    r_m = rand_radii(r_ml)
    sim.add(m=m_ml,r=r_m,a=a,e=e,inc=inc,Omega=0, omega=0,f=f)
sim.testparticle_type = 1
sim.move_to_com()

E0 = sim.calculate_energy()

print sim.N_active
print sim.N

#### COMPLETE REVAMP OF STARTING SIMULATION ####

torb = 2.*np.pi     # 2*pi to convert to regular Gravitational Constant G
time_f = 3.5e1*torb    # integration time, 1.1e6*torb for PCb simulation

N_out = 1000     # Number of status outputs for text file
times = np.geomspace(sim.dt,time_f, N_out)     # Evenly spaced status outputs 

sim.status()

## STARTING SIMULATION / OUTPUTS FOR TEXT AND LINES 
for i, time in enumerate(times):
    sim.integrate(time)
    tst = str(sim.t)
    timeb = float(tst)
    norm_t = (timeb)/(torb)
    active = sim.N_active
    sim.move_to_com()
    print "#######################################################"
    print "### After "+str(norm_t)+" years, "+str(sim.N)+" particles remain. ###"
    print "#######################################################"
    print str(sim.dt)
    for i in range(1,active):
        orbit = sim.particles[i].calculate_orbit(sim.particles[0])
        mass = sim.particles[i].m
        print "Lunar embryo "+str(i)+" has a mass of "+str(mass)+" solar masses."
        print(orbit)
#    if sim.N < 21:
#        for i in range(1,sim.N):
#            orbit = sim.particles[i].calculate_orbit(sim.particles[0])
#            mass = sim.particles[i].m
#            print "Particle "+str(i)+" has a mass of "+str(mass)+" solar masses."
#            print(orbit)
#    else:
#        for i in range(1, 21):
#            orbit = sim.particles[i].calculate_orbit(sim.particles[0])
#            mass = sim.particles[i].m
#            print "Particle "+str(i)+" has a mass of "+str(mass)+" solar masses."
#            #f.write('\n' + "Particle "+str(i)+" has a mass of "+str(mass)+" solar masses.")
#            print(orbit)
#            #f.write('\n' + str(orbit) + '\n')

#####################################################

sim.status() 

#Calculates accuracy/energy loss
dE = abs((sim.calculate_energy() - E0)/E0)
print str(dE) + " energy loss"

fig = rebound.OrbitPlot(sim,slices=True,color=True,unitlabel="Semi-Major Axis (AU)",lim=0.0012, limz=0.001, figsize=(9,8), lw=2.0)
fig.savefig("PCb-com-test-2.png")
