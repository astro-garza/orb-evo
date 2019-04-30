import rebound
import numpy as np
import time

#Supercomment everything!

#### DELETED STAMPS ####

m_p = 1.49e-5
r_p = 6.39e-5
M = 1.037
semi = 1.046

#### ADD TEXT FILE TO IMPORT RESULTS ####
#file = open("452b_esim.txt", "w")
#file.close()
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

def calc_time(r_hill):
    t = r_hill * 100.e6
    return t


r_pl = float(get_roch(r_p))-0.0001

print str(r_p)
print str(calc_hill(semi,m_p,M))

def initcon():
    sim = rebound.Simulation()
    #integrator options
    sim.integrator = "hermes" #hermes, ias15
    sim.ri_hermes.hill_switch_factor = 3.
    sim.ri_hermes.solar_switch_factor = 20.
    
    #### CHANGED TIMESTEP: DELETE 2*PI ####
    sim.dt = 0.00001     # make at least 1/10 period of closest orbit
    #########################################
    
    sim.testparticle_type = 1
    #collision and boundary options
    sim.collision = "direct"
    sim.collision_resolve = "merge"
    sim.collision_resolve_keep_sorted = 1
    sim.boundary = "open"
    boxsize = calc_hill(semi,m_p,M)+0.02 #slightly bigger than hill radius
    sim.configure_box(boxsize)
    sim.track_energy_offset = 1
    #mass in solar and radius in AU
    sim.add(m=m_p,r=r_pl) #Proxima Centauri b (1.5 Me, r=4.6846e-5)

    #power functions for disk
    #Moonlet disk
    Mtot_disk =(1e-2)*m_p   #Mass of disk, 7.97e-6 similar to Pluto-Charon, 2e-4 similar to Jupiter, 1e-2 similar to Earth
    halfdisk = Mtot_disk / 2      #Mass of half disk for 50/50 ratio
    N_e = 20     #Number of lunar embryos
    r_e = 2e-6     # radius of lunar embryos, 300 km
    m_e = halfdisk/N_e
    N_ml = 100                          #Number of moonlets
    r_ml = 6.69e-9                         #Radius of each moonlet, 1km
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
    return sim

sim = initcon()
E0 = sim.calculate_energy()

print sim.N_active
print sim.N

#### COMPLETE REVAMP OF STARTING SIMULATION ####

torb = 2.*np.pi     # 2*pi to convert to regular Gravitational Constant G
time_f = calc_time(calc_hill(semi,m_p,M))*torb    # integration time, 1.1e6*torb for PCb simulation

N_out = 1000     # Number of status outputs for text file
times = np.geomspace(0.00001,time_f, N_out)     # Evenly spaced status outputs 

sim.status()

## STARTING SIMULATION / OUTPUTS FOR TEXT AND LINES 
for i, time in enumerate(times):
    sim.move_to_com()
    sim.integrate(time)
    sim.move_to_com()
    tst = str(sim.t)
    
    print "#######################################################"
    print "### After "+str(sim.t)+" years, "+str(sim.N)+" particles remain. ###"
    print "#######################################################"
    
    if sim.N < 21:
        for i in range(1,sim.N):
            orbit = sim.particles[i].calculate_orbit(sim.particles[0])
            mass = sim.particles[i].m
            print "Particle "+str(i)+" has a mass of "+str(mass)+" solar masses."
            print(orbit)
    else:
        for i in range(1, 21):
            orbit = sim.particles[i].calculate_orbit(sim.particles[0])
            mass = sim.particles[i].m
            print "Particle "+str(i)+" has a mass of "+str(mass)+" solar masses."
            print(orbit)

#####################################################

sim.status() 

#Calculates accuracy/energy loss
dE = abs((sim.calculate_energy() - E0)/E0)
print str(dE) + " energy loss"


