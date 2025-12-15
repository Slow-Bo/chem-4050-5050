#Imports all the functions I can't code as well as a few constants
import numpy as np
import scipy.constants as con
import matplotlib.pyplot as plt

#The constants my functions need to, well function
Kb = con.Boltzmann
eV = con.eV
Amu = con.atomic_mass
#I couldn't figure out the function to bind a value to a certain range so I made my own
def Modulalo(X,Box_Size):
    
    if X > Box_Size:
        #If Xmax is smaller than X it loops back to 0 and what is left of X is the positive displacemnt
        return X % Box_Size
    elif X < 0:
        #If X is negitive it loops back to the maximum value and what reamins of X is the negitive displacemnt
        return X % Box_Size
    else:
        #If nither of the above are true X is valid so it just returns X
        return X
    #Function to determine the neighbors of each cell

#A Function that creates the intital particle chain
def create_chain(length:int,box_size:float, r0: float):
    #The base array is the intial chain of zeros
    Base_array = np.zeros((length, 3))
    #Finds the core of the box
    Current_Position = [box_size/2,box_size/2,box_size/2]
    #spreads out the particles by an initial radius
    for particle in range(len(Base_array)):
        #Ensures the original particle is in the box's core
        if particle == 0:
            Base_array[particle] =  [box_size/2,box_size/2,box_size/2]
        else:
            #calculates the random direction to seperate the particle from the center
            x = np.random.rand() * r0 * np.random.choice([1,-1])
            y = np.random.rand() * np.sqrt(r0**2-x**2) * np.random.choice([1,-1])
            z = np.sqrt(r0**2-x**2-y**2) * np.random.choice([1,-1])
            #Displaces the particle the chosen amount from the previous particle
            new_particle = np.add(Current_Position,[x,y,z])
            #Shoves X Y and Z into the boundry conditions
            for position in new_particle:
                position = Modulalo(position,box_size)
            Current_Position = new_particle
            Base_array[particle] = new_particle
    return Base_array

#A function to assimble our initial velocities
def initial_velocities(length:int,Target_temp:float,mass:float):
    #Converts the mass from Amu to Kg so that units cancel in the scale factor:
    mass_Kg = mass * Amu
    # Give the particles random velocities
    velocities = np.random.uniform(-0.5, 0.5, (length, 3))
    # Calculate the center of mass velocity
    v_com = np.sum(velocities, axis=0) / length
    # Subtract the center of mass velocity from the velocities
    velocities -= v_com
    #A Boltzmann distribution requires the standard devation to be sqrt(kb*T/m) this finds the current std so that can be obtained
    CalcSTD = np.sqrt(np.sum(velocities**2)/(length-1))
    #This Scale factor will ensure that the standard deviation will always be sqrt(kb*T/m), assuming I did my math right
    Scale_factor = np.sqrt(Kb*Target_temp/(mass_Kg*CalcSTD**2))
    velocities *= Scale_factor
    #After running the print a few times it seems to be working however the mass must be in Kg to work I changed the variable to mass Kg so I won't forget that
    #I probably will still forget to do that
    return velocities
    #Judging by the psuedo code there was likely a way to get python to set the standard deviation for me
    #Whatever dispite being as graceful as a crocodile in a tutu, this method works

#A function to compute our harmonic forces
def compute_harmonic_forces(positions, k:float, r0:float, box_size:float):
    #Note to my future self k is the SPRING CONSTANT, not the Boltzmann constant
    #Creates a zero vector the same size as the position vector
    forces = np.zeros_like(positions)
    #Begins the calculate the force on each particle
    for particle in range(len(positions)-1):
        #Calculates the difference in two particles
        displacement = np.subtract(positions[particle+1],positions[particle])
        #This allows for the distance to pass between box boundries.
        #I did not figure this out myself unfortnatly since microsoft's wonderful AI decided to shove the correct answer in my face when I tried to research it
        #For what it's worth I did verify it would work on paper
        displacement = displacement - box_size * np.round(displacement / box_size)
        #Finds the radius of displacemnt by squaring it summing it and squarooting it in that order
        distance = np.sqrt(np.sum(displacement**2))
        #Calculates force magnitude in accordance with hooks law
        force_magnitude = -k * (distance - r0)
        #Calculates the force, the psudeocode does it in two, otherwise I would do it in one
        force = force_magnitude * (displacement/distance)
        #adds new values onto forces
        forces[particle] -= force
        forces[particle+1] += force
    return forces

#Based of this print it seems to be working, if r0 for initialization = r0 it prints values of almost 0
#The values otherwise seem weird to my monkey brain but they do make since based of the logic of everything in the middle being pulled in two directions
#print(compute_harmonic_forces(create_chain(7,7, 2),4,3,7))

#A function to calculate the lenard jones force
def lennard_jones_forces(positions, epsilon, sigma:float, box_size:float, interaction_type:str):
    #note to my future self epsilon 0 is ATTRACTIVE, and epsilon 1 is REPULSIVE
    #Creates a zero vector the same size as the position vector
    forces = np.zeros_like(positions)
    cutoff = 0.0
    #Prepares an e to be the chosen epsilon
    e = float
    #Sorts through all the particles
    for i in range(len(positions)):
        for j in range(i+1,len(positions)):
            #sees if particles are 2 apart or more for attractive and repulsive interactions respectivly
            if interaction_type == 'atr' and abs(i-j) == 2:
                e = epsilon[0]
                #Sets cutoff to box size since the cutoff is irrelevent for attractive forces
                cutoff = box_size
            elif interaction_type == 'rep' and abs(i-j) > 2:
                e = epsilon[1]
                #Sets cutoff to lj_minimum since repulsion dosn't exist beond that
                cutoff = sigma * 2**(1/6)
            else:
                continue
            #Calculates the difference in two particles
            displacement = np.subtract(positions[j],positions[i])
            #This allows for the distance to pass between box boundries.
            #I did not figure this out myself unfortnatly since microsoft's wonderful AI decided to shove the correct answer in my face when I tried to research it
            #For what it's worth I did verify it would work on paper
            displacement = displacement - box_size * np.round(displacement / box_size)
            #Finds the radius of displacemnt by squaring it summing it and squarooting it in that order
            distance = np.sqrt(np.sum(displacement**2))
            
            #Checks if the distance is less than cutoff, should be irrelevent for attractive forces
            if distance < cutoff:
                #Calculates the lenard jones forces and adds them to the forces array
                force_magnitude = 24 * e * ( (sigma / distance)**(12) - 0.5 * (sigma / distance)**6 ) / distance
                force = force_magnitude * (displacement / distance)
                forces[i] -= force
                forces[j] += force
    return forces

def calculate_harmonic_potental(positions,k,r0,box_size):
    Harmonic_Potentials = 0.0
    for particle in range(len(positions)-1):
        #Calculates the difference in two particles
        displacement = np.subtract(positions[particle+1],positions[particle])
        #This allows for the distance to pass between box boundries.
        #I did not figure this out myself unfortnatly since microsoft's wonderful AI decided to shove the correct answer in my face when I tried to research it
        #For what it's worth I did verify it would work on paper
        displacement = displacement - box_size * np.round(displacement / box_size)
        #Finds the radius of displacemnt by squaring it summing it and squarooting it in that order
        distance = np.sqrt(np.sum(displacement**2))
        #Calculates the contribution of this interaction to potental energy
        Harmonic_Potential = 0.5* k * (distance - r0)**2
        Harmonic_Potentials += Harmonic_Potential
    #Returns the total harmonic potential
    return Harmonic_Potential
    
#A function to calculate the lenard jones potential
def lennard_jones_potential(positions, epsilon, sigma:float, box_size:float, interaction_type:str):
    #note to my future self epsilon 0 is ATTRACTIVE, and epsilon 1 is REPULSIVE
    #Creates a zero variable to have the lenard jones energies
    lj_potential = 0.0
    cutoff = 0.0
    #Since the LJ potental has an extra bit for repulsion this acounts for that
    extra_bit = 0.0
    #Prepares an e to be the chosen epsilon
    e = float
    #Sorts through all the particles
    for i in range(len(positions)):
        for j in range(i+1,len(positions)):
            #sees if particles are 2 apart or more for attractive and repulsive interactions respectivly
            if interaction_type == 'atr' and abs(i-j) == 2:
                e = epsilon[0]
                #Sets cutoff to box size since the cutoff is irrelevent for attractive forces
                cutoff = box_size
            elif interaction_type == 'rep' and abs(i-j) > 2:
                e = epsilon[1]
                #Sets cutoff to lj_minimum since repulsion dosn't exist beond that
                cutoff = sigma * 2**(1/6)
                #prepares the extra bit for the LJ equation
                extra_bit = 1/4
            else:
                continue
            #Calculates the difference in two particles
            displacement = np.subtract(positions[j],positions[i])
            #This allows for the distance to pass between box boundries.
            #I did not figure this out myself unfortnatly since microsoft's wonderful AI decided to shove the correct answer in my face when I tried to research it
            #For what it's worth I did verify it would work on paper
            displacement = displacement - box_size * np.round(displacement / box_size)
            #Finds the radius of displacemnt by squaring it summing it and squarooting it in that order
            distance = np.sqrt(np.sum(displacement**2))
            if distance < cutoff:
                #adds the contribution of the potential to toal potential 
                potental_contribution = 4 * e * ((sigma/distance)**12-(sigma/distance)**6+extra_bit)
                lj_potential += potental_contribution
    #returns total potential
    return lj_potential

#A function to shove all our force functions into one for conveninence
def compute_forces(positions,k,r0,epsilon,sigma,box_size):
    #It just runs all my force funcctions, it isn't rocket science
    harmonic_forces_new = compute_harmonic_forces(positions,k,r0,box_size)
    lj_attractive_new = lennard_jones_forces(positions,epsilon,sigma,box_size,'atr')
    lj_repulisve_new = lennard_jones_forces(positions,epsilon,sigma,box_size,'rep')
    forces_new = np.add(np.add(lj_repulisve_new,lj_attractive_new),harmonic_forces_new)
    return forces_new

#A function to shove all our potential functions into one for conveninence
def compute_potential(positions,k,r0,epsilon,sigma,box_size):
    new_harmonics = calculate_harmonic_potental(positions,k,r0,box_size)
    new_LJ_atr = lennard_jones_potential(positions, epsilon, sigma, box_size, 'atr')
    new_LJ_rep = lennard_jones_potential(positions, epsilon, sigma, box_size,'rep')
    new_potential = new_LJ_rep + new_LJ_atr + new_harmonics
    return new_potential

#A function to calculate the verlet velocities
def velocity_verlet(positions, box_size, velocities, forces, dt, mass, k, r0, epsilon, sigma):
    #Sets the velocities to half the new forces over the time step
    velocities += 0.5 * forces / mass * dt
    #changes positions based on velocitiy
    positions += velocities * dt
    #Ensures nobody has left the box
    for position in range(len(positions)):
        for coord in range(len(positions[position])):
            positions[position][coord] = Modulalo(positions[position][coord], box_size)
    #Calculates the new forces based on positions
    forces_new = compute_forces(positions,k,r0,epsilon,sigma,box_size)
    #Sets the velocities to the new forces over the other half of the time step
    velocities += 0.5 * forces_new / mass * dt
    return positions, velocities, forces_new

#A function to do the anderson rescaling
def velocity_rescaling(velocities, target_temp, mass):
    n_particles = len(velocities)
    #converts mass from amu to kilograms so that units cancel
    mass_kg = mass * Amu
    #An empty number to be made into the summed absolute values of the velocities squared
    v_squared = 0.0
    #adds the r^2 value of each velocity to v_squared
    for velocity in velocities:
        v_squared += np.sum(velocity**2)
    kinetic_energy = 0.5 * mass_kg * v_squared
    #For units to cancel mass should be in kilograms if my calculations are correct
    temp_current = (2/3) * kinetic_energy / (n_particles * Kb)
    scaling_factor = np.sqrt(target_temp/temp_current)
    velocities *= scaling_factor
    return velocities

def radius_of_gyration(positions):
    radii = np.zeros_like(positions)
    mean_pot = np.sum(positions,axis=0) / len(positions)
    for position in range(len(positions)):
        radii[position] = (np.subtract(positions[position],mean_pot))**2
    Rg_squared = np.mean(np.sum((radii), axis=1))
    Rg = np.sqrt(Rg_squared)
    return Rg

def calculate_end_to_end_distance(positions):
    Ree = np.sqrt((np.sum(((np.subtract(positions[-1],positions[0]))**2))))
    return Ree

#A new function simply for changing the position for metropolis so there is no string BS, an idea I came up with unfortnatly after the last project
def random_change(positions):
    Valid_Sites = []
    for j, row in enumerate(positions):
        for k, element in enumerate(row):
            Valid_Sites.append((j, k))
    #Choses one of the valid sites at random and gets its X and Y
    Chosen_Site = np.random.choice(len(Valid_Sites))
    Chosen_X = Valid_Sites[Chosen_Site][0]
    Chosen_Y = Valid_Sites[Chosen_Site][1]
    #changes the positions by a small value
    change = np.random.uniform(-0.1,0.1)
    positions[Chosen_X][Chosen_Y] += change
    return positions

#Uses the metropolis algorythom to optimize the chain positions
def optimize_chain(positions,k,r0,epsilon,sigma,box_size,T,steps):
    beta = eV/(T*Kb)
    for i in range(steps):
        OldE = compute_potential(positions,k,r0,epsilon,sigma,box_size)
        new_positions = np.copy(positions)
        new_positions = random_change(new_positions)
        NewE = compute_potential(new_positions,k,r0,epsilon,sigma,box_size)
        Acc = min(1, np.exp(-beta * (NewE - OldE)))
        if np.random.rand() < Acc:
            positions = new_positions
    return positions

#Simulation parameters taken from the sample code, I will need to change mass since currently each paticle is a kilogram, but this will work for now
dt = 0.01  #Time step if I make this in NANOSECONDS I can make length be in NANOMETERS
optimization_steps = 1000 #The number of steps to be used in the optimizer
total_steps = 1000  #Number of steps
box_size = 100.0  #Size of the cubic box
k = 1.0  #Spring constant
#I changed the mass to 3.01E-25 since it is the weight of tyrosine in kilograms, which as a heavier amino acid would be a good stand in for whatever our monomer is.
#Now it dosn't want to work, i'm changing it back
#I figured it out, only the inital velocities function needs mass in kilograms, but everything else needs mass in daltons, so no matter what I put here it won't work,
#Unless I bite the bullet and fix the inital velocities function
#Okay I fixed it, now I will probably still set the mass to that of tyrosine since we aren't creating a chain of hydrogents after all
mass = 181.07  #Particle mass NOT IN KILOGRAMS
r0 = 1.0  #Equilibrium bond length, I am not increasing this since these units should be, no this would be in meters right? if dt is in seconds anyway, what would dt need to be for this to make sense I wonder?
# I figured it out, since the velocity is in meters/second if I say dt is in nanoseconds I can say this is nanometers, which is reasonable enough 10 angstrom, works for me.
target_temperature = 0.1  #Target temperature
rescale_interval = 50  #Steps between velocity rescaling
n_particles = 20  #Number of particles
epsilon_repulsive = 1.0  #Depth of repulsive LJ potential
epsilon_attractive = 0.5  #Depth of attractive LJ potential
sigma = 1.0  #LJ potential parameter
epsilon = [epsilon_attractive,epsilon_repulsive]

#sets the seed for repeatability
np.random.seed(42)
# Arrays to store properties
temperatures = np.linspace(0.1, 1.0, 10)
Rg_values = []
Ree_values = []
potential_energies = []

print(calculate_end_to_end_distance([[50.3049695, 49.76109065, 50.37665302],
 [50.06285129, 50.15102375, 51.26153974],
 [49.68545273, 49.40922583, 51.78757349],
 [48.95969156, 49.97998735, 51.4720916 ],
 [48.33680258, 49.2747441 , 51.13023187],
 [48.60032258, 48.87728036, 50.24523635],
 [49.37245716, 49.25410984, 50.72157495],
 [48.77663822, 49.97584799, 50.36794342],
 [48.53051368, 50.34794417, 49.48617786],
 [49.50468808, 50.46835467, 49.6735492 ],
 [49.27219262, 49.51020925, 49.47240071],
 [49.67999522, 49.95365739, 48.68951136],
 [48.81983376, 49.48548245, 48.42391742],
 [49.46170906, 49.51551732, 47.68056736],
 [48.69645887, 49.09907407, 47.22132909],
 [48.61994052, 48.59530436, 46.32118022],
 [48.42093551, 48.01822829, 47.08871171],
 [48.42637096, 47.47488602, 46.27048052],
 [49.1842839,  47.79114101, 45.73506546],
 [49.43553884, 47.75168136, 46.69224022]]))

for T in temperatures:
    #Sets target temperature to T and creates an empty array
    target_temperature = T
    potential_energy_array = []
    #Creates an array named positions
    init_positions = create_chain(n_particles, box_size, r0)
    #Optimizes said array with metropolis
    optimized_positions = optimize_chain(init_positions,k,r0,epsilon, sigma, box_size,T, optimization_steps)
    #Calculates all velocities
    velocities = initial_velocities(n_particles, target_temperature, mass)
    for step in range(total_steps):
        forces = compute_forces(optimized_positions,k,r0,epsilon,sigma,box_size)
        #Finds the new forces and velocites
        optimized_positions, velocities, forces = velocity_verlet(optimized_positions, box_size, velocities, forces, dt, mass, k, r0, epsilon, sigma)
        #Calculates the potential energy and adds it to the array
        new_potential = compute_potential(optimized_positions,k,r0,epsilon,sigma,box_size)
        potential_energy_array.append(new_potential)
        #Apply thermostat
        if step % rescale_interval == 0:
            velocities = velocity_rescaling(velocities, target_temperature, mass)
    # Compute properties and adds them to the array
    Ree = calculate_end_to_end_distance(optimized_positions)
    Rg = radius_of_gyration(optimized_positions)
    Rg_values.append(Rg)
    Ree_values.append(Ree)
    potential_energies.append(np.mean(potential_energy_array))


# Plotting I took from the project page, it's pretty basic, I'm just lazy and didn't feel like coding it myself.
plt.figure()
plt.plot(temperatures, Rg_values, label='Radius of Gyration')
plt.xlabel('Temperature: (K)')
plt.ylabel('Radius of Gyration: (nm)')
plt.title('Radius of Gyration vs Temperature')
plt.legend()
plt.show()

plt.figure()
plt.plot(temperatures, Ree_values, label='End-to-End Distance')
plt.xlabel('Temperature: (K)')
plt.ylabel('End-to-End Distance: (nm)')
plt.title('End-to-End Distance vs Temperature')
plt.legend()
plt.show()

plt.figure()
plt.plot(temperatures, potential_energies, label='Potential Energy')
plt.xlabel('Temperature: (K)')
plt.ylabel('Potential Energy: (eV)')
plt.title('Potential Energy vs Temperature')
plt.legend()
plt.show()
