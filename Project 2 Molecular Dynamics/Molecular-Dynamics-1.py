#Imports all the functions I can't code as well as a few constants
import numpy as np
import scipy.constants as con

#The constants my functions need to, well function
Kb = con.Boltzmann

#I couldn't figure out the function to bind a value to a certain range so I made my own
def Modulalo(X,Box_Size):
    
    
    if X > Box_Size:
        #If Xmax is smaller than X it loops back to 0 and what is left of X is the positive displacemnt
        return X % Box_Size
    elif X < 0:
        #If X is negitive it loops back to the maximum value and what reamins of X is the negitive displacemnt
        return Box_Size + (X % Box_Size)
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
def initial_velocities(length:int,Target_temp:float,mass_Kg:float):
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
def compute_harmonic_positions(positions, k, r0, box_size):
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
print(compute_harmonic_positions(create_chain(7,7, 2),4,3,7))