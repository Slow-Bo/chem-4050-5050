#Importing all the cool functions I couldn't code myself
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import pandas as pd

#a function to calculate bond length the input is 2 atom's positions
def Bond_Length(Atom1 = [0, 0, 0], Atom2 = [0, 0, 0]):
    #distance starts at 0
    distance = 0
    # The for loop adds the square difference in position in all dimentions together
    for position in range(min(len(Atom1),len(Atom2))):
        diff = (Atom1[position] - Atom2[position])**2
        distance += diff
    return np.sqrt(distance)
    #Checks to see if bond length is valid
    if 0 < np.sqrt(distance) < 1.5:
        # returns the square root of the diffrences squared as per pythagoras' theorm
        return np.sqrt(distance)
    else:
        pass

#A function for subtracting one vector from aonther
def VecSubtract(vector01, vector02):
    #creates a new vector using the first 3 values of two other vectors
    vector03 = [vector01[value] - vector02[value] for value in range(min(len(vector01), len(vector02)))]
    return vector03
#A function for computing the bond angle
def Bond_Angle(Atom1 = [0, 0, 0], Atom2 = [0, 0, 0], Atom3 = [0, 0, 0]):
    #These values will be used in calculations
    Numerator = 0
    #Calls my subtracting function because I am too lazy to write the same thing twice
    Vector_1 = VecSubtract(Atom3, Atom1)
    Vector_2 = VecSubtract(Atom2, Atom1)
    #Figures out the absolute value of each vector by calling the bond length function
    AbsV1 = Bond_Length(Atom1, Atom3)
    AbsV2 = Bond_Length(Atom2, Atom1)
    #Checks for invalid bond lengths
    if AbsV1 == None or AbsV2 == None:
        pass
    else:
    #Figures out the dot product of the vector
        for value in range(min(len(Vector_1),len(Vector_2))):
            Numerator += np.multiply(Vector_1[value], Vector_2[value])
    #Plugs all priorly calculated numbers into the equation to calculate bond angle in degrees
        Result = np.degrees(np.arccos(Numerator/(AbsV1*AbsV2)))

        #I kept getting 0 angles so I added this redundancy to stop it
        if Result < 0.001:
            return
        #Prints what type of angle the bond angle is
        elif Result < 90:
         print("The angle is acute")
        elif Result == 90:
            print("The angle is right")
        else:
            print("The angle is obtuse")

        return (Result)

    #A simple function that takes an input of r (and potentially epsilon and sigma) and outputs the atomic energy
def lennard_jones(r, epsilon=0.01, sigma=3.4):
    
    Atomic_Energy = 4*epsilon*((sigma/r)**12-(sigma/r)**6)
    return Atomic_Energy

distance_between_two = opt.minimize(
    #This is the function
    fun=lennard_jones,  # Objective function to minimize
    #This is the initial guess
    x0=0,                    
    #The arguments of the function scipy is not allowed to change
    args=(0.01, 3.4),
    #This is the method, I don't remeber the difference between them so I'm using the first one taught in class
    method="Nelder-Mead",    
    #This is the tolernace because computers have yet to acheve perfection
    tol=1e-6                 
)
