#Importing all the cool functions I couldn't code myself
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import pandas as pd


#Creating a function to determine the distance between 2 points in cartisian space
#And by create I mean take the bond length function I created last home work and change a few values
def Atomic_Distance(Atom1 = [0, 0, 0], Atom2 = [0, 0, 0]):
    #distance starts at 0
    distance = 0
    # The for loop adds the square difference in position in all dimentions together
    for position in range(min(len(Atom1),len(Atom2))):
        diff = (Atom1[position] - Atom2[position])**2
        distance += diff
    return np.sqrt(distance)

    #A simple function that takes an input of r (and potentially epsilon and sigma) and outputs the atomic energy
def lennard_jones(r, epsilon=0.01, sigma=3.4):
    
    Atomic_Energy = 4*epsilon*((sigma/r)**12-(sigma/r)**6)
    return Atomic_Energy

distance_between_two = opt.minimize(
    #This is the function
    fun=lennard_jones,  # Objective function to minimize
    #This is the initial guess
    x0=0,                    # Initial guess
    #This is the method, I don't remeber the difference between them so I'm using the first one taught in class
    method="Nelder-Mead",    
    #This is the tolernace because computers have yet to acheve perfection
    tol=1e-6                 
)

print(distance_between_two)