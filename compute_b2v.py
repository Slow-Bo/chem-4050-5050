#Importing all the cool functions I couldn't code myself
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
from scipy.integrate import trapezoid
import pandas as pd

#Avogadro's Number
AvN = 6.022140857E23

#Boltzmann Constant (In electron volts)
kB = 8.6173324E-5

#Values to be kept constant for this simulation
Sigma = 3.4
Epsilon = 0.01
WellRange = 1.5

#A simple function that takes an input of r (and potentially epsilon and sigma) and outputs the atomic energy
def lennard_jones(r, epsilon = Epsilon, sigma = Sigma):
    
    Atomic_Energy = 4*epsilon*((sigma/r)**12-(sigma/r)**6)
    return Atomic_Energy

#Another function for hard sphere potential
def Hard_Sphere(r):
    #The trapazoid function runs an array through the functions to be intigrated so this is to allow for that
    if r == float:
        #Just sees where the number would fall in the hard sphere argument
        if r < Sigma:
            return np.inf
        else:
            return 0
    else:
        #Sees where each number in the array would fall in the hard sphere argument and appends the output accordingly
        Output = []
        for n in r:
            if n < Sigma:
                Output.append(np.inf)
            else:
                Output.append(0)
        return Output

#One more function for square well potential
def Square_well(r):
    #The trapazoid function runs an array through the functions to be intigrated so this is to allow for that
    if r == float:
        #Just sees where the number would fall in the square well potential argument
        if r < Sigma:
            return np.inf
        elif Sigma <= r < WellRange * Sigma:
            return -Epsilon
        else:
            return 0
    else:
        #Sees where each number in the array would fall in the square well potential argument and appends the output accordingly
        Output = []
        for n in r:
            if n < Sigma:
                Output.append(np.inf)
            elif Sigma <= n < WellRange * Sigma:
                Output.append(-Epsilon)
            else:
                Output.append(0)
        return Output
        

#The function to be intigrated takes two inputs f (a function for potental) and r (distance between two atoms)
def B2V(f, TempK, r):
    Virtual = -2 * np.pi * AvN * (np.e**((-2 * f(r))/(kB * TempK)) -1) * r**2
    return Virtual

def TheIntagrateInator(f, start, end, n):
    ''' A function that intigrates another function using scipy's Trapazoid function
    
    Parameters:
    f is the function to be integrated
    start is the start bounds of the function
    end is the end bounds of the function
    n is the number of points for the integral evaluation

    Returns:
    Float: The integral of the function on the defined bounds

    '''
    x = np.linspace(start,end,n, endpoint= False)
    dx = (end-start) / n

    integral = trapezoid(f(x),x,dx, axis=0)
    return integral

#I could not figure out how to pass functions with empty variables into other functions, I found this solution online, I will explain how I think it works
#This first part takes the function you need some arguments for and the arguments you want to fill. I think the * makes them optional
def bind_arguments(func, *args, **kwargs):
    #This argument is the function that will be returned when all is done
    def bound_function(*inner_args, **inner_kwargs):
        #This is the part of the function that makes the combined function where the inner args are left blank and the outer args are not
        return func(*args, *inner_args, **kwargs, **inner_kwargs)
    #This returns the combined function
    return bound_function

#Also heres my source because citation is important https://www.geeksforgeeks.org/python/how-to-bind-arguments-to-given-values-in-python-functions/

#This function combines the B2V function with one of the potential functions
B2VHSat100 = bind_arguments(B2V, Hard_Sphere, 100)
Goal = TheIntagrateInator(B2VHSat100,1/1000, 17.0, 1000)
print(Goal)
