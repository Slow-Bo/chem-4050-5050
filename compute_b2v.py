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
def hard_sphere(r):
    #The trapazoid function runs an array through the functions to be intigrated so this is to allow for that
    if isinstance(r, np.ndarray):  
        #Sees where each number in the array would fall in the hard sphere argument and appends the output accordingly
        for n in r:
            #This is needed to target each number in the array, enumerate would be better, but I know how to use this
            #So long as the same n value dosn't appear twice (Which it shouldn't) it will work fine
            Counter = np.where(r == n)
            #Modifys each number in the array to be what it needs to be since transfering the modified numbers to a new array dosn't work
            if n < Sigma:
                r[Counter] = 1000
            else:
                r[Counter] = 0
        return r
    else:
        #Just sees where the number would fall in the hard sphere argument
        if r < Sigma:
            return 1000
        else:
            return 0

#One more function for square well potential
def square_well(r):
    #The trapazoid function runs an array through the functions to be intigrated so this is to allow for that
    if isinstance(r, np.ndarray):
        for n in r:
            #This is needed to target each number in the array, enumerate would be better, but I know how to use this
            #So long as the same n value dosn't appear twice (Which it shouldn't) it will work fine
            Counter = np.where(r == n)
            #Modifys each number in the array to be what it needs to be since transfering the modified numbers to a new array dosn't work
            if n < Sigma:
                r[Counter] = 1000
            elif Sigma <= n < (WellRange * Sigma):
                r[Counter] = -Epsilon
            elif n > (WellRange * Sigma):
                r[Counter] = 0
        return r
    else:
        #Sees where each number in the array would fall in the square well potential argument and appends the output accordingly
        
            #Just sees where the number would fall in the square well potential argument
        if r < Sigma:
            return np.inf
        elif Sigma <= r < WellRange * Sigma:
            return -Epsilon
        else:
            return 0
        

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
#B2VHSat100 = bind_arguments(B2V, Square_well, 100)

#This just prints the arguments to see if they work
#Goal = TheIntagrateInator(B2VHSat100,1/1000, 17.0, 1000)
#print(Goal)

#Creates the range of temperatures over which B2V will be evaluated
x = np.linspace(100, 800, 1000)

#Because of the way I bound functions together the linespace only generates 1 value for y, when I use x but I have a solution
#Defining the arrays to be graphed
yhs = []
ysw = []
ylj = []
#Runs a for loop intigrating each value of x one at a time

for i in range(len(x)):
#Finds the integrals at each value X (assuming it dosn't blow up my laptop)
#also divides all the functions by their order of magnitude so they appear on the same graph
    yhs.append(TheIntagrateInator(bind_arguments(B2V,hard_sphere,x[i]),1/1000, 5*Sigma, 1000)/10**30)
    ysw.append(TheIntagrateInator(bind_arguments(B2V,square_well,x[i]),1/1000, 5*Sigma, 1000)/10**30)
    ylj.append(TheIntagrateInator(bind_arguments(B2V,lennard_jones,x[i]),1/1000, 5*Sigma, 1000)/10**23)


#creates a plot of the polynomial and the original data.
plt.plot(x, yhs, linestyle = '-',marker='', color='blue', label='Hard Sphere Function/10^30')
plt.plot(x, ysw, linestyle = '-',marker='', color='red', label='Square Well Function/10^30')
plt.plot(x, ylj, linestyle = '-',marker='', color='green', label='Lennard Jones Function/10^23')
plt.axhline(0,0,Sigma*5,linestyle = '-',marker='', color='black', label='B2V = 0')
plt.xlabel('Temperature')
plt.ylabel('Second virial coefficient')
plt.legend()
plt.title('The second viral coeffeicent using different functions to calculate potental energy')
plt.show()

"""The first thing to notice is that the second viral coeffeicent for the Square Well and Hard Sphere potentials were much more simmilar
than the lennard jones function. This makes sense as outside of the titular square well they have the same value. 
They were also less efected by temperature, which makes since sense the factor including temperature is always irrelevent for hard sphere,
Since it is always multiplied by either 0 or infinity, and it is almost irrelevent for the square well since it is only important within,
The square well, and even then not super important since epsilon is so small. The function that varies the most with temperature is the Lannard Jones function
This function actually produces a reasonable value at each distance. It is also several orders of magnitude smaller then the other two probably because it
lacks a large range where it goes off to infinity. """