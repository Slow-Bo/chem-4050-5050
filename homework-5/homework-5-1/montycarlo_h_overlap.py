#Imports all the functions I can't code as well as a few constants
import numpy as np
import matplotlib.pyplot as plt
#Importing the trapazoid distribution since beta and exponentnal did not work.
from scipy.stats import trapezoid

#Setting the directory to our file so I can find the rest of my code, and so I save the graphs to the disired spot
#Dang I must have been really frustrated when microsoft's stupid CoPilot AI told me how to do this
import os
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)

#Setting my constants for when I need them
Seed = 42
R = 2
L = 20
BohrRad = 1

#A faunction to calculate the orbitals of the hydrogen atoms. I wanted to name this one orbital, but this is the name I'm suposed to give it.
def psi_2p_z(x, y, z):
        #Converting the coords into r so I don't have to type out this massive expression each time.
        r = np.sqrt(x**2 + y**2 + z**2)
        orbital = 1/(4*np.sqrt(2*np.pi)*BohrRad**(3/2))*(r/BohrRad)*(z/r)*np.e**-(r/(2*BohrRad))
        return orbital

#A function to calculate the monty carlo overlap
def Overlap(R,L,seed,n_point):

    #Set the random seed for reproducibility
    np.random.seed(seed)

    #Sets the integraton limits from 0 to L since the orbital should be semetrical on all axises? axii? Plural of axis!
    a = 0
    b = L

    #Creates a random array of X,Y, and Z points from the bounds 0 to L
    x = np.random.uniform(a, b, n_point)
    y = np.random.uniform(a, b, n_point)
    z = np.random.uniform(a, b, n_point)

    #Defines the integrand as the wave function squared as is how probability density is defined.
    #Thank god there are no imaginary numbers here, since otherwise I would throw my laptop into a lake.
    #Setting the z to +/- R over 2 so the atoms arn't on top of each other.
    integrand = psi_2p_z(x, y, z + R/2) * psi_2p_z(x, y, z - R/2)
    #Takes the integrand, multiplies it by the integration range cubed so it is a volume, then by 8 since it's semetrical on all axies
    integral = 8 * np.mean(integrand) * (b - a)**3
    #The same as before but for variance
    variance = 8 * np.var(integrand) * (b - a)**3

    #Returns the average and standard deviation
    return integral, np.sqrt(variance)

    #Not going to lie, I took a lot of this code from the lecture slides, I modified it to my own means because I personally beileve variables should be dictated outside the function.

#A function to visualize my integrand so I know which optimization function to use
def importance_vision():
        #I only need one axis to graph so the other two are held at 0-
        x = 0
        y = 0
        z = np.linspace(0, L, 100)
        integrand = psi_2p_z(x, y, z + R/2) * psi_2p_z(x, y, z - R/2)
        #Making the graph the new funny number aperently
        plt.figure(figsize=(6, 7))
        #Creating a plot to visualize the integrand. I will only be running this function durring testing.
        plt.plot(z, integrand, label='Integrand')
        plt.xlabel('x')
        plt.ylabel('Value')
        plt.title('The shape of our integrand')
        plt.legend()
        plt.show()

#Running the vison function
#importance_vision()

#Interestingly the function tends to peak around 2.5, is negitive around 0 and flatlines at 0 past 10 I need a function of similar shape
#After trial and error I have found that I can't get beta distribution to work, exponitial distribution would overprioritize the first part of the graph, so now I'll try trapazoid rule.

#A function to calculate the monty carlo overlap, because all the variables are defined outside, I can just copy and paste.
def Better_Overlap(R,L,seed,n_point):

    #Set the random seed for reproducibility
    np.random.seed(seed)

    #Seting the bounds of integration
    a = 0
    b = L

    #The way trapazoid works acording to the scipy docs is that I need to difine what percent of the way to L the trapazoid peaks, and where it starts going down again.
    #Based on my importance_vision() this is what I chose.
    c = 0.01
    d = 0.25

    #Creates a random array of X,Y, and Z points from the bounds 0 to L
    #If I did this right it should follow the beta distribution and be more acurate
    
    #It did not in fact work, as it turns out beta has 0 outside of its domain a lesser man would give up and use another function since it should be close enough.
    #I am a lesser man.
    x = trapezoid.rvs(c,d,loc=0,size=n_point, scale=L)
    y = trapezoid.rvs(c,d,loc=0,size=n_point, scale=L)
    z = trapezoid.rvs(c,d,loc=0,size=n_point, scale=L)

    #Defines the integrand as the wave function squared as is how probability density is defined.
    #Thank god there are no imaginary numbers here, since otherwise I would throw my laptop into a lake.
    #Setting the z to +/- R over 2 so the atoms arn't on top of each other.
    integrand = psi_2p_z(x, y, z + R/2) * psi_2p_z(x, y, z - R/2)

    #Apperently you have to normalize the integrand, don't ask how long it took me to figure that out.
    denominator = trapezoid.pdf(x,c,d,loc=0, scale=L) * trapezoid.pdf(y,c,d,loc=0, scale=L) * trapezoid.pdf(z,c,d,loc=0, scale=L)
    norm_integ = integrand/denominator
    #Takes the integrand, multiplies it by the integration range cubed so it is a volume, then by 8 since it's semetrical on all axies
    integral = 8 * np.mean(norm_integ)
    #The same as before but for variance
    variance = 8 * np.var(norm_integ)

    #Returns the average and standard deviation
    return integral, np.sqrt(variance)

