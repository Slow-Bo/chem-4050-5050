#Imports all the functions I can't code as well as a few constants
import numpy as np
import matplotlib.pyplot as plt
#Importing the expontental function since it is the best function
from scipy.stats import expon

#Setting the seed for the function
seed = 42

#A function for calculating the kenetic wave function using x,y, and z as variables
def psi_1s(x, y, z, Z=1, a0=1):
    r = np.sqrt(x**2+y**2+z**2)
    n_wav = 1/np.sqrt(np.pi*a0**3)*np.e**-(r/a0)
    return n_wav
#A function for calculating the laplacian operator of the wave function using x,y, and z as variables
def laplacian_psi_1s(x, y, z, Z=1, a0=1):
    r = np.sqrt(x**2+y**2+z**2)
    n_wav = psi_1s(x,y,z)
    ddr1 = np.gradient(n_wav,r)
    laplacian = np.gradient(ddr1,r) + 2/r * ddr1
    return laplacian

#The function for calculating the kenetic wave function using importance sampling and three inputs: The distance of integration, the number of disired points, and an optional R array for the off diagonal
def monty_carlo_ken(L,N,R=[0,0,0]):

    #Set the random seed for reproducibility
    np.random.seed(seed)

    #Sets the integraton limits from -L to L since although the orbital is semetrical for the diagonal, it isn't for the off diagonal.
    a = -L
    b = L

    #Creates a random array of X,Y, and Z points from the bounds 0 to L
    x = np.random.uniform(a, b, N)
    y = np.random.uniform(a, b, N)
    z = np.random.uniform(a, b, N)

    #Defines the integrand as the kenetic function how it is defined in the homework
    integrand = psi_1s(x + R[0]/2,y + R[1]/2 ,z + R[2]/2) * (-1/2 * laplacian_psi_1s(x - R[0]/2 ,y - R[1]/2, z - R[2]/2))
    #Takes the integrand, multiplies it by the integration range cubed so it is a volume, then by 8 since it's semetrical on all axies according to the importance vision
    integral = np.mean(integrand) * (b - a)**3
    #The same as before but for variance
    variance = np.var(integrand) * (b - a)**3

    #Returns the average and standard deviation
    return integral, np.sqrt(variance)

def importance_vision(L):
        #I only need one axis to graph so the other two are held at 0-
        x = 0
        y = np.linspace(0, L, 100)
        z = 0
        #I convert the integrand to absolute value since that is what I care most about
        integrand = (psi_1s(x, y, z))**2
        Test_dis = np.e**(-y)/3
        #Making the graph the new funny number aperently
        plt.figure(figsize=(6, 7))
        #Creating a plot to visualize the integrand. I will only be running this function durring testing.
        plt.plot(y, integrand, label='Integrand')
        #the function to be tested
        plt.plot(y,Test_dis, label='Test')
        plt.xlabel('x')
        plt.ylabel('Value')
        plt.title('The shape of our integrand')
        plt.legend()
        plt.show()

#importance_vision(7)
#This importance sampling showed me that I was wasting time mirroring on only two axies since the function is semetrical on all axies
#I'll go back and fix that, but looking at the function it appears that as though the exponental function of scale 3 would be preferable to trapazoid.
#But you do recomend gausian so I will give it a shot
#I tried, I couldn't figure out the probability density

#The function for calculating the kenetic wave function using importance sampling and three inputs: The distance of integration, the number of disired points, and an optional R array for the off diagonal
def better_monty_carlo_ken(L, n_point, R=[0,0,0]):

    #Set the random seed for reproducibility
    np.random.seed(seed)

    #Sets the integraton limits from -L to L since although the orbital is semetrical for the diagonal, it isn't for the off diagonal.
    a = -L
    b = L
    #A single variable to adjust the scale of the exponetntial to match the function best
    S = 3
    #Creates a random array of X,Y, and Z points from the bounds 0 to L
    
    x = expon.rvs(size=n_point, scale=S)
    y = expon.rvs(size=n_point, scale=S)
    z = expon.rvs(size=n_point, scale=S)
    
    #Defines the integrand as the kenetic function how it is defined in the homework, contains an optional input for R if the offdiagonal matrix is disired to be calculated
    integrand = psi_1s(x + R[0]/2,y + R[1]/2 ,z + R[2]/2) * (-1/2 * laplacian_psi_1s(x - R[0]/2 ,y - R[1]/2, z - R[2]/2))

    #Apperently you have to normalize the integrand, don't ask how long it took me to figure that out.
    denominator = expon.pdf(x,scale=S) * expon.pdf(y,scale=S) * expon.pdf(z,scale=S)
    norm_integ = integrand/denominator
    #Takes the integrand, multiplies it by the integration range cubed so it is a volume, then by 8 since it's semetrical on all axies
    integral = np.mean(norm_integ)
    #The same as before but for variance
    variance = np.var(norm_integ)

    #Returns the average and standard deviation
    return integral, np.sqrt(variance)

#What remains of my failure to calculate the integral with the gausian of R, I put a good amount of time into it so I couldn't bring myself to delete it outright
'''def even_better_monty_carlo_lap(L, n_point, R=[0,0,0]):

    #Set the random seed for reproducibility
    np.random.seed(seed)

    #Seting the bounds of integration
    a = 0
    b = L
    
    #Creates a random array of r from 0 t0 sqrt(3)*L
    r = np.random.normal(0,np.sqrt(3)*L,n_point)
    #Finds r's variance and mean which will be needed for the probability density later
    rAvg = np.mean(r)
    rVar = np.var(r)
    #creates x y and z values from r
    x = r/np.sqrt(3)
    y = r/np.sqrt(3)
    z = r/np.sqrt(3)

    #Defines the integrand as the kenetic function how it is defined in the homework
    integrand = psi_1s(x + R[0]/2,y + R[1]/2 ,z + R[2]/2) * (-1/2 * laplacian_psi_1s(x - R[0]/2 ,y - R[1]/2 ,z = R[2]/2))

    #Apperently you have to normalize the integrand, don't ask how long it took me to figure that out.
    #This should be r's probability density according to wikipedia anyway
    denominator = 1/np.sqrt(2*np.pi*rVar) * np.e**-((np.sqrt(x**2+y**2+z**2)-rAvg)**2/(2*rVar))
    norm_integ = integrand/denominator
    #Takes the integrand, since r already acounts for all axies, there is no need for multiplication
    integral = np.mean(norm_integ)
    #The same as before but for variance
    variance = np.var(norm_integ)

    #Returns the average and standard deviation
    return integral, np.sqrt(variance)'''

#A print for testing
print(monty_carlo_ken(7,10**7, [0,0,1.4]))
#print(better_monty_carlo_ken(7,10**7))
