import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Prepares the values for calculation
h = 1
m = 1
L = 40
w = 1
D = 10
Beta = np.sqrt(1 / (2 * D))

#Creates a range of x values from -L/2 to L/2
Real_Space_Grid = np.linspace(-L/2, L/2, 2000)

#Takes the initial X for the anti Harmonic Equation, and deltaX for the Laplacian
x0 = Real_Space_Grid[0]
DeltaX = np.mean(np.diff(Real_Space_Grid))

#Calculates V Harmonic and anti Harmonic at each x
VHarmonic = 1/2 * m * w**2 * Real_Space_Grid ** 2

VAntiHarm = D * (1 - np.e**(-Beta * (Real_Space_Grid - x0)))**2

#Prepares the Matricies for the energy levels

MHarmonic = np.diagflat(VHarmonic)

MAntiHarm = np.diagflat(VAntiHarm)

#Constructs the Off Diagonal Matrix by preparing a 1999 long array of ones and imposing them on the super and sub diagonal
#I mean I just took my code from homework one. It seemed to work there
Ones = np.full(1999, 1)
Sup_Diag_Matrix = np.diagflat(Ones, k=1)
Sub_Diag_Matrix = np.diagflat(Ones, k=-1)

#Combines Off Diagonal Matrix with the Identity matrix times negitive 2 to create the Laplacian Matrix
LP_Matrix = Sup_Diag_Matrix + (np.identity(2000) * -2) + Sub_Diag_Matrix

#This will take the matrix times the one over the distance between points squared as requested
Laplacian = (1/(DeltaX**2)) * LP_Matrix

#This will construct the hamaltonian matracies for the harmonic and antiharmonic equations by first creating the hamiltonian and then adding to both original matricies
Hamiltionian = -h**2/(2*m) * Laplacian

HHarmonic = Hamiltionian + MHarmonic

HAntiHarm = Hamiltionian + MAntiHarm

#A Print for testing
#print(HHarmonic)

#Calculates the Eigenvalues and eigenvector for the Harmonic
HarmonicVals, HarmonicVec = np.linalg.eig(HHarmonic)
#Calculates the Eigenvalues and eigenvector for the Antiharmonic
AntiHarmVals, AntiHarmVec = np.linalg.eig(HAntiHarm)
#Gets the indecies that would sort the Eigenvalues
idxHm = np.argsort(HarmonicVals)
idxAH = np.argsort(AntiHarmVals)
#Sorts the Eigenvalues
Sorted_HarmEig = HarmonicVals[idxHm]
Sorted_AntiEig = AntiHarmVals[idxAH]
#Sorts the Eigenvectors by collumn because if sorted by rows it won't work and hours of time will be wasted
Sorted_AntiEV = AntiHarmVec[:, idxAH]
Sorted_HarmEV = HarmonicVec[:, idxHm]

#A Color array for the plot to look pretty
Colors = ['pink','maroon','red','orange','gold','olive','green','cyan','blue','purple']
#A Label array so I know which color is which
Labels = ['First','Second','Third','Fourth','Fifth','Sixth','Seventh','Eighth','Ninth','Tenth']


#a function to graph everything out
def graph(EigV,Lable: str):
    ''' Graphs a wave fucntion where:
    EigV = the Eignenvector of the wave function
    Lable = The name of the created graph
    '''
    for a in range(10):
        #Picks one wave function out of the array, of course it has to be a collumn, not a row
        Wave_Func = EigV[:, a]
        #Takes aaid function's energy level for later

        #Creates the plot for the wave functions
        plt.plot(Real_Space_Grid,Wave_Func, linestyle = '-', color= Colors[a], label= str(Labels[a]))
        plt.xlabel('x')
        plt.ylabel('\u03C8(x)')
        plt.title(Lable)
    #Sets the legend to the upper right so it's out of the way
    plt.legend().set_loc('upper right')
    plt.show()

def graphlines(Eig,Lable: str):
    #Creates the plot for the energy levels, limiting the X and Y so the graph looks decent
    plt.xlim(-L/2, L/2)
    #Automatically decides where to put the y limiter based off of the energy levles
    plt.ylim(0, Eig[10])
    plt.xlabel('x')
    plt.ylabel('E')
    plt.title(Lable)
    
    #A for loop to put all 10 energy levels on the same graph
    for a in range(10):
        plt.axhline(Eig[a], linestyle = '-', color= Colors[a], label= 'n = ' + str(a + 1))
    plt.legend().set_loc('upper right')
    plt.show()
    

#Calling the function I just wrote to graph the harmonic and antiharmonic wave functions
graph(Sorted_HarmEV,'The first 10 harmonic quantom wave functions')
graph(Sorted_AntiEV, 'The first 10 anti-harmonic quantom wave functions')

#Callling the other function I wrote to graph the energy levels
graphlines(Sorted_HarmEig,'The first ten energy levels of the harmonic quantom wave functions')
graphlines(Sorted_AntiEig,'The first ten energy levels of the antiharmonic quantom wave functions')
