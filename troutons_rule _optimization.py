#Imports all the functions I can't code as well as a few constants
import numpy as np
import scipy.constants as con
import scipy.optimize as opt
from scipy.stats import t
import pandas as pd
import matplotlib.pyplot as plt

#Added this import so file finding works better it turns out that depending on how I boot up VS code it changes the directory to whatever it feels like,
#This holds a gun to VS Code's head and forces it to take the directory I want, though admitildy microsoft's stupid CoPilot AI immediatly told me the answer when I tried to look it up
import os
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)

#Removed all the garbage code held together with hopes and dreams

#The function to be optimized by SciPy, contains the H_v value (y), the T_B value (x), the slope (a), and the intercept (b)
def Square_Error(a= float,b= float,x= [],y= []):
    Error = sum((y-(x*a+b))**2)
    return Error

#This function will actualy call SciPy optimize and retern the optimal slope
#Contains the x values (x), the y values (y), the inital slope guess (a0), and the inital intercept guess (b0)
def Slope_Optimizer(X,Y,a0,b0):
    
    #This is the way you are supposed to pass arguments through functions with optimize, I should have read the SciPy docs last week
    fun=lambda x=[] : Square_Error(x[0], x[1], X, Y)
    
    #Runs the optimize function and retrives the optimal slope and intercept
    Rawdata = opt.minimize(
        fun, #The function to be optimized
        x0= [a0,b0], #The inital guesses
        method= 'Nelder-Mead' #The only method I know for real
        )["x"] #Only returns the optimal value
    return Rawdata[0], Rawdata[1] #returns the slope and intercept
#Calculates the confidence intervals
def Confidence(x_values, y_values, slope, intercept, confidence = 0.95):
    """
    Computes the confidence interval for a set of values and their predicted values.

    Parameters:
    x_values (array): The points of data for the independant variable.
    y_values (array): The points of data for the dependant variable.
    slope (float): The calculated slope of the regression line.
    intercept (float): The calculated intercept of the regression line.
    confidence (float > 1): The desired confidence interval.

    Returns:
    float: The confidance interval for the slope
    float: The confidance interval for the intercept
    """
    #Creates our variables we will need for later
    Freedom = (len(y_values)-2)
    alpha = 1 - confidence
    critical_t_value = t.ppf(1 - alpha/2, Freedom)
    mean = np.mean(x_values)
    #calculates Standard Square Error
    SSR = sum((y_values - (x_values*slope+intercept))**2)
    #I remember you, the math, Variance
    VAR = SSR/Freedom
    
    #calculates the standard Error for both slope and intercept
    SESlope = np.sqrt(VAR/(sum((x_values-mean)**2)))
    SEIntercept = np.sqrt(VAR*(1/len(y_values)+(mean**2/sum((x_values-mean)**2))))

    #Retuns the confindence interval for both values
    return SESlope * critical_t_value, SEIntercept * critical_t_value


#Defines deltaSV
DeltaSV = con.R * 10.5

#Imports our csv as a panda dataframe
df = pd.read_csv('Data/trouton.csv')

#Converts Kcals to Joules
df['H_v (kcal/mol)'] = df['H_v (kcal/mol)'].apply(lambda x: x * con.calorie * 1000)
df.rename(columns={'H_v (kcal/mol)': 'H_v (joules/mol)'}, inplace=True)

T_BValues = df['T_B (K)']
H_VValues = df['H_v (joules/mol)']



#Calculates the slope and interecept using the optimizer function with 10.5 * R and 0 as the inital guesses
slope, intercept = Slope_Optimizer(T_BValues,H_VValues,con.R * 10.5, 0)
#Calculates the confidence intervals for the slope and intercept at 95% Confidence interval
CIS, CII = Confidence(T_BValues,H_VValues,slope,intercept,0.95)

#A print for testing
print(CIS, CII)


#Creates 100 points from 1 to 2500
X = np.linspace(1,2500,100)
#Graphs the calculated slope intercept line
plt.plot(X,X * slope + intercept, linestyle = '-', color= 'black', label= f"H_v = {slope:.3f} +/- {CIS:.3f}(J/mol) + {intercept/1000:.3f} +/- {CII/1000:.3f}(kJ/mol)")

#An array for the for loop to keep track of what classes have been labled
LabledClasses = []
#Runs a for loop for every value in H_V joules/mol
for a in range(len(H_VValues)):
    
    #Gets the x and y from the pd datafarme
    x = T_BValues[a]
    y = H_VValues[a]
    
    #Decides which color the point should be based of the Class
    Color = str
    Marker = str
    Class = df.loc[a, 'Class']
    if Class == 'Perfect liquids':
        Color = 'blue'
        Marker = 'Perfect liquids'
    elif Class == 'Imperfect liquids':
        Color = 'red'
        Marker = 'Imperfect liquids'
    elif Class == 'Liquids subject to quantum effects':
        Color = 'purple'
        Marker = 'Liquids subject to quantum effects'
    else:
        Color = 'Gold'
        Marker = 'Metals'
    #Checks to see if the class has already been labled
    if Marker in LabledClasses:
        #Plots each point, with a color according to its class, but does not add a label
        plt.plot(x,y, marker = 'o', color= Color)
    else:
        #Plots each point, with a color according to its class
        plt.plot(x,y, marker = 'o', color= Color, label= Marker)
        #Adds the Marker to Labeled Classes so classes aren't labled twice
        LabledClasses.append(Marker)
    plt.xlabel('T_B (K)')
    plt.ylabel('H_v (joules/mol)')

#Sets the legend to the upper left so it's out of the way
plt.legend().set_loc('upper left')
#Titles the graph Trouton's Rule
plt.title("Trouton’s Rule Optimization")

#Saves the graph and shows it to ensure accuracy, the arguments bbox_inches and pad_inches ensures that both axis titles are legible
plt.savefig('homework-3-2/Trouton.png', bbox_inches= 'tight', pad_inches= 0.2)
plt.show()

#comments on optimization versus the linear regression equation
'''The difference between the optimization method and the linear regression equation graphs is that there is no difference. They're the same graph. I think the main reson for this
is that the linear regression equation was derived to find the least squares with just one equation. This means that the least squares method is more optimal for your computer
since it just needs to plug in the numbers once to get the disired value. The optimization method simply brute forces the equation untill it reaches a local/global minima.
However the main advantage to optimization is that you can just optimize any equation, but if you wanted to find the minimum value of any equation but the least squares value
than you would need to derive that equation, assuming it can be derived, and find where the derivitive is 0, then determine if it is the global minimum. That is a lot of work.
In conclusion I think for a linear equation just using the linear regression equation is better because it exists and already finds the least squares line with very little effort,
but if you had an equation you didn't have the equivlent of the linear regression equation for optimization would be optimal.'''