#Imports all the functions I can't code as well as a few constants
import numpy as np
import scipy.constants as con
import pandas as pd
import matplotlib.pyplot as plt

#A function I found online that allows the import of IPYNB files
import import_ipynb

#imports the rules for json files
import json

#Pulls the notebook from where I placed it so that we can call upon it
notebook = import_ipynb.find_notebook('Data/lecture-07-regression')

#reads the file as a json file and puts each cell into an array cells
file = open(notebook).read()
cells = json.loads(file)
cells = cells["cells"]

codeCells = []

#Sorts through the cells and adds any code cells to an array codeCells
for i in cells:
    if i["cell_type"] == "code":
        code = ""
        for j in i["source"]:
            code += j
        codeCells.append(code)

#Executes the code cells globally so I can call the functions within from anywhere in the python file
def executeCell(cell_code):
    exec(cell_code, globals())

#These are the code cells corosponding to ols_slope, ols_intercept, and ols respectivly
executeCell(codeCells[2])
executeCell(codeCells[3])
executeCell(codeCells[4])

#This is a terrable way to do this, the code is not even really imported as much as it is pulled straight from the heart of the .ipynb file and shoved into the globals of this file
#I spent 3 hours reaserching nothing but this and this was the only solution I came accross that remotly works
#Here is my source https://www.youtube.com/watch?v=-2Q5ikJbZVU

#Defines deltaSV
DeltaSV = con.R * 10.5

#Imports our csv as a panda dataframe
df = pd.read_csv('Data/trouton.csv')

#Converts Kcals to Joules
df['H_v (kcal/mol)'] = df['H_v (kcal/mol)'].apply(lambda x: x * con.calorie * 1000)
df.rename(columns={'H_v (kcal/mol)': 'H_v (joules/mol)'}, inplace=True)

#A Print Function for testing
print(df)

#Calculates the slope and interecept using the defined functions
slope, intercept = ols(df['T_B (K)'],df['H_v (joules/mol)'])

#A list of colors for the datapoints
Colors = ['blue','red','purple','gold']

#Runns a for loop for every value in H_V joules/mol
for a in range(len(df['H_v (joules/mol)'])):
    
    #Gets the x and y from the pd datafarme
    x = df.loc[a, 'T_B (K)']
    y = df.loc[a, 'H_v (joules/mol)']
    #Decides which color the point should be based of the Class
    Color = str
    Class = df.loc[a, 'Class']
    if Class == 'Perfect liquids':
        Color = 'blue'
    elif Class == 'Imperfect liquids':
        Color = 'red'
    elif Class == 'Liquids subject to quantum effects':
        Color = 'purple'
    else:
        Color = 'Gold'

    #Plots each point, with a color according to its class
    plt.plot(x,y, marker = '.', color= Color)
    plt.xlabel('T_B (K)')
    plt.ylabel('H_v (jules/mol)')
    #Sets the legend to the upper right so it's out of the way
plt.legend().set_loc('upper right')
plt.show()
