#importing mat plot lib and pandas since I want to try to do this from a CSV
#The main reason I'm importing a CSV and not the functions themselves is that I don't want to run that hard drive frier again!
import matplotlib.pyplot as plt
import pandas as pd

#Setting the directory to our file so I can find the rest of my code, and so I save the graphs to the disired spot
#Dang I must have been really frustrated when microsoft's stupid CoPilot AI told me how to do this
import os
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)

#Extracts the CSV file
df = pd.read_csv('Heat_Capcaity_Calculations.csv')

Temperature = df['Temperature']
Internal_Energy = df['Internal Energy'] 
Heat_Capacity = df['Heat Capacity']

#prepares to plot two subplots, the way the matplotlib docs told me to
fig, ax1 = plt.subplots()
# Plot units for U on the left y-axis
ax1.plot(Temperature, Internal_Energy, linestyle = '-', color= 'red')
ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Internal Energy U (eV)', color='red')
ax1.tick_params(axis='y', labelcolor='red')

#Create a second y-axis for the Heat capacity
ax2 = ax1.twinx()
ax2.plot(Temperature, Heat_Capacity, linestyle = '-', color= 'blue')
ax2.set_ylabel('Heat Capacity C_v (eV/K)', color='blue')
ax2.tick_params(axis='y', labelcolor='blue')

#Sets the legend to the upper left so it's out of the way
plt.legend().set_loc('upper left')
#Titles the graph Work done to gas when expanded from 0.1m^3
plt.title("Cv and U as a function of temperature")

#Saves the graph and shows it to ensure accuracy, the arguments bbox_inches and pad_inches ensures that both axis titles are legible
#You didn't ask for this but I'm giving it anyway to be safe
plt.savefig('Cv_and_U_graph.png', bbox_inches= 'tight', pad_inches= 0.2)
plt.show()