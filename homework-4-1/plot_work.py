#importing mat plot lib, I shouldn't need anything else
import matplotlib.pyplot as plt

#Setting the directory to our file so I can find the rest of my code, and so I save the graphs to the disired spot
#Dang I must have been really frustrated when microsoft's stupid CoPilot AI told me how to do this
import os
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)

#Importing the files I wrote
import compute_work_adiabatic as adb
import compute_work_isothermal as iso

#Defining my values based off the values in the original files
R = adb.R
n = iso.n
T = adb.T
V_init = iso.V_init
ADBindx = adb.ADBindx #This one should be the only one that matters since only the adb file has ADBIndex, the rest are just copy pasted
Vf = iso.Vf

#Calculates the work done at every final volume both adiabaticly and isothermicly
work_done_iso = []
work_done_adb = []
#Uses a for loop to ensure theres no wonkyness and it puts one VF in at a time
for v in range(len(Vf)):
    work_done_iso.append(iso.work_calculator_iso(V_init,Vf[v],n,R,T,1000))
    work_done_adb.append(adb.work_calculator_adb(V_init,Vf[v],n,R,T,1000, ADBindx))

#plots the work done with respect to Vf
plt.plot(Vf, work_done_iso, linestyle = '-', color= 'blue', label= "Work done isothermicly")
plt.plot(Vf, work_done_adb, linestyle = '-', color= 'red', label= "Work done adiabaticly")

#lables the X and Y axis
plt.xlabel('Final Volume (m^3)')
plt.ylabel('Work done (Joules)')

#Sets the legend to the upper left so it's out of the way
plt.legend().set_loc('upper left')
#Titles the graph Work done to gas when expanded from 0.1m^3
plt.title("Work done to ideal gas when expanded from 0.1m^3")

#Saves the graph and shows it to ensure accuracy, the arguments bbox_inches and pad_inches ensures that both axis titles are legible
#You didn't ask for this but I'm giving it anyway to be safe
plt.savefig('work-graph.png', bbox_inches= 'tight', pad_inches= 0.2)
plt.show()