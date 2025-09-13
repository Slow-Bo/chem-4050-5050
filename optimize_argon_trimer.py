#Importing all the cool functions I couldn't code myself
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import pandas as pd

#a function to calculate bond length the input is 2 atom's positions
def Bond_Length(Atom1 = [0, 0, 0], Atom2 = [0, 0, 0]):
    #distance starts at 0
    distance = 0
    # The for loop adds the square difference in position in all dimentions together
    for position in range(min(len(Atom1),len(Atom2))):
        diff = (Atom1[position] - Atom2[position])**2
        distance += diff
    return np.sqrt(distance)
    """
    Removed check for valid bond lengths, Argon is not one for close relationships

    if 0 < np.sqrt(distance) < 1.5:
        # returns the square root of the diffrences squared as per pythagoras' theorm
        return np.sqrt(distance)
    else:
        pass
    """

#A function for subtracting one vector from aonther
def VecSubtract(vector01, vector02):
    #creates a new vector using the first 3 values of two other vectors
    vector03 = [vector01[value] - vector02[value] for value in range(min(len(vector01), len(vector02)))]
    return vector03
#A function for computing the bond angle
def Bond_Angle(Atom1 = [0, 0, 0], Atom2 = [0, 0, 0], Atom3 = [0, 0, 0]):
    #These values will be used in calculations
    Numerator = 0
    #Calls my subtracting function because I am too lazy to write the same thing twice
    Vector_1 = VecSubtract(Atom3, Atom1)
    Vector_2 = VecSubtract(Atom2, Atom1)
    #Figures out the absolute value of each vector by calling the bond length function
    AbsV1 = Bond_Length(Atom1, Atom3)
    AbsV2 = Bond_Length(Atom2, Atom1)
    #Checks for invalid bond lengths
    if AbsV1 == None or AbsV2 == None:
        pass
    else:
    #Figures out the dot product of the vector
        for value in range(min(len(Vector_1),len(Vector_2))):
            Numerator += np.multiply(Vector_1[value], Vector_2[value])
    #Plugs all priorly calculated numbers into the equation to calculate bond angle in degrees
        Result = np.degrees(np.arccos(Numerator/(AbsV1*AbsV2)))

        #I kept getting 0 angles so I added this redundancy to stop it
        if Result < 0.001:
            return
        #Prints what type of angle the bond angle is
        elif Result < 90:
         print("The angle is acute")
        elif Result == 90:
            print("The angle is right")
        else:
            print("The angle is obtuse")

        return (Result)

    #A simple function that takes an input of r (and potentially epsilon and sigma) and outputs the atomic energy
def lennard_jones(r, epsilon=0.01, sigma=3.4):
    
    Atomic_Energy = 4*epsilon*((sigma/r)**12-(sigma/r)**6)
    return Atomic_Energy



#This is the function to be optimized added bools to ensure bond angle and atom location will only be returned when wanted
def The_Atom_Space(locals = [0,0,0], calcbondangle = False, giveatomlocation = False):
    """
    Computes the Potental Energy, Bond Lengths, Bond Angles, and atom location in 3D space
    Parameters:

    locals[x1, x2, y2] where:

    x1 (Float): The X position of atom 2
    x2 (Float): The X position of atom 3
    y3 (Float): The Y position of atom 3

    calcbondangle (Bolean): True if you want the Bond Angles and Lengths returned, otherwise false
    giveatomlocation (Bolean): True if you want the Atom Location returned, otherwise false

    The function will ONLY return one piece of information at a time so choose carefully what you want

    Returns:
    Float: The total potental energy and any asked for information
    """
    
    #Locks one atom's position at 0,0,0. We won't be using the since three points will always form a plane, but I need three inputs to run my function
    #Directory so I have a number I can use to prevent re-running the same length. That bug was a non issue in HW 1-2 but it would ruin everything here.
    Atoms = {
        1 : [0,0,0], 
        2 : [locals[0],0,0], 
        3 : [locals[1],locals[2],0]}
    #Setting up our function variables for later
    energies = []
    bond_lengths = []
    bond_angles = []
    atom_locations = []
    for value1, atom1 in Atoms.items():
        for value2, atom2 in Atoms.items():
            #If statement ensures we have 2 different atoms and prevents the same length from being run twice. 
            if value1 < value2:
                #Calls the bond length function and adds the result to the empty array
                BL = Bond_Length(atom1, atom2)
                #Stops the appending of invalid values
                bond_lengths.append(BL)
    #Calculates the energy for each argon distance in the trimer
    for r in bond_lengths:
        V = lennard_jones(r)
        energies.append(V)

    #Sums together all the energies
    total_potental = sum(energies)
    
    #Calls for bond_angle calculations only if asked
    if calcbondangle == True:
        #Nested for loop allows us to grab three atoms from molecules instead of one
        Ammendment: bool = True
        for atom1, value1 in Atoms.items():
            for atom2, value2 in Atoms.items():
                for atom3, value3 in Atoms.items():
                    #If statement ensures we have 3 different atoms
                    if atom1 != atom2 and atom2 != atom3 and atom1 != atom3:
                     #I kept getting duplicate Angles, this method ensures that I won't by only appending every other values
                        if Ammendment == True:
                            #Calls the bond angle function and adds the result to the empty array
                            BA = Bond_Angle(value1, value2, value3)
                            Ammendment = False
                            #Stops the appending of invalid values
                            if BA != None:
                                bond_angles.append(BA)
                        else:
                            Ammendment = True
    #Packages Atom Locations if asked
    if giveatomlocation == True:
        for atom, value in Atoms.items():
            atom_locations.append(value)
    
    #Gives requested information
    if calcbondangle == False and giveatomlocation == False:
        return total_potental
    else:
        return bond_angles, bond_lengths if calcbondangle == True else atom_locations
    #In retrospect I could have just written multiple functions but this is way cooler

#This is the function that does the optimization
distance_between_three = opt.minimize(
    #This is the function
    fun=The_Atom_Space,  
    #This is the initial guess 3.6 was the optimal distance in problem 1
    x0=[3.6,1.8,3.12],                    
    #The arguments of the function scipy is not allowed to change
    args=(False, False),
    #This is the method, I don't remeber the difference between them so I'm using the first one taught in class
    method="Nelder-Mead",    
    #This is the tolernace because computers have yet to acheve perfection
    tol=1e-6                 
)

#A print for testing
#print(The_Atom_Space(distance_between_three["x"], False, True))

#The Print function for completing the assignment
def completing_the_mission():
    print("Optimal Atom Locations in Angstroms")
    #For some reason when appending the atom space it gets a garbage array that I can't figure out how to remove in slot 1, the [1] lets me ignore that quirk
    print(The_Atom_Space(distance_between_three["x"], False, True)[1])
    print("The bond lengths and angles between atoms")
    print(The_Atom_Space(distance_between_three["x"], True, False))
    print("The triangle appears to be equilateral, with the same angles and distances between atoms")

#The function that makes the file
def file_creator():
    #Gets the atom locations, using the [1] to ignore the garbage data
    locals = The_Atom_Space(distance_between_three["x"], False, True)[1]
    
    #Pandas was not taking my sig figs seriously, this secton of code converts my array of arrays into a new array of arrays where all my values are strings
    #This is not a good way to do this, but it was the only way I could figure out on my own
    locals2 = []
    for i in locals:
        blank_array = []
        for n in i:
            n = "{:.6f}".format(n)
            blank_array.append(n)
        locals2.append(blank_array)
    
    #Converts the stringed data to a Pandas Dataframe
    df = pd.DataFrame(locals2, index=['Ar1', 'Ar2', "Ar3"],columns=['','',''])
    
    #adds the lables nessasary for the XYZ file format
    df.loc['Argon Trimer'] = ['', '', '']
    df.loc['3'] = ['','','']
    #rearanges dataframe so that the lable is on top
    df = df.reindex(['3','Argon Trimer','Ar1','Ar2','Ar3'])
    #Makes the file
    df.to_csv('./homework-2-1/ArgonTrimer.xyz', sep=' ',index=True, columns= None, header= False)
    #Prints the Panda Dataframe so I can see I put it together right
    print(df)
    

#Runs all my functions
file_creator()
completing_the_mission()