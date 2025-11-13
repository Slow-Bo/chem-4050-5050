#Imports all the functions I can't code as well as a few constants
import numpy as np

import matplotlib.pyplot as plt

#The variables to run throught all my functions
sizeX = 4
sizeY = 4

#Some filler parameters for testing because of how I coded this any dictionary I feed in as parameters, MUST have this format
Filler_params ={
    'eH': -0.1,
    'eN': -0.1,
    'eNN': 0.05,
    'eHH': 0.05,
    'eNH': -0.05,
    'mu_N': -0.1,
    'mu_H': -0.1,
    'T': 200
}

#I couldn't figure out the function to bind a value to a certain range so I made my own
def Modulalo(X,List_Size):
    #Since python arrays start at 0 the max value for X is one less than list size
    Xmax = List_Size - 1
    if X > Xmax:
        #If Xmax is smaller than X it loops back to 0
        return 0
    elif X < 0:
        #If X is negitive it gives the maximum value for X
        return Xmax
    else:
        #If nither of the above are true X is valid so it just returns X
        return X
    #Function to determine the neighbors of each cell
def mister_rogers(sizeX, sizeY):
    #Prepares a dictionary for later use
    neighbors = {}
    #Evaluates a grid of size X * Y
    for X in range(sizeX):
        for Y in range(sizeY):
            #Places four values in each spot in the dictionary for each negihbor of the index cell
            neighbor = [
            (Modulalo((X-1),sizeX), Y),
            (Modulalo((X+1),sizeX), Y),
            (X, Modulalo((Y-1),sizeY)),
            (X, Modulalo((Y+1),sizeY))
            ]
            neighbors[(X,Y)] = neighbor
    return neighbors

#Creates a grid of size X*Y
def create_grid(sizeX, sizeY):
    raw_grid = np.zeros((sizeX, sizeY))
    #Calls upon the previous function to determine the neighbors of each cell,.
    defined_grid = mister_rogers(sizeX, sizeY)
    #Returns the grid with all neighbors assigned
    return raw_grid, defined_grid

#The function that calculates the energy of each cell
def energy_calculator(grid, neighbors, parameters):
    #Starting energy is 0
    Energy = 0
    #Gets the energy from each of the parameters in parameters
    Nrg = parameters.get('eN')
    Hrg = parameters.get('eH')
    NNrg = parameters.get('eNN')
    HHrg = parameters.get('eHH')
    HNrg = parameters.get('eNH')
    #Evaluates every cell in every row of the grid.
    for i, row in enumerate(grid):
        for j, element in enumerate(row):
            #If the cell is 0 it can be ignored since it has nothing to offer
            if element == 0:
                pass
            #BUT I HAVE TO DO THIS NONSENCE BECAUSE NO STRINGS IN THE INT MATRIX I LOVE PYTHON
            elif element == '0':
                pass
            elif element =='0.0':
                pass
            #If the cell is nitrogen:
            elif element == 'N':
                #creates an array for the neghibors of the current cell
                Adj = neighbors.get((i,j))
                #adds the energy of Nitrogen binding to the grid
                Energy += Nrg
                for position in Adj:
                    X = position[0]
                    Y = position[1]
                    #If a neighbor is 0 it is ignored since it's lazy
                    if grid[X][Y] == 0:
                        pass
                    #Depending on if the neghbor is H or N adds the corosponding energy
                    elif grid[X][Y] == 'N':
                        Energy += NNrg
                    elif grid[X][Y] == 'H':
                        Energy += HNrg
            #If the cell is hydrogen:
            elif element == 'H':
                #creates an array for the neghibors of the current cell
                Adj = neighbors.get((i,j))
                #adds the energy of Hydrogen binding to the grid
                Energy += Hrg
                for position in Adj:
                    X = position[0]
                    Y = position[1]
                    #If a neighbor is 0 it is ignored since it's lazy
                    if grid[X][Y] == 0:
                        pass
                    #Depending on if the neghbor is H or N adds the corosponding energy
                    elif grid[X][Y] == 'N':
                        Energy += HNrg
                    elif grid[X][Y] == 'H':
                        Energy += HHrg
    #Gives back the energy
    return Energy

    
#Function to add a particle to the grid
def gonna_add_one(grid, neighbors, parameters, Empt_Sites, N_Sites, H_Sites):
    """
    Determens Whether and where to add a particle to a grid.

    Parameters:
    grid (matrix): The grid to have a particle added to it.
    neighbors (dictionary): A dictonary tied to the associated grid that contains the neighbors for each cell.
    parameters (dictionary): A dictonary of extremly particular form contaning the needed parameters for simulation.
    Empt_Sites (Integer): The number of empty sites in the grid
    N_Sites (Integer): The number of nitrogen sites in the grid
    H_Sites (Integer): The number of hydrogen sites in the grid

    Returns:
    Empt_Sites (Integer): The number of empty sites in the new grid
    N_Sites (Integer): The number of nitrogen sites in the new grid
    H_Sites (Integer): The number of hydrogen sites in the new grid
    grid (Matrix): The new grid, whether or not it's the same as the old grid
    """
    #Gets the needed parameters from the parameters
    N_Pot = parameters.get('mu_N')
    H_Pot = parameters.get('mu_H')
    #Calculates the energy of the current grid
    OldE = energy_calculator(grid, neighbors, parameters)
    beta = 1/(parameters.get('T'))
    #If there are no empty spots the function terminates immediatly
    if Empt_Sites == 0:
        return Empt_Sites, N_Sites, H_Sites, grid
    #If there are empty spots the function actually prepares to do something
    else:
        #Creates an array for valid sites
        Valid_Sites = []
        #Checks every cell in every row for a valid value
        for i, row in enumerate(grid):
            for j, element in enumerate(row):
                #since I had to make the 0 a string later for the function to work I have to acount for 0 being both a string and not a string addidng a cell to valid values for having any of the three
                if element == 0:
                    Valid_Sites.append((i, j))
                elif element == '0.0':
                    Valid_Sites.append((i, j))
                elif element == '0':
                    Valid_Sites.append((i, j))
        #Choses one of the valid sites at random and gets its X and Y value
        Chosen_Site = np.random.choice(len(Valid_Sites))
        Chosen_X = Valid_Sites[Chosen_Site][0]
        Chosen_Y = Valid_Sites[Chosen_Site][1]
        #Decided to add either a N or H with a 50/50 chance
        N_or_H = np.random.choice(2) > 0.5
        #When the function decides to add Nitrogen
        if N_or_H == True:
            #Converts all the values in the grid to a string because God has a cruel sence of humor and this is the only way it would work
            new_grid = [[str(ele) for ele in k] for k in grid]
            #Adds a Nitrogen to the chosen spot of the new grid
            new_grid[Chosen_X][Chosen_Y] = 'N'
            #Calculates the energy of the new grid
            NewE = energy_calculator(new_grid, neighbors, parameters)
            #Runs the metropolis algrithom and prepares a random value to compare it's result to
            Acc = min(1, (Empt_Sites - N_Sites) / (N_Sites + 1) * np.exp(-beta * ((NewE - OldE) - N_Pot)))
            if np.random.rand() <= Acc:
                #When the grid is accepted by metropolis it sets the old grid equal to it and changes the number of N and empty sites
                grid = new_grid
                N_Sites += 1
                Empt_Sites -= 1
        #When the function decides to add Hydrogen
        elif N_or_H == False:
            #Converts all the values in the grid to a string because God has a cruel sence of humor and this is the only way it would work
            new_grid = [[str(ele) for ele in k] for k in grid]
            #Adds a Hydrogen to the chosen spot of the new grid
            new_grid[Chosen_X][Chosen_Y] = 'H'
            #Calculates the energy of the new grid
            NewE = energy_calculator(new_grid, neighbors, parameters)
            #Runs the metropolis algrithom and prepares a random value to compare it's result to
            Acc = min(1, (Empt_Sites - H_Sites) / (H_Sites + 1) * np.exp(-beta * ((NewE - OldE) - H_Pot)))
            if np.random.rand() < Acc:
                #When the grid is accepted by metropolis it sets the old grid equal to it and changes the number of N and empty sites
                grid = new_grid
                H_Sites += 1
                Empt_Sites -= 1
        #Returns the new amount of emty sites, N, and H on the new grid
        return Empt_Sites, N_Sites, H_Sites, grid
        
#Function to remove things from the grid
def take_away_one(grid, neighbors, parameters, Empt_Sites, N_Sites, H_Sites):
    """
    Determens Whether and where to take a particle from a grid.

    Parameters:
    grid (matrix): The grid to have a particle added to it.
    neighbors (dictionary): A dictonary tied to the associated grid that contains the neighbors for each cell.
    parameters (dictionary): A dictonary of extremly particular form contaning the needed parameters for simulation.
    Empt_Sites (Integer): The number of empty sites in the grid
    N_Sites (Integer): The number of nitrogen sites in the grid
    H_Sites (Integer): The number of hydrogen sites in the grid

    Returns:
    Empt_Sites (Integer): The number of empty sites in the new grid
    N_Sites (Integer): The number of nitrogen sites in the new grid
    H_Sites (Integer): The number of hydrogen sites in the new grid
    grid (Matrix): The new grid, whether or not it's the same as the old grid
    """
    #Gets the needed parameters from the parameters
    N_Pot = parameters.get('mu_N')
    H_Pot = parameters.get('mu_H')
    #Calculates the energy of the current grid
    OldE = energy_calculator(grid, neighbors, parameters)
    beta = 1/(parameters.get('T'))
    #If there are no particles to remove the function immediatly terminates
    if N_Sites + H_Sites == 0:
        return Empt_Sites,N_Sites,H_Sites, grid
    #If the function actually decides to do something:
    else:
        #Creates an array for valid sites
        Valid_Sites = []
        #Checks every cell in every row for a valid value
        for i, row in enumerate(grid):
            for j, element in enumerate(row):
                #Considers a site valid if it has either an N or an H, I just realized I could have programed the other one to accept any value that wasn't N or H and it would have been easier
                if element == 'N':
                    Valid_Sites.append((i, j))
                elif element == 'H':
                    Valid_Sites.append((i, j))
        #Choses one of the valid sites at random and gets its X and Y
        Chosen_Site = np.random.choice(len(Valid_Sites))
        Chosen_X = Valid_Sites[Chosen_Site][0]
        Chosen_Y = Valid_Sites[Chosen_Site][1]
        #Creates a new grid that is the same as the old grid except with the chosen value set to 0
        #Notice how I can use this grid without making everything a strig? Isn't life wonderful when you can do that?
        new_grid = [k[:] for k in grid]
        new_grid[Chosen_X][Chosen_Y] = 0
        #Calculates the energy of the new grid
        NewE = energy_calculator(grid, neighbors, parameters)
        #If the original                grid, contained a nitrogen:
        if grid[Chosen_X][Chosen_Y] == 'N':
            #A failsafe for when the denominator of the first part of the metropolis algrythem == 0
            if Empt_Sites - N_Sites + 1 == 0:
                Denom = 1E-30
            else:
            #If it dosn't the function just calculates the algryhtem and determins if it should accept the new grid
                Denom = Empt_Sites - N_Sites + 1
            Acc = min(1, N_Sites/Denom * np.exp(-beta * ((NewE - OldE) + N_Pot)))
            if np.random.rand() <= Acc:
                #If it does sets the old grid to the new grid and changes the values of N_Sites and Empty_Sites to match
                grid = new_grid
                N_Sites += -1
                Empt_Sites += 1
        #If the original                grid, contained a hydrogen:
        elif grid[Chosen_X][Chosen_Y] == 'H':
            #A failsafe for when the denominator of the first part of the metropolis algrythem == 0
            if Empt_Sites - H_Sites + 1 == 0:
                Denom = 1E-30
            #If it dosn't the function just calculates the algryhtem and determins if it should accept the new grid
            else:
                Denom = Empt_Sites - H_Sites + 1
            Acc = min(1, H_Sites/Denom * np.exp(-beta * ((NewE - OldE) + H_Pot)))
            if np.random.rand() <= Acc:
                #If it does sets the old grid to the new grid and changes the values of N_Sites and Empty_Sites to match
                grid = new_grid
                H_Sites += -1
                Empt_Sites += 1
        #Returns the new grid and the new values of N and H sites
        return Empt_Sites,N_Sites,H_Sites, grid

#Calculates the grand connonical where Steps is the number of steps, sizeX and Y are the dimentions of the grid and parameters are the rules of the grid.
def grand_cannonical(Steps, sizeX, sizeY, parameters):

    #Creates the grid by running the function I coded 200 years ago, hopefully it still works
    grid, neighbors = create_grid(sizeX, sizeY)
    #Calculates the total sites based of the dimentions of the grid and makes them all empty
    Total_Sites = sizeX * sizeY
    Empt_Sites = Total_Sites
    N_Sites = 0
    H_Sites = 0
    #Creates empty arrays to be filled by the grand cannonical
    NCoverage = []
    HCoverage = []
    
    #Every step the function runs
    for i in range(Steps):
        #Choses whether to run the add or subtract function
        #There is no way on God's green earth anyone but me undestands the refrence of the names of these two functions
        add_or_subtract = np.random.choice(2) > 0.5
        if add_or_subtract == True:
            #If the function decides to add feeds all the values into the addition function
            Empt_Sites, N_Sites, H_Sites, grid = gonna_add_one(grid, neighbors, parameters, Empt_Sites, N_Sites, H_Sites)
            
        elif add_or_subtract == False:
            #If the function decides to add feeds all the values into the subtraction function
            Empt_Sites, N_Sites, H_Sites, grid = take_away_one(grid, neighbors, parameters, Empt_Sites, N_Sites, H_Sites)
        
        #Figures out what the coverage of the new grid is and adds it to the Metropolis array
        NCoverage.append(N_Sites/Total_Sites)
        HCoverage.append(H_Sites/Total_Sites)
    #Returns the result of the grand cannonical
    return grid, NCoverage, HCoverage

#Creating a blank grid for testing
grid, grid_data = create_grid(sizeX,sizeY)
#print(grid_data.get((0,1)))

#Creating a not so blank grid for further testing
Test_grid = [['N',0,0,0],
             ['H',0,0,0],
             ['H',0,0,0],
             ['N',0,0,0]
]
#A function that turns all those annoying stings into integers so my code will actually accept them
def string_ripper(lattice,string1:str,string2:str):
    #Creates a matrix of the same dimensions as the original, but with only zeros
    intlattice = np.zeros_like(lattice)
    #Checks every cell in every row for a valid value
    for i, row in enumerate(lattice):
        for j, element in enumerate(row):
            #Checks if a given cell has one of the two valid strings before swaping them with an integer
            if element == string1:
                intlattice[i][j] = 1
            elif element == string2:
                intlattice[i][j] = 2
            else:
                intlattice[i][j] = 0
    #returns a metrix of only integers so python will shut up and accept my freaking matrix
    return intlattice

#A function to plot the latice as a grid
def plot_lattice(lattice, ax, title, particle1:str, particle2:str):

    #Gets the X and Y size of the lattice
    gridsizex = len(lattice[0])
    gridsizey = len(lattice)

    # Create gridlines along the x and y
    ax.set_xticks(np.arange(0, gridsizex + 1, 1))
    ax.set_yticks(np.arange(0, gridsizey + 1, 1))
    #Sets Matplotlib to grid mode
    ax.grid(True)

    # Set the limits of the plot
    ax.set_xlim(0, gridsizex)
    ax.set_ylim(0, gridsizey)

    #Graphs the location of each adsorbant
    for x in range(gridsizex):
        for y in range(gridsizey):
            if lattice[x][y] == 1:
                ax.plot(x + 0.5, y + 0.5, 'o', color='blue', markersize=15, zorder=10)
            elif lattice[x][y] == 2:
                ax.plot(x + 0.5, y + 0.5, 'o', color='red', markersize=15, zorder=10)

    # Add labels for N and H, I found where to put them through trial and error
    ax.text(2, -0.5, title, ha='center')
    ax.text(1, 5.2, 'Adsorbed ' + particle1, color='blue', ha='center')
    ax.text(3.5, 5.2, 'Adsorbed ' + particle2, color='red', ha='center')
    
    # Hide ticks
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    return ax

#Test running the functions, there were more tests sewn throughout, but I got rid of them so I didn't have 1,000,000,000,000 prints
Grid, Junk1, Junk2 = grand_cannonical(100000,sizeX,sizeY,Filler_params) 
print(Grid)

#The arrays containing all the test values in the order eH, eN, eHH, eNN, eNH, muN
Ideal_Mixture = [-0.1,-0.1,0,0,0,-0.1]
Repulsive_Interactions = [-0.1,-0.1,0.05,0.05,0.05,-0.1]
Attractive_Interactions = [-0.1,-0.1,-0.05,-0.05,-0.05,-0.1]
Immisible_HN = [-0.1,-0.1,-0.05,-0.05,0.05,-0.1]
Like_Disolves_Unlike = [-0.1,-0.1,0.05,0.05,-0.05,-0.1]
#Puts all the arrays into a bigger array so I cn run them all at once
The_Big_Kahuna = [Ideal_Mixture,Repulsive_Interactions,Attractive_Interactions,Immisible_HN,Like_Disolves_Unlike]
#These are the basic values we will be testing, they should remain the same in all simulations
n_steps = 10000
mus_H = np.linspace(-0.2, 0, 7)
Ts = np.linspace(1, 201, 7)
#Runs through each array and prepares it for the Grand cannonical
for array in The_Big_Kahuna:
    params = []
    for mu_H in mus_H:
        for T in Ts:
            params.append({
                #Turns our arrays into a dictionary that can be fed into the grand cannoncial, I know it's a bit of a magic numbers solution, which is bad, but it does work.
                'eH': array[0],
                'eN': array[1],
                'eHH': array[2],
                'eNN': array[3],
                'eNH': array[4],
                'mu_H': mu_H,
                'mu_N': array[5],
                'T': T  
                })  #Temperature (in units of k)
    # Run the simulation 5 * 5 grids to be fancy
    sizex = 5
    sizey = 5
    np.random.seed(42)
    #Creates a matrix of matricies of the same size as my lattices
    final_lattice = np.zeros((len(mus_H), len(Ts), sizex, sizey))
    #creates two matricies of the size of the two array values we're changing
    mean_coverage_H = np.zeros((len(mus_H), len(Ts)))
    mean_coverage_N = np.zeros((len(mus_H), len(Ts)))
    for i, param in enumerate(params):
        #runs the grand cannonical and retrives it's results
        lattice, coverage_N, coverage_H = grand_cannonical(n_steps,sizex,sizey,param)
        #Fixed latice is the same as the old lattice but only contains integers so that it dosn't make python go into cardiac arrest when I shove it into the final lattice matrix
        Fixed_lattice = string_ripper(lattice,'H','N')
        final_lattice[i // len(Ts), i % len(Ts)] = Fixed_lattice
        mean_coverage_H[i // len(Ts), i % len(Ts)] = np.mean(coverage_H[-1000:])
        mean_coverage_N[i // len(Ts), i % len(Ts)] = np.mean(coverage_N[-1000:])
    
    # Plot the T-mu_A phase diagram
    fig, axs = plt.subplot_mosaic([[0, 1, 2],[3, 4, 5]], figsize=(6.5, 4.5))

    # Mean coverage of A
    axs[0].pcolormesh(mus_H, Ts, mean_coverage_H.T, cmap='viridis', vmin=0, vmax=1)
    axs[0].set_title(r'$\langle \theta_A \rangle$')
    axs[0].set_xlabel(r'$\mu_A$')
    axs[0].set_ylabel(r'$T$')

    # Mean coverage of B
    axs[1].pcolormesh(mus_H, Ts, mean_coverage_N.T, cmap='viridis', vmin=0, vmax=1)
    axs[1].set_title(r'$\langle \theta_B \rangle$')
    axs[1].set_xlabel(r'$\mu_A$')
    axs[1].set_yticks([])

    # Mean total coverage
    cax = axs[2].pcolormesh(mus_H, Ts, mean_coverage_H.T + mean_coverage_N.T, cmap='viridis', vmin=0, vmax=1)
    axs[2].set_title(r'$\langle \theta_A + \theta_B \rangle$')
    axs[2].set_xlabel(r'$\mu_A$')
    axs[2].set_yticks([])
    fig.colorbar(cax, ax=axs[2], location='right', fraction=0.1)

    # mu_A = -0.2 eV and T = 0.01 / k
    axs[3] = plot_lattice(final_lattice[0, 3], axs[3], r'$\mu_A = -0.2$ eV, $T = 116K$', 'H', 'N')

    # mu_A = -0.1 eV and T = 0.01 / k
    axs[4] = plot_lattice(final_lattice[3, 3], axs[4], r'$\mu_A = -0.1$ eV, $T = 116K$', 'H', 'N')

    # mu_A = 0 eV and T = 0.01 / k
    axs[5] = plot_lattice(final_lattice[6, 3], axs[5], r'$\mu_A = 0$ eV, $T = 116K$', 'H', 'N')

    plt.tight_layout()
    plt.show()