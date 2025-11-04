#Imports all the functions I can't code as well as a few constants
import numpy as np
from scipy.constants import k, eV
import matplotlib.pyplot as plt

#The variables to run throught all my functions
sizeX = 4
sizeY = 4
k_b = k/eV
#Some filler parameters for testing
Filler_params ={
    'eH': -0.1,
    'eN': -0.1,
    'eNN': 0.05,
    'eHH': 0.05,
    'eNH': -0.05,
    'mu_N': -0.1,
    'mu_H': -0.1,
    'T': 298
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
    
def mister_rogers(sizeX, sizeY):
    neighbors = {}
    for X in range(sizeX):
        for Y in range(sizeY):
            neighbor = [
            (Modulalo((X-1),sizeX), Y),
            (Modulalo((X+1),sizeX), Y),
            (X, Modulalo((Y-1),sizeY)),
            (X, Modulalo((Y+1),sizeY))
            ]
            neighbors[(X,Y)] = neighbor
    return neighbors

def create_grid(sizeX, sizeY):
    raw_grid = np.zeros((sizeX, sizeY))
    defined_grid = mister_rogers(sizeX, sizeY)
    return raw_grid, defined_grid

def energy_calculator(grid, neighbors, parameters):
    Energy = 0
    Nrg = parameters.get('eN')
    Hrg = parameters.get('eH')
    NNrg = parameters.get('eNN')
    HHrg = parameters.get('eHH')
    HNrg = parameters.get('eNH')
    for i, row in enumerate(grid):
        for j, element in enumerate(row):
            if element == 0 or '0':
                pass
            elif element == 'N':
                Adj = neighbors.get((i,j))
                Energy += Nrg
                for position in Adj:
                    X = position[0]
                    Y = position[1]
                    if grid[X][Y] == 0:
                        pass
                    elif grid[X][Y] == 'N':
                        Energy += NNrg
                    elif grid[X][Y] == 'H':
                        Energy += HNrg
            elif element == 'H':
                Adj = neighbors.get((i,j))
                Energy += Hrg
                for position in Adj:
                    X = position[0]
                    Y = position[1]
                    if grid[X][Y] == 0:
                        pass
                    elif grid[X][Y] == 'N':
                        Energy += HNrg
                    elif grid[X][Y] == 'H':
                        Energy += HHrg
    return Energy

    

def gonna_add_one(grid, neighbors, parameters, Empt_Sites, N_Sites, H_Sites):
    N_Pot = parameters.get('mu_N')
    H_Pot = parameters.get('mu_H')
    OldE = energy_calculator(grid, neighbors, parameters)
    beta = 1/(parameters.get('T')*k_b)
    if Empt_Sites == 0:
        return Empt_Sites, N_Sites, H_Sites, grid
    else:
        Valid_Sites = []
        for i, row in enumerate(grid):
            for j, element in enumerate(row):
                if element == 0 or '0':
                    Valid_Sites.append((i, j))
        Chosen_Site = np.random.choice(len(Valid_Sites))
        Chosen_X = Valid_Sites[Chosen_Site][0]
        Chosen_Y = Valid_Sites[Chosen_Site][1]
        N_or_H = np.random.choice(2) > 0.5
        if N_or_H == True:
            new_grid = [[str(ele) for ele in k] for k in grid]
            new_grid[Chosen_X][Chosen_Y] = 'N'
            NewE = energy_calculator(new_grid, neighbors, parameters)
            Acc = min(1, (Empt_Sites - N_Sites) / (N_Sites + 1) * np.exp(-beta * ((NewE - OldE) - N_Pot)))
            if np.random.rand() <= Acc:
                grid = new_grid
                N_Sites += 1
                Empt_Sites -= 1
        elif N_or_H == False:
            new_grid = [[str(ele) for ele in k] for k in grid]
            new_grid[Chosen_X][Chosen_Y] = 'H'
            NewE = energy_calculator(new_grid, neighbors, parameters)
            Acc = min(1, (Empt_Sites - H_Sites) / (H_Sites + 1) * np.exp(-beta * ((NewE - OldE) - H_Pot)))
            if np.random.rand() <= Acc:
                grid = new_grid
                H_Sites += 1
                Empt_Sites -= 1
        return Empt_Sites, N_Sites, H_Sites, grid
        

def take_away_one(grid, neighbors, parameters, Empt_Sites, N_Sites, H_Sites):
    N_Pot = parameters.get('mu_N')
    H_Pot = parameters.get('mu_H')
    OldE = energy_calculator(grid, neighbors, parameters)
    beta = 1/(parameters.get('T')*k_b)
    if N_Sites + H_Sites == 0:
        return Empt_Sites,N_Sites,H_Sites, grid
    else:
        Valid_Sites = []
        for i, row in enumerate(grid):
            for j, element in enumerate(row):
                if element == 'N' or 'H':
                    Valid_Sites.append((i, j))
        Chosen_Site = np.random.choice(len(Valid_Sites))
        Chosen_X = Valid_Sites[Chosen_Site][0]
        Chosen_Y = Valid_Sites[Chosen_Site][1]
        new_grid = [k[:] for k in grid]
        new_grid[Chosen_X][Chosen_Y] = 0
        NewE = energy_calculator(grid, neighbors, parameters)
        if grid[Chosen_X][Chosen_Y] == 'N':
            if Empt_Sites - H_Sites + 1 == 0:
                Denom = 1E-30
            else:
                Denom = Empt_Sites - H_Sites + 1
            Acc = min(1, H_Sites/Denom * np.exp(-beta * ((NewE - OldE) + N_Pot)))
            if np.random.rand() <= Acc:
                grid = new_grid
                N_Sites -= 1
                Empt_Sites += 1
        else:
            if Empt_Sites - H_Sites + 1 == 0:
                Denom = 1E-30
            else:
                Denom = Empt_Sites - H_Sites + 1
            Acc = min(1, H_Sites/Denom * np.exp(-beta * ((NewE - OldE) + H_Pot)))
            if np.random.rand() <= Acc:
                grid = new_grid
                H_Sites -= 1
                Empt_Sites += 1
        return Empt_Sites,N_Sites,H_Sites, grid

def grand_cannonical(Steps, sizeX, sizeY, parameters):
    grid, neighbors = create_grid(sizeX, sizeY)
    Total_Sites = sizeX * sizeY
    Empt_Sites = Total_Sites
    N_Sites = 0
    H_Sites = 0
    NCoverage = []
    HCoverage = []
    for i in range(Steps):
        add_or_subtract = np.random.choice(2) > 0.5
        if add_or_subtract == True:
            Empt_Sites, N_Sites, H_Sites, grid = gonna_add_one(grid, neighbors, parameters, Empt_Sites, N_Sites, H_Sites)
        elif add_or_subtract == False:
            Empt_Sites, N_Sites, H_Sites, grid = take_away_one(grid, neighbors, parameters, Empt_Sites, N_Sites, H_Sites)
        NCoverage.append(N_Sites/Total_Sites)
        HCoverage.append(H_Sites/Total_Sites)
    return grid, NCoverage, HCoverage

grid, grid_data = create_grid(sizeX,sizeY)
#print(grid_data.get((0,1)))

Test_grid = [['N',0,0,0],
             ['H',0,0,0],
             ['H',0,0,0],
             ['N',0,0,0]
]

Grid, Junk1, Junk2 = grand_cannonical(100000,sizeX,sizeY,Filler_params) 
print(Grid)