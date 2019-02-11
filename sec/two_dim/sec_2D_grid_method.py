'''
TO DO
------
Fix Movement/collision
    lines 800 and 801
After fix
    look into using Decimal instead of float for calculationsj:w
    Set up histogram function
    set up sacling grid tile
        scale particles down rather than scaling grid up -keep grids 1x1
    I changed the collision function to not use the Values class, this will
        will need to be rememdied 

Other thoughts

'''

from math import *
import random
import time
import copy

# some imports to explore graphing options
#import matplotlib.plyplot as plt
#import numpy as np

class Values:
    diameter = 0.4  # particle diameter (for hard sphere model)
    ncycles = 504 # number of times a particle will be selceted and moved
    nparticles = 17# number of particles
    # particle density - how close the particles will be at simulation start
    rho = 0.5 
    volume_length = ceil((nparticles / rho) ** 0.5 )# length of the box
    # inverse of the length for use in calculations later
    inv_length = 1 / floor(volume_length)
    length = nparticles ** 0.5 # length the particles will be moved per cycle
    n_root = ceil(sqrt(nparticles / rho)) 
    # set the height of the grid_tils relative to the particles diameter
    grid_tile_height = 1 #(2 * diameter) 
    #particle-particle distances to be sorted and counted, per cycle
    distances = [] 
    #sorted distance, ultimately normalized (summed to one)
    histogram = [] 
    # how many cycles to run before taking a histogram measurement
    cycles_to_average = 2 
    time # the time is takes to run the simulation 
         # doesn't include minor set up and file printing

    #control variables
    pause = False
    extra_print = False#True

class Grid_tile:
    '''Declare a class to make object called grid tiles which have the 
    following attributes: a list of neighboring grids, a list of particles 
    within the grid, and a position '''
    
    def __init__(self, i, j):
        # This method builds a grid tile

        # x and y represent the row and column
        # note, the x and y are used to categorize particles as follows:
        #  a tile contains a  particle if tile. x -1 < particle.x <= tile.x
        #  and tile.y - 1 < particle.y <= particle.y 
        self.x = i #grid tile's x position at top left corner
        self.y = j #grid tile's y position at top left corner

        self.neighbor_list = [[None, None] for index in range(0,8)] 
        # down, up, left, right, up left, up right, down left, down right
        # pre-built empty to optimise memory allocation
        

        # create an empty list which will later be changed as the grid
        # is populated with particles
        self.particle_list = [] 

    def __repr__(self):
        # This method tells python how to represent a tile
        return '[{},{}],'.format(self.x, self.y)

def main():
    '''One function to call them all'''
    #construct an empty square grid
    grid = [[None for column in range(0, Values.n_root)] 
            for row in range(0, Values.n_root)] 

    #pre allocate the memory for the particles
    rxy = [[None,None] for xy in range(0, Values.nparticles)]

    #randomly populate particles in ascending order
    define_list(Values.nparticles, Values.volume_length, Values.diameter,
                rxy)
    start_rxy = copy.deepcopy(rxy)
    #construct a grid that is ceiling(sqrt(n)) x ceiling(sqrt(n)) (int x int)
    build_grid(grid, Values.grid_tile_height) 

    #tell each grid which grids it is next to (non-bound box i.e.wrap around)
    populate_grid_neighbors(grid) 

    random.seed(405) #seed the random number generator for result comparison 
    start = time.time() #record the starting time of the simulation
    
    #assign particles to grids, (populate grid tile particle lists)
    assign_particles_to_grid(rxy, grid)

    #moves random particle in random direction a set distance,
    #repeat for n cycles
    single_event_chain(rxy, grid) 

    end = time.time() # stops timer
    Values.time = (end - start)   #calculates 'total' time of experiment
    
    #constructs a histogram for chemical potential calculations
#    build_histogram(Values.rx)  

    #custom file name
#    filename = "{0}_particles_{1}_cycles_{2}_minutes.txt".format(
#        Values.nparticles, Values.ncycles, Values.time/60)
    # write and output file 
#    write_to_file(filename, Values.rx, Values.rx_init) 
    #build_histogram(Values.rx)
#    print(Values.rx) # print final, sorted list to terminal    
    # custom file name
#    filename = "{0}_particles_{1}_cycles_{2}_minutes.txt".format(
#        Values.nparticles, Values.ncycles, Values.time/60) 
    # make a new file using custom file name
#    write_to_file(filename, Values.rx, Values.rx_init) 


    print('\n\n')
    print_particles_in_grid(grid)
    rxy = sorted(rxy)
    print('\n\n')
    print("Start\t\t\tEnd")
    for i in range(0,len(start_rxy)):
        pass
        print("{} {} \t {}".format(i, start_rxy[i], rxy[i]))

    temp_list = [] 
    for i in range(0, len(rxy) - 1):
        for j in range(i + 1, len(rxy)):
            x_dist = abs(rxy[i][0] - rxy[j][0])
            y_dist = abs(rxy[i][1] - rxy[j][1])
            temp_list.append([sqrt(x_dist**2 + y_dist**2),i,j])
    print('\n\n','\nVector Distances < diameter: \n------------------\n')
    temp_list = sorted(temp_list)
    for i in temp_list:
        if (i[0] < Values.diameter):
            print(i)

#    print(rxy) # print final, sorted list to terminal for quick review without opening file

def build_grid(input_grid, grid_dim = 1): 
    """
    O(n^2)
    Creates a grid to hold grid_tiles
    Each tile has height and width equal to Values.grid_tile_height
    """
    # need to update this for scaling particles instead of tiles

    for i in range(0, len(input_grid[0])):
       for j in range(0, len(input_grid[0])):
           input_grid[i][j] = Grid_tile((i * grid_dim), (j * grid_dim))
        #endfor
    #endfor

def populate_grid_neighbors(input_grid): 
    '''
    O(n)
    Determines which grid_tiles are next to each other
    '''
    grid_size = len(input_grid[0]) 
    for i in range(0, grid_size ** 2):
        # indices start at 0 and end at n-1
        # integer division and remainder to interate though rows and columns: 

        if (i == 0): # first grid_tile, bottom left
            input_grid[i][i].neighbor_list = [
                    [0, grid_size - 1], # down
                    [0, grid_size - (grid_size - 1)], # up
                    [grid_size - 1, 0], # left
                    [grid_size - (grid_size - 1), 0], #right
                    [grid_size - 1, grid_size - (grid_size - 1)], #upL
                    [grid_size - (grid_size - 1), grid_size - (grid_size - 1)], #upR
                    [grid_size - 1, grid_size - 1], #downL
                    [grid_size - 1, grid_size - (grid_size - 1)]] #downR
        #endif
        elif (i == ((grid_size ** 2) - 1)): # last grid_tile, top right
            input_grid[grid_size - 1][grid_size - 1].neighbor_list = [
                    [grid_size - 1, grid_size - 2], # down
                    [grid_size - 1, 0], #up
                    [grid_size - 2, grid_size - 1], # left
                    [0, grid_size - 1], #right
                    [grid_size - 2, 0], #upL
                    [0, 0], # upR
                    [grid_size - 2, grid_size - 2], # downL
                    [0, grid_size - 2]] # downR
        #endelse
        elif (i == (grid_size - 1)): # bottom right
            input_grid[i][0].neighbor_list = [
                    [grid_size - 1, grid_size - 1], #down
                    [grid_size - 1, grid_size - (grid_size - 1)], #up
                    [grid_size - 2, 0], #l
                    [0, 0], #R
                    [grid_size - 2, grid_size - (grid_size - 1)], #uL
                    [0, grid_size - (grid_size - 1)],  #uR
                    [grid_size - 2, grid_size - 1], #dl
                    [0, grid_size - 1]] #dR
        #endelse
        elif (i == ((grid_size ** 2) - grid_size)):# top left
            input_grid[0][grid_size - 1].neighbor_list = [
                    [0, grid_size - 2], #down
                    [0, 0], #up
                    [grid_size - 1, grid_size - 1], #L
                    [grid_size - (grid_size - 1), grid_size - 1], #R
                    [grid_size - 1, 0], #uL
                    [grid_size - (grid_size - 1), 0], #uR
                    [grid_size - 1, grid_size - 2], #dL
                    [grid_size - (grid_size - 1), grid_size - 2]] #dR
        #endelse            
        elif ((i // grid_size) == 0 ): # bottom row
            input_grid[i % grid_size][i // grid_size].neighbor_list = [
                    [i % grid_size, grid_size - 1], #down
                    [ i % grid_size, grid_size - (grid_size - 1)], #up
                    [(i - 1) % grid_size, 0], #L
                    [(i + 1) % grid_size, 0], #R
                    [(i - 1) % grid_size, grid_size - (grid_size - 1)], #uL
                    [(i + 1) % grid_size, grid_size - (grid_size - 1)], #uR
                    [(i - 1) % grid_size, grid_size - 1], #dL
                    [(i + 1) % grid_size, grid_size - 1]] #dR
        #endelse
#        elif ((grid_size ** 2 - grid_size) < i < ((grid_size ** 2) - 2)): # top row
        elif ((i // grid_size) == (grid_size - 1)):
            input_grid[i % grid_size][i // grid_size].neighbor_list = [
                    [i % grid_size, grid_size - 2], #up
                    [i % grid_size, 0], #down
                    [(i - 1) % grid_size, grid_size - 1], #L
                    [(i + 1) % grid_size, grid_size - 1], #R
                    [(i - 1) % grid_size, 0], #uL
                    [(i + 1) % grid_size, 0], #uR
                    [(i - 1) % grid_size, grid_size - 2], #dL
                    [(i + 1) % grid_size, grid_size - 2]] #dR
        #endelse
        elif (i % grid_size == 0): #left column
            input_grid[i % grid_size][i // grid_size].neighbor_list = [
                    [i % grid_size, (i // grid_size) - 1], #down
                    [i % grid_size, (i // grid_size) + 1], #up
                    [grid_size - 1, i // grid_size], #L
                    [grid_size - (grid_size - 1), i // grid_size], #R
                    [grid_size - 1, (i // grid_size) + 1], #uL
                    [grid_size - (grid_size - 1), (i // grid_size) + 1], #uR
                    [grid_size - 1, (i // grid_size) - 1], #dL
                    [grid_size - (grid_size - 1), (i // grid_size) - 1]] #dR
        #endelse
        elif (i % grid_size == (grid_size - 1)): #right column
            input_grid[i % grid_size][i // grid_size].neighbor_list = [
                    [i % grid_size, (i // grid_size) - 1], #down
                    [i % grid_size, (i // grid_size) + 1], #up
                    [grid_size - 2, i // grid_size], #L
                    [0, i // grid_size], #R
                    [grid_size - 2, (i // grid_size) + 1], #uL
                    [0, (i // grid_size) + 1], #uR
                    [grid_size - 2, (i // grid_size) - 1], #dL
                    [0, (i // grid_size) - 1]] #dR
        #endelse
        else: # everything else
            input_grid[i % grid_size][i // grid_size].neighbor_list = [
                    [(i % grid_size), (i // grid_size) - 1], #down
                    [i % grid_size, (i // grid_size) + 1], #up
                    [(i - 1) % grid_size, i // grid_size], #L
                    [(i + 1) % grid_size, i // grid_size], #R
                    [(i - 1) % grid_size, (i // grid_size) + 1], #uL
                    [(i + 1) % grid_size, (i // grid_size) + 1], #uR
                    [(i - 1) % grid_size, (i // grid_size) - 1], #dL
                    [(i + 1) % grid_size, (i // grid_size) - 1]] #dR
        #endless
    #endfor

def define_list(nparticles, volume_length, diameter, rxy):
    '''
    automatically, the particles are sorted in ascending order
    one particle per grid space until no more particles

    divide the particles by grid length and start filling the
    grid (it doesn't matter how they start)
    '''
    
    for i in range(0, nparticles):
        rxy[i] = [float(ceil(i % (volume_length - 1))), 
                    float(i // (volume_length - 1))]
    #enddo    

def assign_particles_to_grid(rxy, grid): #O(n) or O(n ** dimension) 
    '''Assigns particles to the grid tile. A particle belongs to the grid for
    which the particle x and y are greater than or equal to the tile x and y 
    but less than the tile (x+1) and (y+1)'''

    # This will need to be updated to reflect scaling the spheres instead
    for i in range(0, len(rxy)):
        x = floor(rxy[i][0] / Values.grid_tile_height)
        y = floor(rxy[i][1] / Values.grid_tile_height)
        grid[x][y].particle_list.append(rxy[i])

def print_particles_in_grid(grid, grid_labels = True):
    """Prints particles in each grid tile in a grid orientation. Origin is in
    bottom left. Prints grid labels by deault and prints and empty grid as 
    'empty'. Output is not (yet) evenly spaced."""

    for j in range(len(grid[0])-1, -1, -1):
        for i in range(0, len(grid[0])):
            if (grid_labels):
                print('|    ',grid[i][j],end ="    | ")
        print()
        for i in range(0, len(grid[0])):
            if (len(grid[i][j].particle_list) == 0):
                print('|    Empty', end ='    | ')
            else:
                temp = [[round(value, 2) for value in pair] for pair in grid[i][j].particle_list]
                print('| ', temp, end =' | ')
        print()
        print("-----------------------------------------------------------------------------------")

def single_event_chain(rxy, grid):
    for i in range(0, Values.ncycles):
        particle, tile, position = select_particle_and_tile(grid)
        direction = random.choice(['UP','DOWN','LEFT','RIGHT'])
        direction = 'UP'

        move(direction, rxy, grid, particle, tile, position)
        if (Values.extra_print):
            if (i%10 == 0):
                print("cycle: ",i)
        
def select_particle_and_tile(grid, override = False, override_num = 0):
    if (override):
        # override for debugging and testing purposes to avoid RNG
        return grid[0][override_num].particle_list[0], grid[0][override_num], 0
    while (True):
        #choose a random grid tile from 2D grid, grid[row][column]
        the_chosen_one = random.choice(random.choice(grid))

        #check that there are particles in the grid
        if (len(the_chosen_one.particle_list) != 0):
        #choose a random particle from chosen grid tile
            random_number = random.randint(0, 
                                    len(the_chosen_one.particle_list) - 1)
            rand_particle = the_chosen_one.particle_list[random_number]
            break
    # returns a particle, a tile, and a position in the list
    return rand_particle,the_chosen_one, random_number

def check_for_collisions(direction, grid, rand_particle, tile, volume_length, inv_length):
    """ grid is the list of tiles, each tile has a list tof particles
    rand_particle is the moving particle
    tile is the tile containing the moving particle
    volume_length, and inv_length were included for debugging reasons (to allow
    different grid sizes more easily
    """
    if (Values.extra_print):
        print("tile: ", tile, "neighbor list: ", tile.neighbor_list)
    x_distance_btwn = 0
    y_distance_btwn = 0
    overlap_count = 0
    next_particle = [[-9.0,[-9.0,-9.0]]] #dummy data to please the compiler
    neighbor_tiles = None
    neighbor_tile = None

    moving = rand_particle

    if (direction == 'UP'):
        #"first check self"
        tile.particle_list
        #"moving up: want to check tile's neighbors up, up-l, up-r, l, and R"
        # up, left, right, up-left, up-right
        neighbor_tiles = [1, 2, 3, 4, 5]
        '''Check collison up'''
        # up, left, right, up-left, up-right
        neighbor_tiles = [1, 2, 3, 4, 5]

        # check self for collisions (in the case of multiple particles fitting
        # in a single grid

        # checking same grid as moving particle
        if (Values.extra_print):
            print("checking within grid")
        for stationary in tile.particle_list:
            if (Values.extra_print):
                print('moving: ', moving, 'stationary: ', stationary)
            x_result = abs(moving[0] - stationary[0])
            if (Values.extra_print):
                print('x result: ', x_result)
            if (x_result <= Values.diameter):
                y_result = moving[1] - stationary[1]
            if (Values.extra_print):
                print('y result: ', y_result)
                if (y_result < 0):
                    distance_to_collision = ( stationary[1] - moving[1] -
                                            sqrt(4*(Values.diameter*0.5)**2 -
                                            abs(stationary[0] - moving[0])**2) )
                    overlap_count += 1
                    next_particle.append([distance_to_collision, stationary,
                                        [tile.x, tile.y], distance_to_collision])


        if (Values.extra_print):
            print("checking above grid")
        if (overlap_count == 0):
            #check neighbors for collisions
            for i in neighbor_tiles:
                # ith neighbor in the neighbor list
                neighbor_tile = tile.neighbor_list[i]
                neighbor = grid[neighbor_tile[0]][neighbor_tile[1]]
                for stationary in neighbor.particle_list:
                    if (Values.extra_print):
                        print('moving: ', moving, 'stationary: ', stationary)
                    x_result = abs(moving[0] - stationary[0])
                    if (Values.extra_print):
                        print('x result: ', x_result)
                    if (x_result > 4*Values.diameter):
                        x_result = abs(x_result - volume_length)
                    if (x_result <= Values.diameter):
                        y_result = moving[1] - stationary[1]
                        if (Values.extra_print):
                            print('y result: ', y_result)
                        if (y_result < 0):
                            distance_to_collision = ( stationary[1] - moving[1] -
                                                    sqrt(4*(Values.diameter*0.5)**2 -
                                                   (stationary[0] - moving[0])**2) )
                            overlap_count += 1
                            next_particle.append([distance_to_collision, stationary,
                                                tile.neighbor_list[i], distance_to_collision])

                        if (y_result > 4*Values.diameter):
                            y_result = moving[1] - (stationary[1] + volume_length) 
                            if (Values.extra_print):
                                print('y resul pbc: ', y_result)
                            distance_to_collision = ( (stationary[1] + volume_length) - moving[1] -
                                                    sqrt(4*(Values.diameter*0.5)**2 -
                                                   (stationary[0] - moving[0])**2) )
                            overlap_count += 1
                            next_particle.append([distance_to_collision, stationary,
                                                tile.neighbor_list[i], distance_to_collision])

        if (overlap_count == 0):
            next_particle.append([Values.diameter, moving,
                                [tile.x, tile.y], Values.diameter])
            



    if (direction == 'UP2'):
        '''Check collison up'''
        # up, left, right, up-left, up-right
        neighbor_tiles = [1, 2, 3, 4, 5]

        # check self for collisions (in the case of multiple particles fitting
        # in a single grid
        for particle in tile.particle_list:
            # don't collide with self
            if (particle == rand_particle):
                continue

            #pbc stands for periodic boundary conditions
            y_distance_btwn = abs(rand_particle[1] - particle[1])
            # the rounding forces either a 0 or 1
            y_distance_pbc = abs(y_distance_btwn - floor(volume_length + 1)*
                            round(y_distance_btwn*inv_length))
            if (Values.extra_print):
                print('ydist: ', y_distance_btwn, 'y pbc: ', y_distance_pbc)

            # no collision: the particle are farther apart than the distance
            # being moved
            if (y_distance_pbc > Values.grid_tile_height):
                continue

            # no collision: the particles are within striking distance but
            # the moving particle is moving away from the stationary particle
            if ((y_distance_pbc <= Values.diameter) and
                 (rand_particle[1] > particle[1])):
                continue

            # particles are within y striking distance: check x distances
            x_distance_btwn = abs(rand_particle[0] - particle[0])
            x_distance_pbc = abs(x_distance_btwn - floor(volume_length + 1)*
                                round(x_distance_btwn*inv_length))
            if (Values.extra_print):
                print('xdist: ', x_distance_btwn, 'x pbc: ', x_distance_pbc)

            # no collision: radius + radius = diameter
            # the x distance between is greater than the diameter
            if (x_distance_pbc >= Values.diameter):
                continue

            # collision!
            overlap_count += 1
            dy = sqrt(abs(Values.diameter**2 - x_distance_pbc**2))
            next_particle.append([sqrt(x_distance_pbc**2 +
                                y_distance_pbc**2), particle,
                                [tile.x//Values.grid_tile_height,
                                tile.y//Values.grid_tile_height],
                                abs(abs(y_distance_pbc) - dy)])

        if (overlap_count == 0):
            #check neighbors for collisions
            for i in neighbor_tiles:
                # ith neighbor in the neighbor list
                neighbor_tile = tile.neighbor_list[i]
                neighbor = grid[neighbor_tile[0]][neighbor_tile[1]]
                for particle in neighbor.particle_list:
                    #don't collide with self (one particle in two tiles)
                    if (particle == rand_particle):
                        continue
                    
                    y_distance_btwn = abs(rand_particle[1] - particle[1])
                    y_distance_pbc = abs(y_distance_btwn - floor(volume_length + 1)*
                                round(y_distance_btwn*inv_length))

                    # no collision
                    if (y_distance_pbc > Values.diameter):
                        continue
                    if ((y_distance_pbc <= Values.diameter) and
                         (rand_particle[1] > particle[1])):# and (y_distance_pbc >  0):
                        continue

                    # if a y collision is detected proceed to check for x
                    x_distance_btwn = abs(rand_particle[0] - particle[0])
                    x_distance_pbc = abs(x_distance_btwn - floor(volume_length + 1)*
                                    round(x_distance_btwn*inv_length))
                    if (x_distance_pbc >= Values.diameter):
                        continue

 #                   print('d: ', Values.diameter, 'x: ', x_distance_pbc)
                    # if collisions exist, record them
                    overlap_count += 1
                    dy = sqrt(abs(Values.diameter**2 - x_distance_pbc**2))
                    next_particle.append([sqrt(x_distance_pbc**2 +
                                        y_distance_pbc**2), particle,
                                        tile.neighbor_list[i],
                                        abs(abs(y_distance_pbc) - dy)])
    elif (direction == 'DOWN'):
        '''Check collison down'''
        # down, left, right, down-left, down-right
        neighbor_tiles = [0, 2, 3, 6, 7]

        #check self for collisions
        for particle in tile.particle_list:
            # don't collide with self
            if (particle == rand_particle):
                continue
            # don't check particles above this particle (no PBC in same tile)
            if (particle[1] > rand_particle[1]):
                continue
            else:
                y_distance_btwn = abs(rand_particle[1] - particle[1])
                y_distance_pbc = (y_distance_btwn - Values.volume_length*
                                round(y_distance_btwn*Values.inv_length))
                

            x_distance_btwn = abs(rand_particle[0] - particle[0])
            x_distance_pbc = (x_distance_btwn - Values.volume_length*
                        round(x_distance_btwn*Values.inv_length))

            # if collisions exist, record them
            if (x_distance_pbc < Values.diameter):
                if (y_distance_pbc < Values.grid_tile_height):
                    overlap_count += 1
                    dy = sqrt(Values.diameter**2 - x_distance_pbc**2)
                    next_particle.append([sqrt(x_distance_pbc**2 +
                                        y_distance_pbc**2), particle,
                                        [tile.y//Values.grid_tile_height,
                                        tile.x//Values.grid_tile_height],
                                        y_distance_pbc - dy])

        #check neighbors for collisions
        for i in neighbor_tiles:
            # ith neighbor in the neighbor list
            neighbor_tile = tile.neighbor_list[i]
            neighbor = grid[neighbor_tile[0]][neighbor_tile[1]]
            for particle in neighbor.particle_list:
                #don't collide with self (one particle in two tiles)
                if (particle == rand_particle):
                    continue
                x_distance_btwn = abs(rand_particle[0] - particle[0])
                # dont check particles above (consider PBC)
                if (particle[1] < rand_particle[1]):
                    y_distance_btwn = abs(rand_particle[1] - particle[1])
                    y_distance_pbc = (y_distance_btwn - Values.volume_length*
                                round(y_distance_btwn*Values.inv_length))

                elif ((particle[1] + Values.volume_length) > 
                        (rand_particle[1] )):
                    y_distance_btwn = abs(rand_particle[1] -
                                        (particle[1] + Values.volume_length))
                    y_distance_pbc = (y_distance_btwn - Values.volume_length*
                                round(y_distance_btwn*Values.inv_length))
                else:
                    continue

                # if collisions exist, record them
                if (x_distance_btwn < Values.diameter):
                    if (y_distance_pbc < Values.grid_tile_height):
                        overlap_count += 1
                        dy = sqrt(Values.diameter**2 - x_distance_btwn**2)
                        next_particle.append([sqrt(x_distance_btwn**2 +
                                            y_distance_pbc**2), particle,
                                            tile.neighbor_list[i],
                                            y_distance_pbc + dy])

    elif (direction == 'LEFT'):
        '''Check collisions left'''
        # down, up, left, up-left, down-left
        neighbor_tiles = [0, 1, 2, 4, 6]
        #check self for collisions
        for particle in tile.particle_list:
            # don't collide with self
            if (particle == rand_particle):
                continue
            # don't check particles to left (no PBC in same tile)
            if (particle[0] < rand_particle[0]):
                continue
            else:
                x_distance_btwn = abs(rand_particle[0] - particle[0])
                x_distance_pbc = (x_distance_btwn - Values.volume_length*
                                round(x_distance_btwn*Values.inv_length))

            y_distance_btwn = abs(rand_particle[1] - particle[1])

            # if collisions exist, record them
            if (y_distance_btwn < Values.diameter):
                if (x_distance_pbc < Values.grid_tile_height):
                    overlap_count += 1
                    dx = sqrt(Values.diameter**2 - y_distance_btwn**2)
                    next_particle.append([sqrt(x_distance_pbc**2 +
                                        y_distance_btwn**2), particle,
                                        [tile.y//Values.grid_tile_height,
                                        tile.x//Values.grid_tile_height],
                                        x_distance_pbc + dx])

        #check neighbors for collisions
        for i in neighbor_tiles:
            # ith neighbor in the neighbor list
            neighbor_tile = tile.neighbor_list[i]
            neighbor = grid[neighbor_tile[0]][neighbor_tile[1]]
            for particle in neighbor.particle_list:
                #don't collide with self (one particle in two tiles)
                if (particle == rand_particle):
                    continue
                y_distance_btwn = abs(rand_particle[1] - particle[1])
                # dont check particles to right (consider PBC)
                if (particle[0] > rand_particle[0]):
                    x_distance_btwn = abs(rand_particle[0] - particle[0])
                    x_distance_pbc = (x_distance_btwn - Values.volume_length*
                                round(x_distance_btwn*Values.inv_length))

                elif ((particle[0] + Values.volume_length) > 
                        (rand_particle[0] )):
                    x_distance_btwn = abs(rand_particle[0] -
                                        (particle[0] + Values.volume_length))
                    x_distance_pbc = (x_distance_btwn - Values.volume_length*
                                round(x_distance_btwn*Values.inv_length))
                else:
                    continue

                # if collisions exist, record them
                if (y_distance_btwn < Values.diameter):
                    if (x_distance_btwn < Values.grid_tile_height):
                        overlap_count += 1
                        dx = sqrt(Values.diameter**2 - y_distance_btwn**2)
                        next_particle.append([sqrt(x_distance_pbc**2 +
                                            y_distance_btwn**2), particle,
                                            tile.neighbor_list[i],
                                            x_distance_pbc - dx])

    elif (direction == 'RIGHT'):
        '''Check collisions right'''
        # down, up, right, up-right, down-right
        neighbor_tiles = [0, 1, 3, 5, 7]
        #check self for collisions
        for particle in tile.particle_list:
            # don't collide with self
            if (particle == rand_particle):
                continue
            # don't check particles to the right (no PBC in same tile)
            if (particle[0] > rand_particle[0]):
                continue
            else:
                x_distance_btwn = abs(rand_particle[0] - particle[0])
                x_distance_pbc = (x_distance_btwn - Values.volume_length*
                                round(x_distance_btwn*Values.inv_length))

            y_distance_btwn = abs(rand_particle[1] - particle[1])

            # if collisions exist, record them
            if (y_distance_btwn < Values.diameter):
                if (x_distance_btwn < Values.grid_tile_height):
                    overlap_count += 1
                    dx = sqrt(Values.diameter**2 - y_distance_btwn**2)
                    next_particle.append([sqrt(x_distance_pbc**2 +
                                        y_distance_btwn**2), particle,
                                        [tile.y//Values.grid_tile_height,
                                        tile.x//Values.grid_tile_height],
                                        x_distance_pbc - dx])

        #check neighbors for collisions
        for i in neighbor_tiles:
            # ith neighbor in the neighbor list
            neighbor_tile = tile.neighbor_list[i]
            neighbor = grid[neighbor_tile[0]][neighbor_tile[1]]
            for particle in neighbor.particle_list:
                #don't collide with self (one particle in two tiles)
                if (particle == rand_particle):
                    continue
                y_distance_btwn = abs(rand_particle[1] - particle[1])
                # dont check particles to left (consider PBC)
                if (particle[0] > rand_particle[0]):
                    x_distance_btwn = abs(rand_particle[0] - particle[0])
                    x_distance_pbc = (x_distance_btwn - Values.volume_length*
                                round(x_distance_btwn*Values.inv_length))

                elif ((particle[0] + Values.volume_length) < 
                        (rand_particle[0] )):
                    x_distance_btwn = abs(rand_particle[0] -
                                        (particle[0] + Values.volume_length))
                    x_distance_pbc = (x_distance_btwn - Values.volume_length*
                                round(x_distance_btwn*Values.inv_length))
                else:
                    continue

                # if collisions exist, record them
                if (y_distance_btwn < Values.diameter):
                    if (x_distance_btwn < Values.grid_tile_height):
                        overlap_count += 1
                        dx = sqrt(Values.diameter**2 - y_distance_btwn**2)
                        next_particle.append([sqrt(x_distance_pbc**2 +
                                            y_distance_btwn**2), particle,
                                            tile.neighbor_list[i],
                                            x_distance_btwn - dx])


    # returns the number of overlaps, the list of overlapping particles,
    # and the tile in which the particles were found
    return overlap_count, sorted(next_particle)

def move(direction, rxy, grid, particle, tile, position, volume_length = Values.volume_length):
    '''Legacy Comment'''
    #Rather than adding up, I decided to subtract - am I pessimistic?

    if (Values.extra_print):
        print("pre move")
        print_particles_in_grid(grid, True)
    distance = 1*Values.length
    x_distance_btwn = None
    y_distance_btwn = None
    neighbor_tile = None
    start_rxy = copy.deepcopy(rxy)
    if (Values.extra_print):
        print('\n')
        print("list{}: {} moving: {}".format([tile.x, tile.y], tile.particle_list, tile.particle_list[position]))

    if direction == 'UP':

        '''check particle for collisions with grid above and left/right
        move particle from grid list below to grid list above
        remove particle from original list
        update particle position and movement'''
        while (distance > 0):
#        for counter in range(11):
            if (Values.extra_print):
                print("distance: ",distance)
                print('\n\n')
                print("Start\t\t\tEnd")
                for i in range(0,len(start_rxy)):
                    print("{} {} \t {}".format(i, start_rxy[i], rxy[i]))

            # do some calculations for printing debugging information?
            # no, this is to build the histogram
            if (Values.extra_print):
                temp_list = [] 
                for i in range(0, len(rxy) - 1):
                    for j in range(i + 1, len(rxy)):
                        x_dist = abs(rxy[i][0] - rxy[j][0])
                        y_dist = abs(rxy[i][1] - rxy[j][1])
                        temp_list.append([sqrt(x_dist**2 + y_dist**2),i,j])

            # check for collisions
            if (Values.extra_print):
                print('in overlap check')
            overlap_count, next_particle = check_for_collisions(
                direction, grid, particle, tile, volume_length,
                Values.inv_length)

            # no collisions
            if (overlap_count == 0):

                if (Values.extra_print):
                    print("No particle collision")
                # remove particle from starting grid
                #temp_tile = next_particle[1][2]
                temp = next_particle[1][1]
                if (Values.extra_print):
                    print("particle moving: ", temp)
                
                # start modifying the y-value of that particle
                if (Values.extra_print):
                    print("start temp[1]: ", temp[1])
                # update particle (moving up)
                if (distance > Values.diameter):
                    temp[1] += Values.diameter
        #            if (temp[1] > volume_length):
        #                temp[1] -= volume_length
                    if (Values.extra_print):
                        print("updated temp[1]: ", temp[1])
                else:
                    temp[1] += distance
        #            if (temp[1] > volume_length):
        #                temp[1] -= volume_length
                    if (Values.extra_print):
                        print("update 2 temp[1]: ", temp[1])
                # periodic boundary condition
                # determine destination grid (up is second in the list)
                y = floor(temp[1] / Values.grid_tile_height)
                if (y - tile.y > 0):
                    destination = (grid[tile.neighbor_list[1][0]]
                                     [tile.neighbor_list[1][1]])
                else:
                    destination = tile

                if (temp[1] >= ((volume_length) )):
                    temp[1] = temp[1] - floor(volume_length)
                    if (Values.extra_print):
                        print("pbc temp[1]: ", temp[1])
                # put particle in destination grid
                destination.particle_list.append(temp)
                # update cycle movement
                distance -=  Values.diameter  # pfc chan
                # move the next particle (last particle in the new list)
                particle = destination.particle_list[-1]
                position = len(destination.particle_list) - 1
                tile = destination

                if (Values.extra_print):
                    print("\n\npost move-no collision")
                    print_particles_in_grid(grid, grid_labels = True)

            # collision
            # next_particle should = [distance to nearest particle, particle,
            #                           grid that particle is in]
            if (overlap_count > 0):
                # the closest collision
                # second value because the first is broken, dummy data
                distance_to_particle, particle_hit, neighbor_tile, delta_y = next_particle[1]
                if (Values.extra_print):
                    print("particle_hit: {} particle moving: {}".format(
                        particle_hit, particle))
                #the index of the particle that will be hit
                neigh_x = neighbor_tile[0]
                neigh_y = neighbor_tile[1]
                if (Values.extra_print):
                    print('neighbor tile: ', neighbor_tile)
                    print('neigh x: ', neigh_x, 'neigh y: ', neigh_y)
                    print(grid[neigh_x][neigh_y].particle_list)
                index = grid[neigh_x][neigh_y].particle_list.index(particle_hit)

                # update particle (moving up)
                if (distance > delta_y):
                    particle[1] += delta_y
                else:
                    particle[1] += distance
                # determine destination grid (up is second in the list)
                y = floor(particle[1] / Values.grid_tile_height)
                if (Values.extra_print):
                    print('y: ', y, "tile.y: ", tile.y)
                if (y - tile.y > 0):
                    destination = (grid[tile.neighbor_list[1][0]]
                                     [tile.neighbor_list[1][1]])
                else:
                    destination = tile
                # periodic boundary condition
                if (particle[1] >= (floor(volume_length) )):
                    particle[1] = particle[1] - floor(volume_length )
                # put particle in destination grid
                tile.particle_list.pop()
                destination.particle_list.append(particle)
                # update cycle movement
                if (Values.extra_print):
                    print("distance: {} delta_y: {}".format(distance, delta_y))
                distance -= delta_y
                # move the next particle (particle being hit)
                particle = particle_hit
                position = index
                tile = grid[neigh_x][neigh_y]

                if (Values.extra_print):
                    print("\n\npost move-collision")
                    print_particles_in_grid(grid, grid_labels = True)



                if (Values.pause):
                    input()

    if direction == 'DOWN':
        '''check particle for collisions with grid below and left/right
        move particle from grid list below to grid list above
        remove particle from original list
        update particle position and movement'''
        while (distance > 0):
            if (Values.extra_print):
                print("distance: ",distance)

            # check for collisions
            overlap_count, next_particle = check_for_collisions(
                direction, grid, particle, tile)

            # no collisions
            if (overlap_count == 0):
                if (Values.extra_print):
                    print("particles: {} position: {}".format(
                            tile.particle_list, position))

                # remove particle from starting grid
                temp = tile.particle_list.pop(position)
                # update particle (moving up)
                # perhaps this should be plus: look into how the grid is implemented
                temp[1] -= Values.grid_tile_height
                # periodic boundary condition
                if (temp[1] < 0 ):
                    temp[1] = temp[1] + Values.volume_length
                # determine destination grid (down is first in the list)
                destination = (grid[tile.neighbor_list[0][0]]
                                 [tile.neighbor_list[0][1]])
                # put particle in destination grid
                destination.particle_list.append(temp)
                # update cycle movement
                distance -=  Values.grid_tile_height
                # move the next particle (last particle in the new list)
                particle = destination.particle_list[-1]
                position = len(destination.particle_list) - 1
                tile = destination

            # collision
            # next_partcile should = [distance to nearest particle, particle,
            #                           grid that particle is in]

            if (overlap_count > 0):
                # the closest collision
                # second value because the first is broken, dummy data
                distance_to_particle, particle_hit, neighbor_tile, delta_y = next_particle[1]
                #the index of the particle that will be hit
                neigh_x = neighbor_tile[0]
                neigh_y = neighbor_tile[1]
                index = grid[neigh_x][neigh_y].particle_list.index(
                    particle_hit)

                # update particle (moving up)
                # perhaps this should be minus: look into how the grid is implemented
                particle[1] += delta_y
                # periodic boundary condition
                if (particle[1] < 0):
                    particle[1] = particle[1] + Values.volume_length

                # update cycle movement
                distance -= delta_y
                # move the next particle (particle being hit)
                particle = particle_hit
                position = index
                tile = grid[neigh_x][neigh_y]

                if (Values.pause):
                    input()

    if direction == 'LEFT':
        '''check particle for collisions with grid left and up/down
        move particle from grid list below to grid list above
        remove particle from original list
        update particle position and movement'''
        while (distance > 0):
            if (Values.extra_print):
                print("distance: ",distance)

            # check for collisions
            overlap_count, next_particle = check_for_collisions(
                direction, grid, particle, tile)

            # no collisions
            if (overlap_count == 0):
                if (Values.extra_print):
                    print("particles: {} position: {}".format(
                            tile.particle_list, position))

                # remove particle from starting grid
                temp = tile.particle_list.pop(position)
                # update particle (moving up)
                # perhaps this should be minus: look into how the grid is implemented
                temp[0] -= Values.grid_tile_height
                # periodic boundary condition
                if (temp[0] < 0):
                    temp[0] = temp[0] + Values.volume_length
                # determine destination grid (left is second in the list)
                destination = (grid[tile.neighbor_list[3][0]]
                                 [tile.neighbor_list[3][1]])
                # put particle in destination grid
                destination.particle_list.append(temp)
                # update cycle movement
                distance -=  Values.grid_tile_height
                # move the next particle (last particle in the new list)
                particle = destination.particle_list[-1]
                position = len(destination.particle_list) - 1
                tile = destination

            # collision
            # next_partcile should = [distance to nearest particle, particle,
            #                           grid that particle is in]

            if (overlap_count > 0):
                # the closest collision
                # second value because the first is broken, dummy data
                distance_to_particle, particle_hit, neighbor_tile, delta_y = next_particle[1]
                if (Values.extra_print):
                    print("particle_hit: {} particle moving: {}".format(
                        particle_hit, particle))
                #the index of the particle that will be hit
                neigh_x = neighbor_tile[0]
                neigh_y = neighbor_tile[1]
                if (Values.extra_print):
                    print("particle list{}: {}".format([neigh_x, neigh_y], grid[neigh_x][neigh_y].particle_list))
                index = grid[neigh_x][neigh_y].particle_list.index(
                    particle_hit)

                # update particle (moving left)
                particle[0] -= delta_y
                # periodic boundary condition
                if (particle[0] < 0):
                    particle[0] = particle[0] + Values.volume_length
                # update cycle movement
                if (Values.extra_print):
                    print("distance: {} delta_y: {}".format(distance, delta_y))
                distance -= delta_y
                # move the next particle (particle being hit)
                particle = particle_hit
                position = index
                tile = grid[neigh_x][neigh_y]

                if (Values.pause):
                    input()

    if direction == 'RIGHT':
        '''check particle for collisions with grid right and up/down
        move particle from grid list below to grid list above
        remove particle from original list
        update particle position and movement'''
        while (distance > 0):
            if (Values.extra_print):
                print("distance: ",distance)

            # check for collisions
            overlap_count, next_particle = check_for_collisions(
                direction, grid, particle, tile)

            # no collisions
            if (overlap_count == 0):
                if (Values.extra_print):
                    print("particles: {} position: {}".format(
                            tile.particle_list, position))

                # remove particle from starting grid
                temp = tile.particle_list.pop(position)
                # update particle (moving right)
                temp[0] += Values.grid_tile_height
                # periodic boundary condition
                if (temp[0] > Values.volume_length):
                    temp[0] = temp[0] - Values.volume_length
                # determine destination grid (left is second in the list)
                destination = (grid[tile.neighbor_list[3][0]]
                                 [tile.neighbor_list[3][1]])
                # put particle in destination grid
                destination.particle_list.append(temp)
                # update cycle movement
                distance -=  Values.grid_tile_height
                # move the next particle (last particle in the new list)
                particle = destination.particle_list[-1]
                position = len(destination.particle_list) - 1
                tile = destination

            # collision
            # next_partcile should = [distance to nearest particle, particle,
            #                           grid that particle is in]

            if (overlap_count > 0):
                # the closest collision
                # second value because the first is broken, dummy data
                distance_to_particle, particle_hit, neighbor_tile, delta_y = next_particle[1]
                if (Values.extra_print):
                    print("particle_hit: {} particle moving: {}".format(
                        particle_hit, particle))
                #the index of the particle that will be hit
                neigh_x = neighbor_tile[0]
                neigh_y = neighbor_tile[1]
                if (Values.extra_print):
                    print("particle list{}: {}".format([neigh_x, neigh_y], grid[neigh_x][neigh_y].particle_list))
                index = grid[neigh_x][neigh_y].particle_list.index(
                    particle_hit)

                # update particle (moving left)
                particle[0] += delta_y
                # periodic boundary condition
                if (particle[0] > Values.volume_length):
                    particle[0] = particle[0] - Values.volume_length
                # update cycle movement
                if (Values.extra_print):
                    print("distance: {} delta_y: {}".format(distance, delta_y))
                distance -= delta_y
                # move the next particle (particle being hit)
                particle = particle_hit
                position = index
                tile = grid[neigh_x][neigh_y]

                if (Values.pause):
                    input()

def build_histogram(rx): #compute all particle-particle distances and average
    temp_list = [] # this list will hold distances
    total = 0
    for i in range(0, len(rx) -1):
        for j in range(i+1, len(rx)):
            temp_list.append(rx[i]- rx[j])

    for i in range(0, len(temp_list) - 1):
        rlower = (i + 1) * 0.01
        rupper = rlower + 0.01
        volshell = (rupper - rlower)
        temp = (0.01 * (2.0 * i + 1.0) / 2.0) * Values.volume_length / volshell   
        Values.histogram.append(temp)

def write_to_file(file_name, rx, rx_init):
    output_file = open(file_name, 'w')
    
    output_file.write(
    """     System parameters
    -------------------------------------

    Number of particles:        {0:>7.4f}
    Number of Cycles:           {1:>7.4f}
    System time (seconds):      {2:>7.4f}
    Rho:                        {3:>7.4f}
    Volume Length:              {4:>7.4f}
    Lemgth:                     {5:>7.4f}

    -------------------------------------
    
    Initial      Final  

    """.format(
            Values.nparticles, Values.ncycles, Values.time,
            Values.rho, Values.volume_length, Values.nparticles ** 0.5))
    
    for i in range(len(rx)):
        output_file.write("{0:>7.4f} {1:>7.4f} \n".format(rx_init[i], rx[i]))

    if (len(Values.histogram) > 1): #print the histrogram if it has values
        for i in range(len(Values.histogram) - 1):
            ravg = 0.01 * (2.0 * i + 1.0) / 2.0
            output_file.write("ravg: {0:7.4f}\thistogram: {1:7.4f}\n".format(
                                ravg, Values.histogram[i]))


if __name__ == "__main__":
    main()
