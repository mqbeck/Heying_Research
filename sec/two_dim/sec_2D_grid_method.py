'''
Possible cause: not implementing move function correctly: update particle,
    check collisions, but particle new and particle old have the same x or y 
    because collision is called before particle new is moved and the smalled 
    distance change is the (potentially) overlapping new and old particles
Solution: particle new = particle old +/- some distance (radius?) to prevent
    overlap

same Problem different words: particle new == particle old
possible cause: move function updating new particle improperly
    updates particle to have new x or y position, then checks for collisions,
    but the new particle and old particle have the same x or y 
    coordinate so the change in that direction (x or y) is zero
solution: properly implement a move funtion with a correct distance
another solution: ignore particles that won't collide because they're in the
    opposite direction of movement (i.e. for up: if ynew - yold > 0 ignore)
'''


from math import *
import random
import time
import copy

class Values:
    diameter = 1  # particle diameter (for hard sphere model)
    ncycles = 200 # number of times a particle will be selceted and moved
    nparticles = 17 # number of particles
    # particle density - how close the particles will be at simulation start
    rho = 0.5 
    volume_length = (nparticles / rho) ** 0.5 # length of the box
    print('length: ', volume_length)
    # inverse of the length for use in calculations later
    inv_length = 1 / volume_length 
    length = nparticles ** 0.5 # length the particles will be moved per cycle
    n_root = ceil(sqrt(nparticles)) 
    # set the height of the grid_tils relative to the particles diameter
    grid_tile_height = (2 * diameter) 
    #particle-particle distances to be sorted and counted, per cycle
    distances = [] 
    #sorted distance, ultimately normalized (summed to one)
    histogram = [] 
    # how many cycles to run before taking a histogram measurement
    cycles_to_average = 2 
    time # the time is takes to run the simulation 
         # doesn't include minor set up and file printing

class Grid_tile:
    # declare a class to make object called grid tiles which have the 
    # following attributes a list of neighboring grids, a list of particles 
    # within the grid, and a position
    
    def __init__(self, i, j):
    # constructor method to build an instance of the class called an object

        self.neighbor_list = [[None, None] for i in range(0,8)] 
        # down, up, left, right, up left, up right, down left, down right
        # pre-built empty to optimise memory allocation
        
        # x and y represent the row and column
        self.x = i #grid tile's x position at top left corner
        self.y = j #grid tile's y position at top left corner
        #maybe better to define this as the center

        self.particle_list = [] 
        # create an empty list which will later be changed as the grid
        # is populated with particles
    
        # list.append(value) adds value to the list as the last element
        # list.pop(i) removes and returns the element at position i, 
        # default to the last element

    def __repr__(self):
        return '[{},{}]'.format(self.x, self.y)

def main():
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

    #Values.rx = define_list(Values.nparticles,Values.volume_length)
    random.seed(405) #seed the random number generator for result comparison 
    start = time.time() #record the starting time of the simulation
    
    #assign particles to grids, (populate grid tile particle lists)
    assign_particles_to_grid(rxy, grid)

    #moves random particle in random direction a set distance,
    #repeat for n cycles
    single_event_chain(rxy, grid) 

    #sorts resultant list for easy manipulation and data analysis
#    sort_list(Values.rxy) 

    end = time.time() # stops timer
    Values.time = (end - start)   #calculates 'total' time of experiment
    
    #constructs a histogram for chemical potential calculations
#    build_histogram(Values.rx)  

#    filename = "{0}_particles_{1}_cycles_{2}_minutes.txt".format(Values.nparticles, Values.ncycles, Values.time/60) #custom file name
#   write_to_file(filename, Values.rx, Values.rx_init) # write and output file
#    print(Values.rx) # print final, sorted list to terminal    #build_histogram(Values.rx)
  #  filename = "{0}_particles_{1}_cycles_{2}_minutes.txt".format(Values.nparticles, Values.ncycles, Values.time/60) # custom file name
  #  write_to_file(filename, Values.rx, Values.rx_init) # make a new file using custom file name
    rxy = sorted(rxy)
    print('\n\n\n')
    print("Start\tEnd")
    for i in range(0,len(start_rxy)):
        print("{}\t{}".format(start_rxy[i], rxy[i]))
#    print(rxy) # print final, sorted list to terminal for quick review without opening file

def build_grid(input_grid, grid_dim): 
    #O(n**2) to build the grid
    # populates the grid with grid tiles with dimensions 2*particle diameter 
    #(each grid can hold 4 particles)

    for i in range(0, len(input_grid[0])):
       for j in range(0, len(input_grid[0])):
           input_grid[i][j] = Grid_tile((j * grid_dim), (i * grid_dim))
        #endfor
    #endfor

def populate_grid_neighbors(input_grid): 
    # O(n) to populate nearest neighbor list

    grid_size = len(input_grid[0]) 
    for i in range(0, grid_size ** 2):
    # indices start at 0 and end at n-1
    # integer division and remainder to interate though rows and columns: 

        if (i == 0): # first grid_tile
            input_grid[i][i].neighbor_list = [
                    [grid_size - (grid_size - 1), 0], #down
                    [grid_size - 1, 0], # up
                    [0, grid_size - 1], # left
                    [0, grid_size - (grid_size - 1)], #  right
                    [grid_size - 1, grid_size - 1], #UpL
                    [grid_size - 1, grid_size - (grid_size - 1)], #upR
                    [grid_size - (grid_size - 1), grid_size - 1], #downL
                    [grid_size - (grid_size - 1), grid_size - (grid_size - 1)]] #downR
        #endif
        elif (i == ((grid_size ** 2) - 1)): # last grid_tile
            input_grid[grid_size - 1][grid_size - 1].neighbor_list = [
                    [0, grid_size - 1], #down
                    [grid_size - 2, grid_size - 1], # up
                    [grid_size - 1, grid_size - 2], # left
                    [grid_size - 1, 0], #right
                    [grid_size - 2, grid_size - 2], # upL
                    [grid_size - 2, 0], #upR
                    [0, grid_size - 2], # downL
                    [0, 0]] # downR
        #endelse
        elif (i == (grid_size - 1)): # top right corner
            input_grid[0][i].neighbor_list = [
                    [grid_size - (grid_size - 1), grid_size - 1], #D
                    [grid_size - 1, grid_size - 1], #U
                    [0, grid_size - 2], #L
                    [0, 0], #R
                    [grid_size - 1, grid_size - 2], #UL
                    [grid_size - 1, 0], #UR
                    [grid_size - (grid_size - 1), grid_size - 2], # DL
                    [grid_size - (grid_size - 1), 0]] #DR
        #endelse
        elif (i == ((grid_size ** 2) - grid_size)):# bottom left
            input_grid[grid_size - 1][0].neighbor_list = [
                    [0, 0], #D
                    [grid_size - 2, 0], #U
                    [grid_size - 1, grid_size - 1], #L
                    [grid_size - 1, grid_size - (grid_size - 1)], #R
                    [grid_size - (grid_size - 1), grid_size - 1], #UL
                    [grid_size - 2, grid_size - (grid_size - 1)], #UR
                    [0, grid_size - 1], #DL
                    [0, grid_size - (grid_size - 1)]] #DR
        #endelse            
        elif ((i // grid_size) == 0 ): # top row
            input_grid[i // grid_size][i % grid_size].neighbor_list = [
                    [grid_size - (grid_size - 1), i % grid_size], #D
                    [grid_size - 1, i % grid_size], #U
                    [0, (i - 1) % grid_size], #L
                    [0, (i + 1) % grid_size], #R
                    [grid_size - 1, (i - 1) % grid_size], #UL
                    [grid_size - 1, (i + 1) % grid_size], #UR
                    [grid_size - (grid_size - 1), (i - 1) % grid_size], #DL
                    [grid_size - (grid_size - 1), (i + 1) % grid_size]] #DR
        #endelse
#        elif ((grid_size ** 2 - grid_size) < i < ((grid_size ** 2) - 2)): # bottom row
        elif ((i // grid_size) == (grid_size - 1)):
            input_grid[i // grid_size][i % grid_size].neighbor_list = [
                    [0, i % grid_size], #D
                    [grid_size - 2, i % grid_size], #U
                    [grid_size - 1, (i - 1) % grid_size], #L
                    [grid_size - 1, (i + 1) % grid_size], #R
                    [grid_size - 2, (i - 1) % grid_size], #UL
                    [grid_size - 2, (i + 1) % grid_size], #UR
                    [0, (i - 1) % grid_size], #DL
                    [0, (i + 1) % grid_size]] #DR
        #endelse
        elif (i % grid_size == 0 and not(i // grid_size != 0 and
                i // grid_size == (grid_size - 1))): #left column
            input_grid[i // grid_size][i % grid_size].neighbor_list = [
                    [(i // grid_size) + 1, i % grid_size], #D
                    [(i // grid_size) - 1, i % grid_size], #U
                    [i // grid_size, grid_size - 1], #L
                    [i // grid_size, grid_size - (grid_size - 1)], #R
                    [(i // grid_size) - 1, grid_size - 1], #UL
                    [(i // grid_size) - 1, grid_size - (grid_size - 1)], #UR
                    [(i // grid_size) + 1, grid_size - 1], #DL
                    [(i // grid_size) + 1, grid_size - (grid_size - 1)]] #DR
        #endelse
        elif (i % grid_size == (grid_size - 1) and not(i // grid_size != 0 and
                i // grid_size == (grid_size - 1))): #right column
            input_grid[i // grid_size][i % grid_size].neighbor_list = [
                    [(i // grid_size) + 1, i % grid_size], #D
                    [(i // grid_size) - 1, i % grid_size], #U
                    [i // grid_size, grid_size - 2], #L
                    [i // grid_size, 0], #R
                    [(i // grid_size) - 1, grid_size - 2], #UL
                    [(i // grid_size) - 1, 0], #UR
                    [(i // grid_size) + 1, grid_size - 2], #DL
                    [(i // grid_size) + 1, 0]] #DR
        #endelse
        else: # everything else
            input_grid[i // grid_size][i % grid_size].neighbor_list = [
                    [(i // grid_size) + 1, i % grid_size], #D
                    [(i // grid_size) - 1, i % grid_size], #U
                    [i // grid_size, (i - 1) % grid_size], #L
                    [i // grid_size, (i + 1) % grid_size], #R
                    [(i // grid_size) - 1, (i - 1) % grid_size], #UL
                    [(i // grid_size) - 1, (i + 1) % grid_size], #UR
                    [(i // grid_size) + 1, (i - 1) % grid_size], #DL
                    [(i // grid_size) + 1, (i + 1) % grid_size]] #DR
        #endless
    #endfor

def define_list(nparticles, volume_length, diameter, rxy):
    # automatically, the particles are sorted in ascending order
    '''set this up on the grid using the grid size'''
    # one particle per grid space until no more particles

    # divide the particles by grid length and start filling the
    # grid (it doesn't matter how they start)
    
    for i in range(0, nparticles):
        rxy[i] = [i // Values.volume_length,
                    i % Values.volume_length]

    for i in rxy:
        print(i)
    #enddo    

def assign_particles_to_grid(rxy, grid): #O(n) or O(n ** dimension) 
    
    for i in range(0, len(rxy)):
        x = floor(rxy[i][0] / Values.grid_tile_height)
        y = floor(rxy[i][1] / Values.grid_tile_height)
        grid[x][y].particle_list.append(rxy[i])
        
def single_event_chain(rxy, grid):
    for i in range(0, Values.ncycles):
        particle, tile, position = select_particle_and_tile(grid)
#        direction = random.choice(['UP','DOWN','LEFT','RIGHT'])
        direction = "UP"
        move(direction, rxy, grid, particle, tile, position)
        if (i//10):
            print("cycle: ",i)
        
def select_particle_and_tile(grid):
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

def check_for_collisions(direction, grid, rand_particle, tile):
    x_distance_btwn = 0
    y_distance_btwn = 0
    overlap_count = 0
    next_particle = [[-9.0,[-9.0,-9.0]]] #dummy data to please the compiler
    neighbor_tiles = None
    neighbor_tile = None

    if (direction == 'UP'):
        # up, left, right, up-left, up-right
        neighbor_tiles = [1, 2, 3, 4, 5]

        #check self for collisions
        for particle in tile.particle_list:
            if (particle == rand_particle):
                continue

            x_distance_btwn = abs(rand_particle[0] - particle[0])

            if (particle[1] > rand_particle[1] or
                particle[1] < (rand_particle[1] - Values.grid_tile_height)):

                y_distance_btwn = abs(rand_particle[1] - particle[1])

            # if collisions exist, record them
            if (x_distance_btwn < Values.grid_tile_height):
                if (y_distance_btwn < Values.grid_tile_height):
                    overlap_count += 1
                    next_particle.append([sqrt(x_distance_btwn**2 +
                                        y_distance_btwn**2), particle,
                                        [tile.y//Values.grid_tile_height,
                                        tile.x//Values.grid_tile_height],
                                        x_distance_btwn])
                    return overlap_count, sorted(next_particle)

        #check neighbors for collisions
        for i in neighbor_tiles:
            # ith neighbor in the neighbor list
            neighbor_tile = tile.neighbor_list[i]
            neighbor = grid[neighbor_tile[0]][neighbor_tile[1]]
            for particle in neighbor.particle_list:
                if (particle == rand_particle):
                    continue
                x_distance_btwn = abs(rand_particle[0] - particle[0])
                # dont check particle behind rand particle
                if (particle[1] > rand_particle[1] or
                    particle[1] < (rand_particle[1] -
                    Values.grid_tile_height)):

                    y_distance_btwn = abs(rand_particle[1] - particle[1])

                # if collisions exist, record them
                if (x_distance_btwn < Values.grid_tile_height):
                    if (y_distance_btwn < Values.grid_tile_height):
                        overlap_count += 1
                        next_particle.append([sqrt(x_distance_btwn**2 +
                                            y_distance_btwn**2), particle,
                                            tile.neighbor_list[i],
                                            x_distance_btwn])
    elif (direction == 'DOWN'):
        # up, left, right, up-left, up-right
        neighbor_tiles = [0, 2, 3, 6, 7]
        #check neighbors for collisions
        for i in neighbor_tiles:
            # ith neighbor in the neighbor list
            neighbor_tile = tile.neighbor_list[i]
            neighbor = grid[neighbor_tile[0]][neighbor_tile[1]]
            for particle in neighbor.particle_list:
                x_distance_btwn = abs(rand_particle[0] - particle[0])
                if (particle[1] < rand_particle[1] or
                    particle[1] > (rand_particle[1] - Values.grid_tile_height)):
                    y_distance_btwn = abs(rand_particle[1] - particle[1])
                # if collisions exist, record them
                if (x_distance_btwn < Values.diameter):
                    if (y_distance_btwn < Values.diameter):
                        overlap_count += 1
                        next_particle.append([sqrt(x_distance_btwn**2 +
                                            y_distance_btwn**2), particle])
            #check self for collisions
            for particle in tile.particle_list:
                x_distance_btwn = abs(rand_particle[0] - particle[0])
                if (particle[1] < rand_particle[1] or
                    particle[1] > (rand_particle[1] - Values.grid_tile_height)):
                    y_distance_btwn = abs(rand_particle[1] - particle[1])
                # if collisions exist, record them
                if (x_distance_btwn < Values.diameter):
                    if (y_distance_btwn < Values.diameter):
                        overlap_count += 1
                        next_particle.append([sqrt(x_distance_btwn**2 +
                                            y_distance_btwn**2), particle])
    elif (direction == 'LEFT'):
        # up, left, right, up-left, up-right
        neighbor_tiles = [0, 1, 2, 4, 6]
        #check neighbors for collisions
        for i in neighbor_tiles:
            # ith neighbor in the neighbor list
            neighbor_tile = tile.neighbor_list[i]
            neighbor = grid[neighbor_tile[0]][neighbor_tile[1]]
            for particle in neighbor.particle_list:
                if (particle[0] < rand_particle[0] or
                    particle[0] > (rand_particle[0] - Values.grid_tile_height)):
                    x_distance_btwn = abs(rand_particle[0] - particle[0])
                y_distance_btwn = abs(rand_particle[1] - particle[1])
                # if collisions exist, record them
                if (x_distance_btwn < Values.diameter):
                    if (y_distance_btwn < Values.diameter):
                        overlap_count += 1
                        next_particle.append([sqrt(x_distance_btwn**2 +
                                            y_distance_btwn**2), particle])
            #check self for collisions
            for particle in tile.particle_list:
                if (particle[0] < rand_particle[0] or
                    particle[0] > (rand_particle[0] - Values.grid_tile_height)):
                    x_distance_btwn = abs(rand_particle[0] - particle[0])
                y_distance_btwn = abs(rand_particle[1] - particle[1])
                # if collisions exist, record them
                if (x_distance_btwn < Values.diameter):
                    if (y_distance_btwn < Values.diameter):
                        overlap_count += 1
                        next_particle.append([sqrt(x_distance_btwn**2 +
                                            y_distance_btwn**2), particle])
    elif (direction == 'RIGHT'):
        # up, left, right, up-left, up-right
        neighbor_tiles = [0, 1, 3, 5, 7]
        #check neighbors for collisions
        for i in neighbor_tiles:
            # ith neighbor in the neighbor list
            neighbor_tile = tile.neighbor_list[i]
            neighbor = grid[neighbor_tile[0]][neighbor_tile[1]]
            for particle in neighbor.particle_list:
                if (particle[0] > rand_particle[0] or
                    particle[0] < (rand_particle[0] - Values.grid_tile_height)):
                    x_distance_btwn = abs(rand_particle[0] - particle[0])
                y_distance_btwn = abs(rand_particle[1] - particle[1])
                # if collisions exist, record them
                if (x_distance_btwn < Values.diameter):
                    if (y_distance_btwn < Values.diameter):
                        overlap_count += 1
                        next_particle.append([sqrt(x_distance_btwn**2 +
                                            y_distance_btwn**2), particle])
            #check self for collisions
            for particle in tile.particle_list:
                if (particle[0] > rand_particle[0] or
                    particle[0] < (rand_particle[0] - Values.grid_tile_height)):
                    x_distance_btwn = abs(rand_particle[0] - particle[0])
                y_distance_btwn = abs(rand_particle[1] - particle[1])
                # if collisions exist, record them
                if (x_distance_btwn < Values.diameter):
                    if (y_distance_btwn < Values.diameter):
                        overlap_count += 1
                        next_particle.append([sqrt(x_distance_btwn**2 +
                                            y_distance_btwn**2), particle])
    else:
        pass

    # returns the number of overlaps, the list of overlapping particles,
    # and the tile in which the particles were found
    return overlap_count, sorted(next_particle)

def move(direction, rxy, grid, particle, tile, position):
    '''Legacy Comment'''
    #Rather than adding up, I decided to subtract - am I pessimistic?
    distance = Values.length
    x_distance_btwn = None
    y_distance_btwn = None
    neighbor_tile = None
    print("list{}: {} position: {}".format([tile.x, tile.y], tile.particle_list, position))
    for g_row in grid:
        for g_column in g_row:
            print('list{}: {}'.format([g_column], g_column.particle_list))

    if direction == 'UP':
#        while (distance > 0):
        for counter in range(11):
            print("distance: ",distance)
            #check particle diameter for collision  with grid above and left/right
            #move particle from grid list below to grid list above
            #remove particle from original list
            # update particle position and movement

            # check for collisions
            overlap_count, next_particle = check_for_collisions(
                direction, grid, particle, tile)
            # no collisions
            # next_particle should = None
            if (overlap_count == 0):
                print("particles: {} position: {}".format(tile.particle_list,
                                                            position))
                # remove particle from starting grid
                temp = tile.particle_list.pop(position)
                # update particle (moving up)
                # perhaps this should be minus: look into how the grid is implemented
                temp[1] += Values.grid_tile_height
                # determine destination grid (up is second in the list)
                destination = (grid[tile.neighbor_list[1][0]]
                                 [tile.neighbor_list[1][1]])
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
                distance_to_particle, particle_hit, neighbor_tile, delta_x = next_particle[1]
                print("particle_hit: {} particle moving: {}".format(
                        particle_hit, particle))
                #the index of the particle that will be hit
                neigh_x = neighbor_tile[0]
                neigh_y = neighbor_tile[1]
                print("particle list{}: {}".format([neigh_x, neigh_y], grid[neigh_x][neigh_y].particle_list))
                index = grid[neigh_x][neigh_y].particle_list.index(
                    particle_hit)

                # update particle (moving up)
                # perhaps this should be minus: look into how the grid is implemented
                delta_y = (abs(sqrt(Values.diameter - delta_x**2)))
                particle[1] += delta_y
                # update cycle movement
                print("distance: {} delta_y: {}".format(distance, delta_y))
                distance -= delta_y
                # move the next particle (particle being hit)
                print('Counter: {0}, Hit: {1}, Start: {2}, Delta: {3}'.format(
                        counter, particle_hit, particle, delta_y))
                particle = particle_hit
                position = index
                tile = grid[neigh_x][neigh_y]



    if direction == 'DOWN':
        while (distance > 0):
            print(distance)
            #check particle diameter for collision  with grid below and left/right
            #move particle from grid list above to grid list below
            #remove particle from original list
            # update particle position and movement

            # check for collisions
            overlap_count, next_particle = check_for_collisions(direction, grid,
                                                                particle, tile)
            # no collisions
            if (overlap_count == 0):
                # remove particle from starting grid
                temp = tile.particle_list.pop(position)
                # update particle (moving down)
                # perhaps this should be plus:  look into how the grid is implemented
                temp[1] -= Values.grid_tile_height
                # determine destination grid (down is first in the list)
                destination = tile.neighbor_list[0]
                # put particle in destination grid
                destination.particle_list.append(temp)
                # update cycle movement
                distance -=  Values.grid_tile_height
                # move the next particle (last particle in the new list)
                particle = destination.particle_list[-1]
                position = len(destination.particle_list[-1])

            # collision
            if (overlap_count > 0):
                # the closest collision
                distance_to_particle, particle_being_hit = next_particle[1]
                #the index of the particle that will be hit
                index = tile.particle_list.index(particle_being_hit)

                # update particle (moving down)
                # perhaps this should be plus:  look into how the grid is implemented
                delta_y = abs(particle_being_hit[1] - particle[1])
                particle[1] -= delta_y
                # update cycle movement
                distance -= delta_y
                # move the next particle (particle being hit)
                particle = particle_being_hit
                position = index

    if direction == 'LEFT':
        while (distance > 0):
            print(distance)
            #check particle diameter for collision  with grid up, down, and left
            #move particle from to grid list to the left
            #remove particle from original list
            # update particle position and movement

            # check for collisions
            overlap_count, next_particle = check_for_collisions(direction, grid,
                                                                particle, tile)
            # no collisions
            if (overlap_count == 0):
                # remove particle from starting grid
                temp = tile.particle_list.pop(position)
                # update particle (moving down)
                # perhaps this should be minus: look into how the grid is implemented
                temp[0] -= Values.grid_tile_height
                # determine destination grid (left is third in the list)
                destination = tile.neighbor_list[2]
                # put particle in destination grid
                destination.particle_list.append(temp)
                # update cycle movement
                distance -=  Values.grid_tile_height
                # move the next particle (last particle in the new list)
                particle = destination.particle_list[-1]
                position = len(destination.particle_list[-1])

            # collision
            if (overlap_count > 0):
                # the closest collision
                distance_to_particle, particle_being_hit = next_particle[1]
                #the index of the particle that will be hit
                index = tile.particle_list.index(particle_being_hit)

                # update particle (moving left)
                # perhaps this should be minus: look into how the grid is implemented
                delta_x = abs(particle_being_hit[0] - particle[0])
                particle[0] += delta_x
                # update cycle movement
                distance -= delta_x
                # move the next particle (particle being hit)
                particle = particle_being_hit
                position = index

    if direction == 'RIGHT':
        while (distance > 0):
            print(distance)
            #check particle diameter for collision with grid below, above, right
            #move particle from left to right
            #remove particle from original list
            # update particle position and movement

            # check for collisions
            overlap_count, next_particle = check_for_collisions(direction, grid,
                                                                particle, tile)
            # no collisions
            if (overlap_count == 0):
                # remove particle from starting grid
                temp = tile.particle_list.pop(position)
                # update particle (moving down)
                # perhaps this should be plus:  look into how the grid is implemented
                temp[0]-= Values.grid_tile_height
                # determine destination grid (right is fourth in the list)
                destination = tile.neighbor_list[3]
                # put particle in destination grid
                destination.particle_list.append(temp)
                # update cycle movement
                distance -=  Values.grid_tile_height
                # move the next particle (last particle in the new list)
                particle = destination.particle_list[-1]
                position = len(destination.particle_list[-1])

            # collision
            if (overlap_count > 0):
                # the closest collision
                distance_to_particle, particle_being_hit = next_particle[1]
                #the index of the particle that will be hit
                index = tile.particle_list.index(particle_being_hit)

                # update particle (moving down_
                # perhaps this should be plus:  look into how the grid is implemented
                delta_x = abs(particle_being_hit[0] - particle[0])
                particle[0] -= delta_x
                # update cycle movement
                distance -= delta_x
                # move the next particle (particle being hit)
                particle = particle_being_hit
                position = index

def sort_list(rx): # O(nlogn)
    # sort the list, rx, for simpler implementation 
    # of Single Event Chain in one-dimension
    
    # under the advisement of of Joe DeLuca this sort will
    # follow the merge sort algorithm until n <= 20 at 
    # which point it will become an insertion sort
    
    if len(rx) < 20:
        insertion_sort(rx) 
    else:
        middle = len(rx) // 2
        lefthalf = rx[:middle]  # this is Python's "slice" notation
        righthalf = rx[middle:] # string[start:stop:step]

        sort_list(lefthalf)
        sort_list(righthalf)

        i = 0
        j = 0
        k = 0
        while i < len(lefthalf) and j < len(righthalf):
            if lefthalf[i] < righthalf[j]:
                rx[k] = lefthalf[i]
                i = i + 1
            else:
                rx[k]=righthalf[j]
                j = j + 1
            k = k + 1

        while i < len(lefthalf):
            rx[k] = lefthalf[i]
            i = i + 1
            k = k + 1

        while j < len(righthalf):
            rx[k] = righthalf[j]
            j = j + 1
            k = k + 1
    Values.rx = rx

def insertion_sort(rx):
    position = None
    value = None
               
    for i in range(1, len(rx)):
        value = rx[i]
        position = i
                
        while position > 0 and rx[position - 1] > value:
            rx[position] = rx[position - 1]       
            position -= 1
        rx[position] = value
        
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
"""System parameters

Number of particles:        {0:>7.4f}
Number of Cycles:           {1:>7.4f}
System time (seconds):      {2:>7.4f}
Rho:                        {3:>7.4f}
Volume Length:              {4:>7.4f}
Lemgth:                     {5:>7.4f}

Initial      Final  

""".format(
            Values.nparticles, Values.ncycles, Values.time,
            Values.rho, Values.volume_length, Values.nparticles ** 0.5))
    
    for i in range(len(rx)):
        output_file.write("{0:>7.4f} {1:>7.4f} \n".format(rx_init[i], rx[i]))

    if (len(Values.histogram) > 1): #only print the histrogram if it has values
        for i in range(len(Values.histogram) - 1):
            ravg = 0.01 * (2.0 * i + 1.0) / 2.0
            output_file.write("ravg: {0:7.4f}    histogram: {1:7.4f}\n".format(ravg, Values.histogram[i]))


if __name__ == "__main__":
    main()
