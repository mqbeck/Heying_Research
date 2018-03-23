from math import *
import random
import time
import copy

class Values:
    diameter = 1  # particle diameter (for hard sphere model)
    ncycles = 2000 # number of times a particle will be selceted and moved
    nparticles = 17 # number of particles
    rho = 0.5 # particle density - how close the particles will be at simulation start
    volume_length = (nparticles / rho) ** 0.5 # length of the box
    inv_length = 1 / volume_length # inverse of the length for use in calculations later
    length = nparticles ** 0.5 # length the particles will be moved per cycle
    rx = [None] * nparticles # initial, empty list
    ry = [None] * nparticles 
    n_root = ceil(sqrt(nparticles)) 
    grid_tile_height = (2 * diameter) # set the height of the grid_tils relative to the particles diameter
    distances = [] #particle-particle distances to be sorted and counted, per cycle
    histogram = [] #sorted distance, ultimately normalized (summed to one)
    cycles_to_average = 2 # how many cycles to run before taking a histogram measurement
    time  # the time is takes to run the simulation - doesn't include minor set up and file printing

class Grid_tile:
    # declare a class to make object called grid tiles which have the following attributes
    # a list of neighboring grids, a list of particles within the grid, and a position
    
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

def main():
    #construct an empty square grid
    grid = [[None for column in range(0, Values.n_root)] 
            for row in range(0, Values.n_root)] 
    #pre allocate the memory for the particles
    rxy = [[None,None] for xy in range(0, Values.n_root)]

    #construct a grid that is ceiling(sqrt(n)) x ceiling(sqrt(n)) (int x int)
    build_grid(grid, Values.grid_tile_height) 

    #tell each grid which grids it is next to (non-bound box i.e.wrap around)
    populate_grid_neighbors(grid) 

    #Values.rx = define_list(Values.nparticles,Values.volume_length)
    random.seed(405) #seed the random number generator for result comparison 
    start = time.time() #record the starting time of the simulation
    
    #assign particles to grids, (populate grid tile particle lists)
    assign_particles_to_grid(rxy, grid)
    for i in grid:
        print(i)

    #moves random particle in random direction a set distance, repeat for n cycles
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
    print('\n\n\n')
    print(rxy) # print final, sorted list to terminal for quick review without opening file

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

def define_list(nparticles, volume_length, diameter,rx, ry, rxy):
    # automatically, the particles are sorted in ascending order
    vector = volume_length / nparticles - diameter
    xpos = 0
    for i in range(0, nparticles, ceil(nparticles / 10)):
        ypos = 0
        for j in range(i, i + ceil(nparticles / 10)):
            xpos += random.random() * vector
            rx[j] = xpos
            xpos += diameter * (i // (volume_length))

            ypos += random.random() * vector
            ry[j] = ypos
            ypos += diameter
        
            rxy[j] = [rx[j],ry[j]]
    #enddo    

def assign_particles_to_grid(rxy, grid): #O(n) or O(n ** dimension) 
    
    for i in range(0, len(rxy)):
        x = floor(rxy[i][0] // Values.grid_tile_height)
        y = floor(rxy[i][1] // Values.grid_tile_height)
        grid[x][y].particle_list.append(rxy[i])
        
def single_event_chain(rxy, grid):
    for i in range(0, Values.ncycles):
        direction = random.choice(['UP','DOWN','LEFT','RIGHT'])
        move(direction, rxy)
        
def move(direction, rxy, grid):
    #Rather than adding up, I decided to subtract - am I pessimistic?
    distance = Values.length
    
    the_chosen_one = random.randint(0, grid)
    random_number = random.randint(0, len(the_chosen_one.particle_list))
    rand_particle = the_chosen_one.particle_list(random_number)
    overlap_count = 0
    if direction == 'UP':
            #check particle diameter for collision  with grid above and left/right
            #move particle from grid list below to grid list above
            #remove partilce from original list
            for i in range(1,6):
                neighbor = the_chosen_one.neighbor_list[i]
                for particle in neighbor.particle_list:
                    if(abs(rand_particle[0] - particle[0]) < Values.diameter):
                        if(abs(rand_particle[1] - particle[1]) < Values.diameter):
                            overlap_count += 1
            if (overlap_count == 0):
                temp = the_chosen_one.particle_list.pop(random_number)
                the_chosen_one.neighbor_list[1].particle_list.append(temp)

    if direction == 'UP':
            #check particle diameter for collision  with grid above and left/right
            #move particle from grid list below to grid list above
            #remove partilce from original list
            for i in [1, 2, 3, 4, 5]:
                neighbor = the_chosen_one.neighbor_list[i]
                for particle in neighbor.particle_list:
                    if(abs(rand_particle[0] - particle[0]) < Values.diameter):
                        if(abs(rand_particle[1] - particle[1]) < Values.diameter):
                            overlap_count += 1
            if (overlap_count == 0):
                temp = the_chosen_one.particle_list.pop(random_number)
                the_chosen_one.neighbor_list[1].particle_list.append(temp)

    if direction == 'DOWN':
            #check particle diameter for collision  with grid above and left/right
            #move particle from grid list below to grid list above
            #remove partilce from original list
            for i in [0, 2, 3, 6, 7]:
                neighbor = the_chosen_one.neighbor_list[i]
                for particle in neighbor.particle_list:
                    if(abs(rand_particle[0] - particle[0]) < Values.diameter):
                        if(abs(rand_particle[1] - particle[1]) < Values.diameter):
                            overlap_count += 1
            if (overlap_count == 0):
                temp = the_chosen_one.particle_list.pop(random_number)
                the_chosen_one.neighbor_list[1].particle_list.append(temp)

    if direction == 'LEFT':
            #check particle diameter for collision  with grid above and left/right
            #move particle from grid list below to grid list above
            #remove partilce from original list
            for i in [0, 1, 2, 4, 6]:
                neighbor = the_chosen_one.neighbor_list[i]
                for particle in neighbor.particle_list:
                    if(abs(rand_particle[0] - particle[0]) < Values.diameter):
                        if(abs(rand_particle[1] - particle[1]) < Values.diameter):
                            overlap_count += 1
            if (overlap_count == 0):
                temp = the_chosen_one.particle_list.pop(random_number)
                the_chosen_one.neighbor_list[1].particle_list.append(temp)

    if direction == 'RIGHT':
            #check particle diameter for collision  with grid above and left/right
            #move particle from grid list below to grid list above
            #remove partilce from original list
            for i in [0, 1, 3, 5, 7]:
                neighbor = the_chosen_one.neighbor_list[i]
                for particle in neighbor.particle_list:
                    if(abs(rand_particle[0] - particle[0]) < Values.diameter):
                        if(abs(rand_particle[1] - particle[1]) < Values.diameter):
                            overlap_count += 1
            if (overlap_count == 0):
                temp = the_chosen_one.particle_list.pop(random_number)
                the_chosen_one.neighbor_list[1].particle_list.append(temp)

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
