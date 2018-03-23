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
    rxy_init = None # initial value of particles, needed to make python happy
    distances = [] #particle-particle distances to be sorted and counted, per cycle
    histogram = [] #sorted distance, ultimately normalized (summed to one)
    cycles_to_average = 2 # how many cycles to run before taking a histogram measurement
    time  # the time is takes to run the simulation - doesn't include minor set up and file printing

class Grid_tile:
    # declare a class to make object called grid tiles which have the following attributes
    # a list of neighboring grids, a list of particles within the grid, and a position
    

    def __init__(self, i, j):
    # constructor method to build an instance of the class called an object
                                             # use: x = Grid_tile(i,j)

        self.neighbor_list = [[None, None] for i in range(0,8)] 
        # down, up, left, right, up left, up right, down left, down right
        # pre-built empty to optimise memory allocation
        
        # x and y represent the row and column
        self.x = i #grid tile's x position at top left corner
        self.y = j #grid tile's y position at top left corner
        #maybe better to define this as the center

        particle_list = [] 
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
    rxy = [[None,None] for xy in range(0, Values.nroot)]

    #construct a grid that is ceiling(sqrt(n)) x ceiling(sqrt(n)) (int x int)
    build_grid(grid, Values.grid_tile_height) 

    #tell each grid which grids it is next to (non-bound box i.e.wrap around)
    populate_grid_neighbors(grid) 

    #Values.rx = define_list(Values.nparticles,Values.volume_length)
    random.seed(405) #seed the random number generator for result comparison 
    start = time.time() #record the starting time of the simulation
    



     #define x_position of particles
#    define_list_x(Values.nparticles,Values.volume_length,Values.diameter,Values.rx)
    #define y_position of particles
#    define_list_y(Values.nparticles,Values.volume_length,Values.diameter,Values.ry)

    #take a list of particles, rxy, and based on location in grid builds a
    #list of particles in that grid
 #   assign_particles_to_grid(rxy, grid)

    #moves random particle in random direction a set distance, repeat for n cycles
#     single_event_chain() 

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
  #  print(Values.rx) # print final, sorted list to terminal for quick review without opening file

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
        x = floor(rxy[i][0] + Values.diameter)
        y = floor(rxy[i][1] + Values.diameter)
        grid[x][y].particle_list.append(rxy[i])
    #psuedo code:
 #   for particle in rxy:
 #       if particle[x] is within grid[x] +  2*Values.diameter and particle[y] is within grid[y] + 2*diamater
 #       then grid.particle_list.append(particle)
 #       
#
#    # maybe use some sort of division:
#        x = (n_root * 2) / particle[x]
#        grid[x]/particle_list.append(particle)

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
        
def single_event_chain():
    for i in range(0, Values.ncycles):
        direction = random.choice([-1, 1])
        move(direction)
        
def move(direction):
    rx = Values.rx
    #Rather than adding up, I decided to subtract - am I pessimistic?
    distance = Values.length
    
    particle = random.randint(0, Values.nparticles - 1)
    if direction == 1:
        for i in range(particle, particle + Values.nparticles, direction):
            #I think that you need to check whether the particles
            #are closer within the box or within the PBCs
            j = i
            k = i + 1
            if (j > (Values.nparticles - 1)):
                j -= Values.nparticles
            if (k > (Values.nparticles - 1)):
                k -= Values.nparticles
            #endif
            vector = abs(rx[k] - rx[j])
            if rx[j] > rx[k]:
                traveled = Values.volume_length - vector - Values.diameter
            #endif
            #vector -= Values.volume_length*round(vector*Values.inv_length)
            else:
                traveled = abs(rx[k] - rx[j]) - Values.diameter
            
            if (distance > traveled):
                rx[j] += traveled
                distance -= traveled
            else:
                rx[j] += distance
                distance = 0
            #endif
            
            #Since you are moving them to the right, you should
            #not have to check for the rx[i] < volume_length condition
            if rx[j] > Values.volume_length:
                rx[j] -= Values.volume_length
            #endif
            if distance == 0:
                break
            #endif
        #enddo
    elif direction == -1: # negative ouptut implies problem here
        for i in range(particle, particle - Values.nparticles, direction):
            j = i
            k = i - 1
       #     print("j: ", j, "k: ", k)
            
            if (j < 0):
                j += Values.nparticles
            #endif
            if (k < 0):
                k += Values.nparticles
            #endif
            
            vector = abs(rx[k] - rx[j])
            if rx[j] < rx[k]:
                traveled = Values.volume_length - vector - Values.diameter
            #endif
            #vector -= Values.volume_length*round(vector*Values.inv_length)
            else:
                traveled = abs(rx[k] - rx[j]) - Values.diameter
            #vector -= Values.volume_length*round(vector*Values.inv_length)
            #print("vector: ", vector, "traveled: ", traveled)
            
            if (distance > traveled):
                rx[j] -= traveled
                distance -= traveled
            else:
                rx[j] -= distance
                distance = 0
            #endif    
            
            if rx[j] < 0:
                rx[j] += Values.volume_length
            #endif
            if distance == 0:
                break
            #endif
        #enddo
    #endif

    Values.rx = rx

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
