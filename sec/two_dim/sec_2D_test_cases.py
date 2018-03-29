#test file for sec_2D_grid_method.py

from sec_2D_grid_method import *

#test the class constructor for Grid_tile
test_tile = Grid_tile(1,2) #build an object of the class Grid_tile called test_tile
assert(test_tile.x == 1) #test the x value of test_tile
assert(test_tile.y == 2) #stops program if false otherwise nothing is printed to termnial
assert(test_tile.neighbor_list == [[None, None],[None, None],[None, None],[None, None],
                   [None, None],[None, None],[None, None],[None, None]])

# test build_grid(input_grid, grid_dim) function
test_grid = [[None for column in range(0,3)] for row in range(0,3)] #create an empty 3x3 grid
        
build_grid(test_grid, 2)
assert(test_grid[0][0].x == 0)
assert(test_grid[2][2].x == 4)
assert(test_grid[0][0].y == 0)
assert(test_grid[2][2].y == 4)

# test cases for the populate_grid(input_grid) function
# populates a nearest neighbor list
populate_grid_neighbors(test_grid)

assert(test_grid[0][0].neighbor_list ==
					[[1,0],[2,0],[0,2],[0,1],[2,2],[2,1],[1,2],[1,1]])
assert(test_grid[0][1].neighbor_list ==
					[[1,1],[2,1],[0,0],[0,2],[2,0],[2,2],[1,0],[1,2]])
assert(test_grid[0][2].neighbor_list ==
					[[1,2],[2,2],[0,1],[0,0],[2,1],[2,0],[1,1],[1,0]])
assert(test_grid[1][0].neighbor_list ==
					[[2,0],[0,0],[1,2],[1,1],[0,2],[0,1],[2,2],[2,1]])
assert(test_grid[1][1].neighbor_list ==
					[[2,1],[0,1],[1,0],[1,2],[0,0],[0,2],[2,0],[2,2]])
assert(test_grid[1][2].neighbor_list ==
					[[2,2],[0,2],[1,1],[1,0],[0,1],[0,0],[2,1],[2,0]])
assert(test_grid[2][0].neighbor_list ==
					[[0,0],[1,0],[2,2],[2,1],[1,2],[1,1],[0,2],[0,1]])
assert(test_grid[2][1].neighbor_list ==
					[[0,1],[1,1],[2,0],[2,2],[1,0],[1,2],[0,0],[0,2]])
assert(test_grid[2][2].neighbor_list == 
					[[0,2],[1,2],[2,1],[2,0],[1,1],[1,0],[0,1],[0,0]])

nparticles = 14
volume_length = (nparticles / 0.5) 
diameter = 1
rx = [None] * nparticles
ry = [None] * nparticles
rxy = [[[None],[None]] for xy in range(0, nparticles)]
define_list(nparticles, volume_length, diameter, rx, ry, rxy)
assert(rxy[0][0] < rxy[1][0]) #first x is smaller than second x
assert(rxy[0][1] < rxy[1][1]) #first y is smalller than second y
#for i in rxy:
#    print(i)

# test cases for the assign_particles_to_grid(rxy, grid) function
test_particles = [[0,0],[0,1],[1,0],[2,1]]
assign_particles_to_grid(test_particles, test_grid)
assert(test_particles[0] in test_grid[0][0].particle_list)
assert(test_particles[1] in test_grid[0][0].particle_list) 
assert(test_particles[2] in test_grid[0][0].particle_list)
assert(test_particles[3] in test_grid[1][0].particle_list)


