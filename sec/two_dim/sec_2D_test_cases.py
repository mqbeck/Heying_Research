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
        
build_grid(test_grid, 1)
assert(test_grid[0][0].x == 0)
assert(test_grid[2][2].x == 2)
assert(test_grid[0][0].y == 0)
assert(test_grid[2][2].y == 2)

# test cases for the populate_grid(input_grid) function
# populates a nearest neighbor list
populate_grid_neighbors(test_grid)

# neighbor list order: down, up, left, right, up left, up right, down left, down right
assert(test_grid[0][0].neighbor_list ==
					[[0,2],[0,1],[2,0],[1,0],[1,2],[1,1],[2,2],[2,1]])
assert(test_grid[0][1].neighbor_list ==
					[[0,0],[0,2],[1,1],[2,1],[2,2],[1,2],[2,0],[1,0]])
assert(test_grid[0][2].neighbor_list ==
					[[0,1],[0,0],[2,2],[1,2],[2,0],[1,0],[2,1],[1,1]])
assert(test_grid[1][0].neighbor_list ==
					[[1,2],[1,1],[0,0],[2,0],[0,1],[2,1],[0,2],[2,2]])
assert(test_grid[1][1].neighbor_list ==
					[[1,0],[1,2],[0,1],[2,1],[0,2],[2,2],[0,0],[2,0]])
assert(test_grid[1][2].neighbor_list ==
					[[1,1],[1,0],[0,2],[2,2],[0,0],[2,0],[0,1],[2,1]])
assert(test_grid[2][0].neighbor_list ==
					[[2,2],[2,1],[1,0],[0,0],[1,1],[0,1],[1,2],[0,2]])
assert(test_grid[2][1].neighbor_list ==
					[[2,0],[2,2],[1,1],[0,1],[1,2],[0,2],[1,0],[0,0]])
assert(test_grid[2][2].neighbor_list == 
					[[2,1],[2,0],[1,2],[0,2],[1,0],[0,0],[1,1],[0,1]])

nparticles = 14
volume_length = (nparticles / 0.5) 
diameter = 1
rxy = [[[None],[None]] for xy in range(0, nparticles)]

# test cases for the define list (nparticles, volume length, diamete, rxy) func
define_list(nparticles, volume_length, diameter, rxy)

# rxy is an array, rxy[i] is an specfic particle, rxy[i][0] is the x position
# so the following is the distance between the first and second particle
assert(((rxy[0][0] - rxy[1][0])**2 + (rxy[0][1] - rxy[1][1])**2)**(0.5) >= diameter)
for i in range(0, len(rxy)-1):
    assert(((rxy[i][0] - rxy[i+1][0])**2 + (rxy[i][1] - rxy[i+1][1])**2)**(0.5) >= diameter)

# test cases for the assign_particles_to_grid(rxy, grid) function
test_particles = [[0,0],[0,1],[1,0],[2,1]]
assign_particles_to_grid(test_particles, test_grid)
assert(test_particles[0] in test_grid[0][0].particle_list)
assert(test_particles[1] in test_grid[0][1].particle_list) 
assert(test_particles[2] in test_grid[1][0].particle_list)
assert(test_particles[3] in test_grid[2][1].particle_list)

# Test cases for select_particle_and_tile(grid) function
chosen_particle, chosen_tile, position = select_particle_and_tile(test_grid)
assert(chosen_particle in chosen_tile.particle_list)
assert(chosen_tile.particle_list[position] == chosen_particle)


# test cases for the check_for_collisions(direction, grid, rand_particle, tile)

# test_grid with particles in list form:
# [0,0] : (0,0)
# [1,0] : (1,0)
# [2,0] : 
# [0,1] : (0,1)
# [1,1] : 
# [2,1] : (2,1)
# [0,2] :
# [1,2] : 
# [2,2] : 

direction = 'UP'
grid = test_grid







