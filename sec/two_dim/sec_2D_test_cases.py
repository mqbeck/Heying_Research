#test file for sec_2D_grid_method.py

import unittest

from sec_2D_grid_method import *

class TestGridMethod(unittest.TestCase):

    def test_grid_tile_constructor(self):
        '''verifies that the grid_tile class constructor method creates grid 
        tiles objects with x and y coordinates in the orer x,y and formats the
        neighbor list as 8 none 2-tuples'''

        #build an object of the class Grid_tile called test_tile
        test_tile = Grid_tile(1,2) 

        #test the x and y value of test_tile
        self.assertEqual(test_tile.x, 1) 
        self.assertEqual(test_tile.y, 2) 

        #tests the format of the neighbor list
        self.assertEqual(test_tile.neighbor_list, 
                          [[None, None],[None, None],[None, None],[None, None],
                          [None, None],[None, None],[None, None],[None, None]])

    def test_build_grid(self):
        '''verifies the build_grid function takes a two dimensional list of 
        none 2-tuples and converts it into a grid with x and y coordinates by 
        filling the grid with grid_tile objects'''

        #create an empty 3x3 grid and assign x and y coordinates to it
        test_grid = [[None for column in range(0,3)] for row in range(0,3)] 
        build_grid(test_grid, 1)

        # tests that the x and y are in the order x,y
        self.assertEqual(test_grid[0][0].x, 0)
        self.assertEqual(test_grid[0][0].y, 0)

        self.assertEqual(test_grid[0][1].x, 0)
        self.assertEqual(test_grid[0][1].y, 1)

        self.assertEqual(test_grid[1][0].x, 1)
        self.assertEqual(test_grid[1][0].y, 0)

        self.assertEqual(test_grid[2][2].x, 2)
        self.assertEqual(test_grid[2][2].y, 2)

    def test_populate_grid_neighbors(self):
        '''verifies that the populate_grid_neighbors function associates grid 
        tiles neighbors in the order: down, up, left, right, up left, up right,
        down left, down right'''

        #create an empty 3x3 grid and assign x and y coordinates to it
        test_grid = [[None for column in range(0,3)] for row in range(0,3)] 
        build_grid(test_grid, 1)
        
        # populates a nearest neighbor list
        populate_grid_neighbors(test_grid)

        # neighbor list order: down, up, left, right, up left, up right, down left, down right
        self.assertEqual(test_grid[0][0].neighbor_list,
                            [[0,2],[0,1],[2,0],[1,0],[2,1],[1,1],[2,2],[2,1]])
        self.assertEqual(test_grid[0][1].neighbor_list,
                            [[0,0],[0,2],[2,1],[1,1],[2,2],[1,2],[2,0],[1,0]])
        self.assertEqual(test_grid[0][2].neighbor_list,
                            [[0,1],[0,0],[2,2],[1,2],[2,0],[1,0],[2,1],[1,1]])
        self.assertEqual(test_grid[1][0].neighbor_list,
                            [[1,2],[1,1],[0,0],[2,0],[0,1],[2,1],[0,2],[2,2]])
        self.assertEqual(test_grid[1][1].neighbor_list,
                            [[1,0],[1,2],[0,1],[2,1],[0,2],[2,2],[0,0],[2,0]])
        self.assertEqual(test_grid[1][2].neighbor_list,
                            [[1,1],[1,0],[0,2],[2,2],[0,0],[2,0],[0,1],[2,1]])
        self.assertEqual(test_grid[2][0].neighbor_list,
                            [[2,2],[2,1],[1,0],[0,0],[1,1],[0,1],[1,2],[0,2]])
        self.assertEqual(test_grid[2][1].neighbor_list,
                            [[2,0],[2,2],[1,1],[0,1],[1,2],[0,2],[1,0],[0,0]])
        self.assertEqual(test_grid[2][2].neighbor_list, 
                            [[2,1],[2,0],[1,2],[0,2],[1,0],[0,0],[1,1],[0,1]])

    def test_define_list(self):
        '''verifies that: the  define_list function creates particles with 
        no overlap starting at the origin and moving from left to right, bottom
        to top'''

        # variables small enough to track by hand
        nparticles = 3
        volume_length = (nparticles / 0.5) 
        diameter = 1
        rxy = [[[None],[None]] for xy in range(0, nparticles)]

        # run the function
        define_list(nparticles, volume_length, diameter, rxy)

        # rxy is an array, rxy[i] is an specfic particle, rxy[i][0] is the x position
        # so the following is the distance between the first and second particle
        self.assertTrue(((rxy[0][0] - rxy[1][0])**2 + (rxy[0][1] -
                        rxy[1][1])**2)**(0.5) >= diameter)

        #distance between every particle is greater than diameter (no overlap)
        for i in range(0, len(rxy)-1):
            self.assertTrue(((rxy[i][0] - rxy[i+1][0])**2 + (rxy[i][1] -
                            rxy[i+1][1])**2)**(0.5) >= diameter)

    def test_assign_particles_to_grid(self):
        '''verifies that the assign_particles_to_grid function places the
        particles in the grid tiles so that the particle x and y are greater
        than or equal the tile x and y but less than the tile (x+1) and (y+1)'''

        #create an empty 3x3 grid and assign x and y coordinates to it
        test_grid = [[None for column in range(0,3)] for row in range(0,3)] 
        build_grid(test_grid, 1)

        # pre-constructed particles and grid
        test_particles = [[0.0,0.0],[0.0,1.0],[1.0,0.0],[2.0,1.0]]
        assign_particles_to_grid(test_particles, test_grid)

        self.assertTrue(test_particles[0] in test_grid[0][0].particle_list)
        self.assertTrue(test_particles[1] in test_grid[0][1].particle_list) 
        self.assertTrue(test_particles[2] in test_grid[1][0].particle_list)
        self.assertTrue(test_particles[3] in test_grid[2][1].particle_list)

    def test_select_particle_and_tile(self):
        '''verifies that a tile is selected and that a particle in the selected
        tile is selected. if the tile is empty then a new tile is selected'''

        grid = [[None, None], [None, None]]
        build_grid(grid)
        grid[0][0].particle_list.append([0,0])
        
        chosen_particle, chosen_tile, position = select_particle_and_tile(grid)

        self.assertTrue(chosen_particle in chosen_tile.particle_list)
        self.assertTrue(chosen_tile.particle_list[position] == chosen_particle)

        self.assertEqual(chosen_particle, [0,0])
        self.assertEqual(chosen_tile, grid[0][0]) 
        self.assertEqual(position, 0)

        # using the ultra-secret override command
        chosen_particle, chosen_tile, position = select_particle_and_tile(grid, True)
        self.assertEqual(chosen_particle, [0,0])
        self.assertEqual(chosen_tile, grid[0][0]) 
        self.assertEqual(position, 0)

    def test_check_for_collisions_small_grid(self):
        '''verifies that the check_for_collisions function detects potentional 
        overlap and records them for use in the move function. it first checks
        collisions in its own tile and if none are found checks the neighbors
        in the direction traveled'''
        
        #test cases designed with diamter = 1 for the pass/fail conditions
        Values.diameter = 1
        # smallest possible case (because of list indexes)
        # [1,0] | [1,1]
        # empty   empty
        # [0,0] | [1,0]
        #   x     empty
        # no overlap

        grid = [[None, None], [None, None]]
        build_grid(grid)

        populate_grid_neighbors(grid)
        grid[0][0].particle_list.append([0,0])

        chosen_particle, chosen_tile, position = select_particle_and_tile(grid)
        overlap_count, overlapping_particles = check_for_collisions('UP', grid,
                                                chosen_particle, chosen_tile,
                                                2, 0.5)

        # no collisions
        self.assertEqual(overlap_count, 0)
        self.assertEqual(overlapping_particles, [[-9.0,[-9.0,-9.0]]])

        # [0,1] | [1,1]
        # empty   empty
        # [0,0] | [1,0]
        #   x        x
        # two overlaps (left and right)
        grid[1][0].particle_list.append([1.2,0])
        overlap_count, overlapping_particles = check_for_collisions('UP', grid,
                                                chosen_particle, chosen_tile,
                                                2, 0.5)
        
        # a 2 implies that it's finding both left and right collisions
        # a 1 implies its finding up and thus everything is rotated
        self.assertEqual(overlap_count, 2)
        

    def test_check_for_collisions_large_grid_non_PBC(self):
        '''verifies that the check_for_collisions function detects potentional 
        overlap and records them for use in the move function. it first checks
        collisions in its own tile and if none are found checks the neighbors
        in the direction traveled'''

        #building a four by four grid
        grid = [[None, None, None, None],
                [None, None, None, None],
                [None, None, None, None],
                [None, None, None, None]]

        build_grid(grid)
        populate_grid_neighbors(grid)

        # collision path
        grid[0][0].particle_list.append([0,0])
        grid[0][1].particle_list.append([0,1])

        # select particle at origin and check for collision
        chosen_particle, chosen_tile, position = select_particle_and_tile(grid, True)
        overlap_count, overlapping_particles = check_for_collisions('UP', grid,
                                                chosen_particle, chosen_tile,
                                                4, 0.25)

        # one collision
        self.assertEqual(overlap_count, 1)
        self.assertNotEqual(overlapping_particles, [[-9.0,[-9.0,-9.0]]])


        # remove particle and put it in the tile above
        # and repeat the test
        grid[0][1].particle_list.pop()
        grid[0][2].particle_list.append([0,2])

        # select particle at origin and check for collision
        chosen_particle, chosen_tile, position = select_particle_and_tile(grid, True)
        overlap_count, overlapping_particles = check_for_collisions('UP', grid,
                                                chosen_particle, chosen_tile,
                                                4, 0.25)
        # no collisions
        self.assertEqual(overlap_count, 0)
        self.assertEqual(overlapping_particles, [[-9.0,[-9.0,-9.0]]])

        # remove particle and put it one tile over and some distance up
        # and repeat the test
        grid[0][2].particle_list.pop()
        grid[0][1].particle_list.append([0.4,1])
        
        # select particle at origin and check for collision
        chosen_particle, chosen_tile, position = select_particle_and_tile(grid, True)
        overlap_count, overlapping_particles = check_for_collisions('UP', grid,
                                                chosen_particle, chosen_tile,
                                                4, 0.25)
        self.assertEqual(overlap_count, 1)
        self.assertNotEqual(overlapping_particles, [[-9.0,[-9.0,-9.0]]])

    def test_check_for_collisions_large_grid_PBC(self):
        '''verifies that the check_for_collisions function detects potentional 
        overlap and records them for use in the move function. it first checks
        collisions in its own tile and if none are found checks the neighbors
        in the direction traveled'''

        #building a four by four grid
        grid = [[None, None, None, None],
                [None, None, None, None],
                [None, None, None, None],
                [None, None, None, None]]

        build_grid(grid)
        populate_grid_neighbors(grid)

        # collision path
        grid[0][0].particle_list.append([0,0])
        grid[0][3].particle_list.append([0,3])

        # select particle at origin and check for collision
        chosen_particle, chosen_tile, position = select_particle_and_tile(grid, True, 3)
        overlap_count, overlapping_particles = check_for_collisions('UP', grid,
                                                chosen_particle, chosen_tile,
                                                4, 0.25)
        # one collision
        self.assertEqual(overlap_count, 1)
        self.assertNotEqual(overlapping_particles, [[-9.0,[-9.0,-9.0]]])


        # remove particle and put it in the tile above
        # and repeat the test
        grid[0][0].particle_list.pop()
        grid[0][1].particle_list.append([0,1])

        # select particle at origin and check for collision
        chosen_particle, chosen_tile, position = select_particle_and_tile(grid, True, 3)
        overlap_count, overlapping_particles = check_for_collisions('UP', grid,
                                                chosen_particle, chosen_tile,
                                                4, 0.25)
        # no collisions
        self.assertEqual(overlap_count, 0)
        self.assertEqual(overlapping_particles, [[-9.0,[-9.0,-9.0]]])

        # remove particle and put it one tile over and some distance up
        # and repeat the test
        grid[0][1].particle_list.pop()
        grid[0][0].particle_list.append([0,0.4])
        
        # select particle at origin and check for collision
        chosen_particle, chosen_tile, position = select_particle_and_tile(grid, True, 3)
        overlap_count, overlapping_particles = check_for_collisions('UP', grid,
                                                chosen_particle, chosen_tile,
                                                4, 0.25)
        self.assertEqual(overlap_count, 1)
        self.assertNotEqual(overlapping_particles, [[-9.0,[-9.0,-9.0]]])

    def test_check_for_collisions_large_grid_special_case(self):
        '''verifies that the check_for_collisions function detects potentional 
        overlap and records them for use in the move function. it first checks
        collisions in its own tile and if none are found checks the neighbors
        in the direction traveled'''

        #building a four by four dimentional grid
        grid = [[None, None, None, None],
                [None, None, None, None],
                [None, None, None, None],
                [None, None, None, None]]

        build_grid(grid)
        populate_grid_neighbors(grid)

        # collision path
        grid[0][0].particle_list.append([0,0])
        grid[0][1].particle_list.append([0,1])
        grid[0][2].particle_list.append([0,2])
        grid[1][0].particle_list.append([1,0])
        grid[1][1].particle_list.append([1,1])
        grid[1][2].particle_list.append([1,2])
        grid[1][3].particle_list.append([1,3])
        grid[2][0].particle_list.append([2,0])
        grid[2][1].particle_list.append([2,1])
        grid[2][2].particle_list.append([2,2])
        grid[3][0].particle_list.append([3,0])
        grid[3][1].particle_list.append([3,1])
        grid[3][2].particle_list.append([3,2])

#        print_particles_in_grid(grid)

        # select particle at origin and check for collision
        chosen_particle, chosen_tile, position = select_particle_and_tile(grid, True, 0)
        overlap_count, overlapping_particles = check_for_collisions('UP', grid,
                                                chosen_particle, chosen_tile,
                                                4, 0.25)

#        print(chosen_particle, chosen_tile, position)
        # one collision
        self.assertEqual(overlap_count, 3)
        self.assertNotEqual(overlapping_particles, [[-9.0,[-9.0,-9.0]]])

    def test_move(self):
        '''Verifies that the move function calls the collision function and 
        acts on the information provided. if a collision is found then it moves
        a particle until collision. if no collision then it moves one grid tile
        '''

        direction = "UP"
        rxy  = [[0.0, 2], [0, 3]]
        grid = [[None, None, None, None],
                [None, None, None, None],
                [None, None, None, None],
                [None, None, None, None]]

        build_grid(grid)
        populate_grid_neighbors(grid)
        assign_particles_to_grid(rxy, grid)

#        print_particles_in_grid(grid)

        chosen_particle, chosen_tile, position = select_particle_and_tile(grid, True, 2)
        overlap_count, overlapping_particles = check_for_collisions('UP', grid,
                                                chosen_particle, chosen_tile,
                                                4, 0.25)
#        print(chosen_particle, chosen_tile, position)

        move(direction, rxy, grid, chosen_particle, chosen_tile, position, 4)    
        

        



if __name__ == '__main__':
    unittest.main()
