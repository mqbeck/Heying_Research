from math import *
import random
import time
import copy

class Values:
    diameter = 1
    ncycles = 2000
    nparticles = 250
    rho = 0.5
    volume_length = nparticles / rho
    inv_length = 1 / volume_length
    length = nparticles ** 0.5
    rx = [None] * nparticles
    rx_init = None
    distances = [] #particle-particle distances to be sorted and counted, per cycle
    histogram = [] #sorted distance, ultimately normalized (summed to one)
    cycles_to_average = 2
    time

def main():
    #Values.rx = define_list(Values.nparticles,Values.volume_length)
    random.seed(405)
    start = time.time()
    define_list(Values.nparticles,Values.volume_length,Values.diameter,Values.rx)
   # sort_list(Values.rx)
    Values.rx_init = copy.deepcopy(Values.rx)
    #print(Values.rx)
    single_event_chain()
    sort_list(Values.rx)
    #print(Values.rx)
    end = time.time()
    Values.time = (end - start)
    #build_histogram(Values.rx)
    filename = "{0}_particles_{1}_cycles_{2}_minutes.txt".format(Values.nparticles, Values.ncycles, Values.time/60)
    write_to_file(filename, Values.rx, Values.rx_init)
    print(Values.rx)
    
def define_list(nparticles,volume_length,diameter,rx):
    #return [(random.randint(0,(Values.volume_length))) for i in range(nparticles)]
    vector = volume_length/nparticles - diameter
    xpos = 0
    for i in range(1,nparticles + 1):
        xpos += random.random()*vector
        rx[i-1]=xpos
        xpos += diameter
    #enddo    
    
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