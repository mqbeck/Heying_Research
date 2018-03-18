#Name:                  Max Beck
#Instructor:            Dr. Heying
#Project:               Particles
#Term:                  Winter 2017

# This program aims to convert a particle simulator in fortran into
#`python for speed and language feature comparison
# python 3 is being used

from math import *
import random
import time

class Values:
    ncycles = 10000 # number of cycles where 'measurements are made'
    nequil_1 = 1000 # the first equilibration cycles
    nequil_2 = 1000 # the second equilibration cycles where the displacement vector is found
    n = 1500 # number of particles
    print("particles: ", n)
    sz = 1000 # bins holding particle distances
    rho = 0.5 # particle density
    q = 1.0 # particle ratio
    drmax = 0.5
    vol = n / rho
    delta = 0.0
    xA = 0.750 #mole fraction
    nswap = 20
    cavrad = 0.50
    grtest = 10
    factor = ncycles/grtest
    vAsum = None #necessary
    vBsum = None #necessary
    wsuma = None #necessary
    wsumb = None #necessary
    start = None
    end = None
    run_time = None
    
def main():
    start = time.time()
    
    # uses Particles and Coords classes
    # uses ran1, distrib, and initial functions
    nequil_1 = Values.nequil_1
    nequil_2 = Values.nequil_2
    ncycles = Values.ncycles
    print("nequil_1: ", nequil_1, "   nequil_2: ", nequil_2, "   ncycles: ", ncycles, "rho: ", Values.rho)
    nadjst = 100
    nswap = Values.nswap
    grtest = Values.grtest
    gmax = Values.sz
    delgr = 0.01
    gcaa = [None] * Values.sz # necessary
    gcab = [None] * Values.sz # necessary
    gcbb = [None] * Values.sz # necessary
    rxij = None #necessary
    rxnewi = None #necessary
    rxj = None # necessary
    ratio = None # necessary else declaration in a for loop
    va = None #necessary
    vb = None # necessary
    gbin = None # necesary
    plcd = None #necesary
    rid = None #necessary
    lid = None #necessary
    m = None
  
    rho = Values.rho
    cavrad = Values.cavrad

    q = Values.q

    delta = Values.delta
    Coords.nA = Values.xA * Values.n
    Coords.nB = Values.n - Coords.nA
    qhalf = q / 2.0
  
    drmax = Values.drmax
    accpt = 0
    box = Values.vol
    boxinv = 1.0 / box
    sigmaab = (q + 1.0) * (1.0 + delta) / 2.0
    
    print ("rho: ", rho)
    print ("into initial")#debugging
    initial(rho)
    print ("Out of initial")#debugging
    print ("nA & nB",Coords.nA,Coords.nB)
        

    #This first function checks the initial configuration
    #to alert you to any overlap
    overlap()
    
    #Initialize the random number generator - its default is to reference system time
    random.seed()
    
    print("Into main")
    ntot = nequil_1 + nequil_2 + ncycles
    
    #The big loop; each time through is a cycle
    for h in range(1, ntot):
        
        #A message to terminal to assure you that program is moving forward
        if h % 1000 == 0:
            print(h,"cycle")
        #endif
        
        #This loop repeats n + nswap times every cycle
        #On average, each particles gets moved and a nswaps are attempted
        for k in range(1, nswap + Values.n + 1):
            #This random number decideds where to attempt to move or swap a particle
            hh = random.randint(0,(nswap + Values.n - 1))
            #This random number decides which particle will be involved
            i = random.randint(0,(Values.n - 1))
            rxi = Coords.rx[i]
            pid = Coords.coordsid[i]
            
            #Particle movement attempt
            if hh < (Values.n + 1):
                #Move particle i a random distance in either direction
                rxnewi = rxi + random.uniform(-1,1) * drmax
                rxnewi = rxnewi - box * int(round(rxnewi * boxinv))
     
                #Check for overlap with particle i's new position
                noovlp = True
                for j in range(0, Values.n):
                    if j == i:
                        continue
                    #endif
                    rxij = rxnewi - Coords.rx[j]
                    rxij = rxij - box * int(round(rxij * boxinv))
                    rxij = abs(rxij)
                    if pid + Coords.coordsid[j] == 2:
                        if rxij < 1.0:
                            noovlp = False
                            break
                        #endif
                    elif pid + Coords.coordsid[j] == 3: 
                        if rxij < sigmaab:
                            noovlp = False
                            break
                        #endif
                    else:
                        if rxij < q:
                            noovlp = False
                            break
                        #endif
                    #endif
                #enddo

                #If there was no overlap, j should be n - 1 from loop above,
                #Unfortunately, it will still be that if it overlaps with particle n-1
                #if j == (Values.n - 1):
                if noovlp: 
                    Coords.rx[i] = rxnewi 
                    accpt = accpt + 1
                #endif
            else:
                #Attempt to swap particle i with particle j
                j = i
                #Ensure i&j are not the same particle
                while j==i:
                    j = random.randint(0,(Values.n - 1))
                #enddo
                
                rxj = Coords.rx[j]
                Coords.rx[i] = rxj
                Coords.rx[j] = rxi

                #Check for overlap with particle i's new position
                noovlp = True
                for m in range(0, Values.n):
                    if i == m:
                        continue
                    #endif
                    rxij = Coords.rx[i] - Coords.rx[m]
                    rxij = rxij - box * int(round(rxij * boxinv))
                    rxij = abs(rxij)
                    if pid + Coords.coordsid[m] == 2:
                        if rxij < 1.0:
                            noovlp = False
                            break
                        #endif
                    elif pid + Coords.coordsid[m] == 3:
                        if rxij < sigmaab:
                            noovlp = False
                            break
                        #endif
                    else:
                        if rxij < q:
                            noovlp = False
                            break
                        #endif
                    #endif
                #enddo
                #If there was no overlap, check for overlap with particle j's new position
                #if m == (Values.n - 1):
                if noovlp:
                    for m in range(0, Values.n):
                        if j == m:
                            continue
                        #endif
                        rxij = Coords.rx[j] - Coords.rx[m]
                        rxij = rxij - box * int(round(rxij * boxinv))
                        rxij = abs(rxij)
                        if Coords.coordsid[j] + Coords.coordsid[m] == 2:
                            if rxij < 1.0:
                                noovlp = False
                                break
                            #endif
                        elif Coords.coordsid[j] + Coords.coordsid[m] == 3:
                            if rxij < sigmaab:
                                noovlp = False
                                break
                            #endif
                        else:
                            if rxij < q:
                                noovlp = False
                                break
                            #endif
                        #endif
                    #enddo
                    #Swap was unsuccessful, due to overlap with j; move particles back
                    if noovlp==False:
                        Coords.rx[i] = rxi
                        Coords.rx[j] = rxj
                    #endif
                #Swap was unsuccessful, due to overlap with i; move particles back
                else:
                    Coords.rx[i] = rxi
                    Coords.rx[j] = rxj
                #endif
            #endif
        #enddo

        #end of particle movement
        #Time for adjustments and sampling of configuration
        if h >= nequil_1:
            if  h % nadjst == 0:
                ratio = accpt/(Values.n*nadjst)
                if ratio > 0.4:
                    drmax = drmax * 1.05
                else:
                    drmax = drmax * 0.95
                #endif
                accpt = 0
                if drmax > (box/2.0):
                    drmax = box/2.0
                #endif
            #endif
        #endif

        #Reset all the bins and counters
        if h == (nequil_1 + nequil_2):
            #for i in range(0, gmax + 1):
            giaa = [0] * gmax
            giab = [0] * gmax
            gibb = [0] * gmax
            gcaa = [0] * gmax
            gcab = [0] * gmax
            gcbb = [0] * gmax
        
            Values.vAsum = 0.0
            Values.vBsum = 0.0
            Values.wsuma = 0.0
            Values.wsumb = 0.0
        #endif

        if h > (nequil_1 + nequil_2) and (h % grtest) == 0:
            #To find particle-particle distances, you only need to go to second to last
            #since the next loop with find the distance between it and the last
            for i in range(0, Values.n - 1):
                pid = Coords.coordsid[i]
                rxi = Coords.rx[i]
                
                for j in range(i + 1, Values.n):
                    rxij = rxi - Coords.rx[j]
                    rxij = rxij - box * int(round(rxij * boxinv))
                    rxij = abs(rxij)
                    gbin = int(rxij/delgr) 
                    if gbin >= gmax:
                        gbin = gmax - 1
                    #endif
                    if pid + Coords.coordsid[j] == 2:
                        giaa[gbin] = giaa[gbin] + 1
                    elif pid + Coords.coordsid[j] == 3:
                        giab[gbin] = giab[gbin] + 1
                    else:
                        gibb[gbin] = gibb[gbin] + 1
                    #endif
                #enddo
            #enddo
    # the section below is used to determine the
    # excess chemical potential of the fluid

            plcd = 0
            while True:
                #This routine is based on knowing the number of successful cavity insertions
                #so first, a cavity has to be inserted such that it contains no centers
                
                rxi = random.uniform(-1,1) * box
                rmin = 10.0
                rid = 0
                lmin = 10.0
                lid = 0
                #Check to see whether cavity includes any centers
                #Also, find the closest particle on the left and right for later
                noovlp = True
                for j in range(0, Values.n):
                    rxij = rxi - Coords.rx[j]
                    rxij = rxij - box * int(round(rxij * boxinv))
                    if rxij > 0.0:
                        if rxij < rmin:
                            rmin = rxij
                            rid = Coords.coordsid[j]
                        #endif
                    else:
                        rxij = abs(rxij)
                        if rxij < lmin:
                            lmin = rxij
                            lid = Coords.coordsid[j]
                        #endif
                    #endif
                    if rxij < cavrad:
                        noovlp = False
                        break
                    #endif
                #enddo
                #If cavity location is OK, find available space within cavity
                #if j == (Values.n - 1):
                if noovlp:
                    if rid == 1:
                        if rmin > 1.5:
                            rvol = 0.5
                        else:
                            rvol = rmin - 1.0
                        #endif
                    else:
                        if rmin > (qhalf + 1.0):
                            rvol = 0.5
                        else:
                            rvol = rmin - (qhalf + 0.5)
                        #endif
                    #endif
                    if lid == 1:
                        if lmin > 1.5:
                            lvol = 0.5
                        else:
                            lvol = lmin - 1.0
                        #endif
                    else:
                        if lmin > (qhalf + 1.0):
                            lvol = 0.5
                        else:
                            lvol = lmin - (qhalf + 0.5)
                        #endif
                    #endif
                    #Available space in cavity for another A particle
                    ssgesum = rvol + lvol
                    if ssgesum > 0.0:
                        vA = ssgesum
                    else:
                        vA = 0.0
                    #endif
                    Values.vAsum = Values.vAsum + vA

                    if rid == 1:
                        if rmin > (qhalf + 1.0):
                            rvol = 0.5
                        else:
                            rvol = rmin - (qhalf + 0.5)
                        #endif
                    else:
                        if rmin > (q + 0.5):
                            rvol = 0.5
                        else:
                            rvol = rmin - q
                        #endif
                    #endif
                    if lid == 1:
                        if lmin > (qhalf + 1.0):
                            lvol = 0.5
                        else:
                            lvol = lmin - (qhalf + 0.5)
                        #endif
                    else:
                        if lmin > (q + 0.5):
                            lvol = 0.5
                        else:
                            lvol = lmin - q
                        #endif
                    #endif
                    #Available space in cavity for another B particle
                    ssgesum = rvol + lvol
                    if ssgesum > 0.0:
                        vB = ssgesum
                    else:
                        vB = 0.0
                    #endif
                    Values.vBsum = Values.vBsum + vB
                    plcd = plcd + 1
                    
                #endif
                if plcd == 100:
                    break
                #endif
            #enddo
      
            #Widom insertion technique, 100 insertions, resulting in a 0 or 1
            for i in range(0, 100):
                rxi = random.uniform(-1,1) * box
                aok = True
                bok = True
                
                for j in range(0, Values.n):
                    rxij = rxi - Coords.rx[j]
                    rxij = rxij - box * int(round(rxij * boxinv))
                    rxij = abs(rxij)

                    #Potential for A
                    if Coords.coordsid[j] == 1:
                        if rxij < 1.0:
                            aok = False
                        #endif
                    else:
                        if rxij < sigmaab:
                            aok = False
                        #endif
                    #endif
                    
                    #Potential for B
                    if Coords.coordsid[j] == 1:
                        if rxij < sigmaab:
                            bok = False
                        #endif
                    else:
                        if rxij < q:
                            bok = False
                        #endif
                    #endif
                    if aok == False:
                        break
                    #endif
                #enddo
                if aok:
                    Values.wsuma = Values.wsuma + 1.0
                #endif
                if bok:
                    Values.wsumb = Values.wsumb + 1.0
                #endif
            #enddo
        #endif
    #enddo

# finish some sort of timing here
    end = time.time()
    Values.run_time = end - start
# wrtie to file instead of printing

    distrib(gmax, rho, Values.factor, delgr, box, giaa, giab, gibb, gcaa, gcab, gcbb)
#    print("ravgi \t", "gaa(r)i \t", "gab(r)i \t ", "gbb(r)")
#    for i in range(gmax - 1):
#        ravg = delgr * (2.0 * i + 1.0) / 2.0
#        print('ravg: {0:7.4f}    graa: {1:7.4f}    grab: {2:7.4f}   grbb: {3:7.4f}'.format(
#                ravg, Distrib.graa[i], Distrib.grab[i], Distrib.grbb[i]))
        
    output_file = "output_{0}_{1}_{2}".format(ncycles, 
    nequil_1 + nequil_2, Values.n)
    
    write_to_file(output_file, gmax, delgr)
    
    #These need to be converted to chemical potentials, but that involves ln; if it is not
    #working ln 0 or infinity is a common output but not as helpful in terms of debugging
#    print("drmax: ",drmax)
    
#    print("Available A space: ",Values.vAsum/Values.factor/100)
#    print("Available B space: ",Values.vBsum/Values.factor/100)
    
#    print("Widom A space: ",Values.wsuma/Values.factor/100)
#    print("Widom B space: ",Values.wsumb/Values.factor/100)
    
       # print('ravg: {0:7.4f}    giaa: {1:7.4f}    giab: {2:7.4f}   gibb: {3:7.4f}'.format(
        #        ravg, giaa[i], giab[i], gibb[i]))
        
        # print('ravg: {0:7.4f}    gcaa: {1:7.4f}    gcab: {2:7.4f}   gcbb: {3:7.4f}'.format(
        #        ravg, gcaa[i], gcab[i], gcbb[i]))

# check the scope of variables in theses classes
# usable anywhere with value retention as Class.variable


    
def overlap():

#This first double loop checks the initial configuration
#to alert you to any overlap

    for i in range(0, Values.n - 1):
        pid = Coords.coordsid[i]
        rxi = Coords.rx[i]
        for j in range(i + 1, Values.n):
            rxij = rxi - Coords.rx[j]
            rxij = rxij - Values.vol * int(round(rxij * (1 / Values.vol)))
            rxij = abs(rxij)
       
            if pid + Coords.coordsid[j] == 2:
                if rxij < 1.0:
                    print ("Overlap!", i, j, rxij)
                #endif
            elif pid + Coords.coordsid[j] == 3:
                sigmaab = (Values.q + 1.0) * (1.0 + Values.delta) / 2.0
                if rxij < sigmaab:
                    print ("Overlap!", i, j, rxij)
                #endif
            else:
                if rxij < Values.q:
                    print ("Overlap!", i, j, rxij)
                #endif
            #endif
        #enddo
    #enddo
    
def write_to_file(file_name, gmax, delgr):
    output_file = open(file_name, 'w')
    
    output_file.write(
"""System parameters
Number of particles:        {0:>7.4f}
Reduced density:            {1:>7.4f}
Diameter Ratio:             {2:>7.4f}
Mole Fraction:              {3:>7.4f}
Packing fraction:           {4:>7.4f}
Non-additivity:             {5:>7.4f}

Equilibration Cycles:       {6:>7.4f}
Number of Cycles:           {7:>7.4f}
Number of MC moves:         {8:>7.4f}
Attempted swaps per cycles: {9:>7.4f}
Drmax:                      {10:>7.4f}
 
System time (minutes):      {11:>7.4f}    

Cavity Radius:              {12:>7.4f}
Available A space:          {13:>7.4f}   
Available B space:          {14:>7.4f}

Widom A space:              {15:>7.4f}
Widom B space:              {16:>7.4f}\n""".format(
            Values.n, 
            Values.rho,
            Values.q,
            Values.xA, 
            Values.rho*(Values.xA*(1.0  - Values.q) + Values.q),
            Values.delta,
            Values.nequil_1 + Values.nequil_2, 
            Values.ncycles, 
            Values.n * Values.ncycles, 
            Values.nswap, 
            Values.drmax,
            Values.run_time / 60, 
            Values.cavrad, 
            Values.vAsum/Values.factor/100.0,
            Values.vBsum/Values.factor/100.0,
            Values.wsuma/Values.factor/100.0,
            Values.wsumb/Values.factor/100.0))
  
    for i in range(gmax - 1):
        ravg = delgr * (2.0 * i + 1.0) / 2.0
        output_file.write(
            'ravg: {0:7.4f}    graa: {1:7.4f}    grab: {2:7.4f}      grbb: {3:7.4f}\n'.format(ravg,Distrib.graa[i], 
                         Distrib.grab[i], Distrib.grbb[i]))
class Particles:
    n = Values.n
    print("particles: ", n)
    sz = 100

class Coords:
    rx = [None] * Values.n   # initalize an empty list of length n
    q = None
    coordsid = [None] * Values.n
    nA = None
    nB = None

class Distrib:
    gdaa = None
    gdab = None
    gdbb = None
    
    graa = None
    grab = None
    grbb = None

    ravg = None

# check how scope affects these variables
def initial(param_rho):
    rho = param_rho
    rxref = 0.0
    spc = 1.0 / rho
            
    for i in range(0, Values.n):
        Coords.rx[i] = rxref + spc * (i)
        if i > Coords.nA:
            Coords.coordsid[i] = 2
        else:
            Coords.coordsid[i] = 1
        #endif
    #enddo

def distrib(gmax, rho, factor, delgr, box, gdaa, 
            gdab, gdbb, graa, grab, grbb):
    
    graa = [None] * Values.sz
    grab = [None] * Values.sz
    grbb = [None] * Values.sz

    volshell = None
    
    for i in range(0, gmax - 1):
        rlower = (i + 1) * delgr
        rupper = rlower + delgr
        volshell = (rupper - rlower)
        Distrib.ravg = gdaa[i] / (factor * (Coords.nA * Coords.nA))   
        graa[i] = Distrib.ravg * box / volshell
        Distrib.ravg = gdab[i] / (factor * (2*Coords.nA * Coords.nB))  
        grab[i] = Distrib.ravg * box / volshell
        Distrib.ravg = gdbb[i] / (factor * (Coords.nB * Coords.nB))   
        grbb[i] = Distrib.ravg * box / volshell


    Distrib.graa = graa
    Distrib.grab = grab
    Distrib.grbb = grbb
    
    
if __name__ == "__main__":
    main()    
