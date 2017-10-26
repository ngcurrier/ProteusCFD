#!/usr/bin/python

import sys
import collections

def main(filename):
    fin = open(filename, 'r')

    i = 0
    lineFound = -999
    alllines = fin.readlines()
    for line in alllines:
        # find the list of output mole fractions
        if line.find('MOLE FRACTION') == 1:
            # save the line number and return
            lineFound = i
            break
        i = i + 1

    # jump to line in file, it is always followed by two
    # carriage returns, load dictionary with values
    molFrac = {}
    print alllines[lineFound]
    for id in range(lineFound+2, len(alllines)):
        if len(alllines[id].strip()) == 0:
            break
        else:
            tokens = alllines[id].split()
            molFrac[tokens[0]] = float(tokens[1])
            
  
    molFrac = sorted(molFrac.items(), key=lambda x:x[1], reverse=True)
    # Keep only the first N number of items
    Species = {}
    N = 6
    for i in range(0,N):
        Species[molFrac[i][0]] = molFrac[i][1]

    print Species

    # ensure mol fractions add up to one
    # adjust last species to take up the remaineder
    sum = 0.0
    nspecies = len(Species)
    i = 0
    for specie in Species:
        if i == nspecies -1:
            Species[specie] = 1.0 - sum
        sum = sum + Species[specie]
        i = i + 1

    print "Sum of mole fractions: " + str(sum)
    print Species
        
if __name__ == "__main__":
    print 'USAGE: ' + str(sys.argv[0]) + ' <CEA output filename>'
    main(sys.argv[1])
