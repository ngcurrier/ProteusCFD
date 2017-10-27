#!/usr/bin/python

import sys
import collections
import subprocess

def runCEA(inputFilename, pressure):
    fout = open('tmp.txt','w')
    fout.write('autoinput') # we will write autoinput.inp for CEA to run
    fout.close()

    fout = open('autoinput.inp','w')
    fin =  open(inputFilename,'r')
    for line in fin:
        if line.find('{PressureKey}')  > 0:
            # replace the key with out value
            line2 = line.replace('{PressureKey}', str(pressure))
            fout.write(line2)
        else:
            # echo the line
            fout.write(line)

    fin.close()
    fout.close()
            
    p = subprocess.Popen('./FCEA2 < tmp.txt', stdout=subprocess.PIPE, shell=True)
    print str(p.communicate()[0])

def readMassFractions(filename):
    maxNumberSpecies = 6
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
    for id in range(lineFound+2, len(alllines)):
        if len(alllines[id].strip()) == 0:
            break
        else:
            tokens = alllines[id].split()
            molFrac[tokens[0]] = float(tokens[1])
            
  
    molFrac = sorted(molFrac.items(), key=lambda x:x[1], reverse=True)
    # Keep only the first N number of items
    Species = collections.OrderedDict()
    for i in range(0,maxNumberSpecies):
        Species[molFrac[i][0]] = molFrac[i][1]

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
    # Dump some output for Proteus input files
    print 'moleFractions = [',
    i = 0
    for specie in Species:
        print str(Species[specie]),
        if i != len(Species) - 1:
            print ',',
        i = i + 1
    print ']\n'
    return Species
    
    
if __name__ == "__main__":
    if(len(sys.argv) != 2):
        print 'USAGE: ' + str(sys.argv[0]) + ' <CEA output filename>'
        exit()

    pressurePts = range(100,4000,100)
    pressures = collections.OrderedDict()
    for press in pressurePts:
        print "Pressure: " +  str(press)
        runCEA(sys.argv[1], press)
        species = readMassFractions('autoinput.out')
        pressures[str(press)] = species


    # dump out the pressure dictionary
    print 'Pressure,',
    for press in pressures:
        for specie in pressures[press]:
            print specie + ',',
        break
    print ''

    # output the table in CSV format
    for press in pressures:
        print press + ',',
        for specie in pressures[press]:
            print str(pressures[press][specie]) + ',',

        print ''
