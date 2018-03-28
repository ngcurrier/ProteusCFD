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
    found = False
    for line in fin:
        if line.find('{PressureKey}')  > 0:
            # replace the key with out value
            line2 = line.replace('{PressureKey}', str(pressure))
            fout.write(line2)
            found = True
        else:
            # echo the line
            fout.write(line)

    fin.close()
    fout.close()
            
    if found == False:
        raise("WARNING: did not find key {PressureKey} in *.inp file")

    p = subprocess.Popen('./FCEA2 < tmp.txt', stdout=subprocess.PIPE, shell=True)
    print str(p.communicate()[0])

def readMassFractions(filename, option):
    # Keep only the first N number of items
    maxNumberSpecies = 20

    #indxval = 1 # chamber
    #indxval = 2 # throat
    indxval = 3 # exit

    indxval = option

    fin = open(filename, 'r')

    i = 0
    lineFound = -999
    alllines = fin.readlines()
    for line in alllines:
        # find the list of output mole fractions
        if line.find('MOLE FRACTION') == 1:
            # save the line number and return
            lineFound = i
            print "Mole fraction key on line: " + str(i)
            break
        i = i + 1

    # echo the mole fraction block
    for id in range(lineFound+2, len(alllines)):
        if len(alllines[id].strip()) == 0:
            break
        print alllines[id]

    # jump to line in file, it is always followed by two
    # carriage returns, load dictionary with values
    molFrac = {}
    for id in range(lineFound+2, len(alllines)):
        if len(alllines[id].strip()) == 0:
            break
        else:
            tokens = alllines[id].split()
            # Tokens[0] contains element ID, tokens 1 contains mole fraction (ill-formatted)
            # Take constant width output from first column of mole fractions
            tokens[1] = alllines[id][17:25].replace(" " , "")  # chamber value
            tokens[2] = alllines[id][26:35].replace(" " , "")  # throat value
            tokens[3] = alllines[id][35:44].replace(" " , "")  # exit pressure 1 value
            # modify token since CEA sometimes outputs dumb stuff that's not really scientific notation
            # sometimes it skips the E- E+ notation
            fr = tokens[indxval].find('-')
            if(fr != -1 and fr != 0):
                modified = tokens[indxval][0:fr] + 'e' + tokens[indxval][fr:]
                print tokens[0] + ':\t' + tokens[indxval] + ' converted to readable ' + modified
                tokens[indxval] = modified
            molFrac[tokens[0]] = float(tokens[indxval])
    
    molFrac = sorted(molFrac.iteritems(), key=lambda(k,v):(v,k), reverse=True)
    Species = collections.OrderedDict()
    i = 0
    for frac in molFrac:
        Species[frac[0]] = frac[1]
        i = i + 1
        if i > maxNumberSpecies:
            break

    # ensure mol fractions add up to one
    # adjust last species to take up the remaineder
    sum = 0.0
    nspecies = len(Species)
    print 'Number of species read: ' + str(nspecies)
    i = 0
    for specie in Species:
        if i == nspecies -1:
            print 'Ajusting mole fraction to sum to 1.0 of ' + specie + ' ' +  str(Species[specie]) + ' to ' + str(1.0 - sum)
            if(1.0 - sum > 1.0e-3):
                print 'WARNING'
                print 'WARNING: correction balance > 1.0e-12, include more species in balance'
                print 'WARNING'
                exit()
            Species[specie] = 1.0 - sum
            if(Species[specie] < 0.0):
                print 'WARNING'
                print 'WARNING: mole fraction balance results in mole fraction < 0.0'
                print 'WARNING: clipping'
                Species[specie] = 0.0
        sum = sum + Species[specie]
        i = i + 1

    print "Sum of mole fractions: " + str(sum)
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
    beginP = 10
    if(len(sys.argv) != 3):
        print 'USAGE: ' + str(sys.argv[0]) + ' <CEA output filename> < 1-chamber, 2-throat, 3-exit>'
        exit()

    pressurePts = range(beginP,4000,20)
    pressures = collections.OrderedDict()
    for press in pressurePts:
        print "Pressure: " +  str(press)
        runCEA(sys.argv[1], press)
        species = readMassFractions('autoinput.out', int(sys.argv[2]))
        pressures[str(press)] = species


    # Since we are ordering the list, and the order may change
    # we need to use the keys directly to do the lookup
    specieOrder = pressures[str(beginP)].keys()
    print specieOrder

    # dump out the pressure dictionary
    print 'Pressure,',
    for specie in specieOrder:
        print specie + ',',
    print ''

    # output the table in CSV format
    for press in pressures:
        print press + ',',
        for specie in specieOrder:
            try:
                print str(pressures[press][specie]) + ',',
            except:
                # If species didn't make the list, it's not a driver, set it to zero
                print str(0.0) + ',',

        print ''
