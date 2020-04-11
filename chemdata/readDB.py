#!/usr/bin/env python3
import os
import sys
import h5py
import numpy
import re
import math

R_UNIV = float(8.31446261815324) # J/K.mol

# h5 - hdf5 db object
# speciesList - list of species names (ie. ['H2O', 'H2'])
# massFractions - massfractions list in order of species, must sum to one
#                 NOTE: for best results pass largest massfraction species last
# T_K - temperature in Kelvin to evaluate
def evaluateMixtureCp(h5, speciesList, massFractions, T_K):
        
    spdb = []
    for specie in speciesList:
        spdb.append(retrieveDBSpecies(h5, specie))

    Cplist = []
    for specie in spdb:
        Cplist.append(evaluateCEACp(specie, T_K))

    #compute mixture
    cp = 0.0
    for imf, icp in zip(massFractions, Cplist):
        cp += imf*icp

    return cp # J/kg.K

def evaluateNASA7Cp(spdb, T_K):
    T = T_K
    MW = spdb.get('MW')[()]
    Rsp = R_UNIV/MW*1000.0

    Cp_R = a + b*T + c*T*T + d*T*T*T + e*T*T*T*T
    Cp = Cp_R*Rsp
    
    return Cp # J/kg.K
    
def evaluateNASA7H(spdb, T_K):
    T = T_K
    MW = spdb.get('MW')[()]
    Rsp = R_UNIV/MW*1000.0

    H_RT = a + b*T/2.0 + c*T*T/3.0 + d*T*T*T/4.0 + e*T*T*T*T/5.0 + f/T
    H = H_RT*T*Rsp

    return H # J/kg
    
def evaluateCEACp(spdb, T_K):
    # look for interval in which T_K lives
    Tintervals = spdb.get('CEA_Tintervals')[()]
    MW = spdb.get('MW')[()]
    Rsp = R_UNIV/MW*1000.0

    found = False
    for iT in range(0, Tintervals):
        T1key = 'CEA_T'+str(iT+1)
        T2key = 'CEA_T'+str(iT+2)
        T1 = spdb.get(T1key)[()]
        T2 = spdb.get(T2key)[()]
        if T_K > T1 and T_K <= T2:
            found = iT+1
    if not found:
        raise ValueError('evaluate CEA CP could not find temperature range for %f' % T_K)
    exp = spdb.get('CEA_T'+str(found)+'_exp')[()]
    coeff = spdb.get('CEA7_T'+str(found))[()]

    Cp_R = coeff[0]*T_K**exp[0] + coeff[1]*T_K**exp[1] + coeff[2]*T_K**exp[2] \
           + coeff[3]*T_K**exp[3] + coeff[4]*T_K**exp[4] + coeff[5]*T_K**exp[5] \
           + coeff[6]*T_K**exp[6]

    return Cp_R*Rsp # J/kg.K

def retrieveDBSpecies(h5, species, verbose=False):
    # strip off asterisk if they exist (output common from CEA)
    species = species.strip('*')
    
    spList = h5.get('species')
    speciesObject = spList.get(species)

    #species = species.upper()
    if verbose:
        print('\n-----------------------------------------------------')
        print('Retrieving species: %s' % species)
        print('-----------------------------------------------------')
        avail = (list(speciesObject.items()))
        
        # loops over all keys like 'MW', etc
        for item in avail:
            key = item[0]
            array = speciesObject.get(key)
            array = array[()]
            print()
            print(key + ': ' + str(array))

    return speciesObject
        
def listSpecies(h5):
    spList = h5.get('species')
    spList = sorted(spList)
    f = open('speciesList.txt','w')
    for item in spList:
        f.write(item+'\n')
    f.close()

def getSpeciesMW(h5, species):
    sp = retrieveDBSpecies(h5, species)
    MW = float(sp.get('MW')[()])
    return MW # g/mol
    
if __name__ == "__main__":

    if len(sys.argv) != 2:
        print('USAGE: %s <mixture CSV filename>' % sys.argv[0])
        exit()
    
    fileName = sys.argv[1]

    f = open(fileName, 'r')
    print('Reading mixture file: %s' % fileName)

    # this is the standard location after installation
    dbfile = os.path.expanduser('~/.proteusCFD/database/chemdb.hdf5')
    print('Reading: %s' % dbfile)
    
    h5 = h5py.File(dbfile, 'r')

    speciesList = []
    spdbs = []
    massFractions = []
    MWs = []
    #first line is header
    lines = f.readlines()
    for iline in range(1,len(lines)):
        line = lines[iline]
        # first value is species, second is massfraction
        # expect the file to be a csv
        tokens = line.replace(' ', '').split(',')
        speciesList.append(tokens[0])
        spdbs.append(retrieveDBSpecies(h5, speciesList[-1], verbose=True))
        massFractions.append(float(tokens[1]))
        MWs.append(getSpeciesMW(h5, speciesList[-1]))
  
    # sanity check on mass fractions
    sum = 0.0
    for i in range(0,len(massFractions)-1):
        sum = sum + massFractions[i]

    # check for closeness of summation to unity
    if abs((1.0 - massFractions[-1]) - sum) > 1.0e-12:
        print('WARNING: last massfraction does not sum to one, it will be modified')

    massFractions[-1] = 1.0 - sum
    mixtureMW = 0.0
    print('\nSpecies\t\tMassFractions \tMW')
    print('--------------------------------------------')
    for isp, imf, imw in zip(speciesList, massFractions, MWs):
        mixtureMW += imf/imw
        print('%s\t\t%.4f\t\t%.4f' % (isp, imf, imw))
    mixtureMW = 1.0/mixtureMW
    print('\n*** Mixture MW (g/mol): %f ***' % mixtureMW)
    print('*** Mixture Rsp (J/kg.K): %f ***' % (R_UNIV/mixtureMW*1000.0))
        

    listSpecies(h5)
    fout = open('resultsCP.csv', 'w')


    for spdb,isp in zip(spdbs,speciesList):
        print('\n*************** %s ********************' % isp)
        print('T(K)\t\t Cp(J/kg.K)')
        for i in range(0,80):
            T_K = 250 + i*50.0
            cp = evaluateCEACp(spdb, T_K)
            print('%f\t %f' % (T_K, cp))
    
    
    fout.write('T(K),Cp(J/kg.K)\n')
    print('\n*************** MIXTURE ************************')
    print('T(K)\t\t Cp(J/kg.K)')
    for i in range(0,80):
        T_K = 250 + i*50.0
        cp = evaluateMixtureCp(h5, speciesList, massFractions, T_K)
        print('%f\t %f' % (T_K, cp))
        fout.write('%f,%f\n' %(T_K, cp))
        
    fout.close()
    f.close()
    h5.close()
