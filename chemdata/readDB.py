#!/usr/bin/env python3
import os
import sys
import h5py
import numpy
import re
import math

R_UNIV = 8.31446261815324 # J/K.mol

def evaluateCEACp(spdb, T_K):
    # look for interval in which T_K lives
    Tintervals = spdb.get('CEA_Tintervals')[()]
    MW = spdb.get('MW')[()]

    found = False
    for iT in range(0, Tintervals):
        T1key = 'CEA_T'+str(iT+1)
        T2key = 'CEA_T'+str(iT+2)
        T1 = spdb.get(T1key)[()]
        T2 = spdb.get(T2key)[()]
        if T_K > T1 and T_K < T2:
            found = iT+1
    if not found:
        raise ValueError('evaluate CEA CP could not find temperature range for %f' % T_K)
    exp = spdb.get('CEA_T'+str(found)+'_exp')[()]
    coeff = spdb.get('CEA7_T'+str(found))[()]

    Cp_R = coeff[0]*T_K**exp[0] + coeff[1]*T_K**exp[1] + coeff[2]*T_K**exp[2] \
           + coeff[3]*T_K**exp[3] + coeff[4]*T_K**exp[4] + coeff[5]*T_K**exp[5] \
           + coeff[6]*T_K**exp[6]

    Rsp = R_UNIV/MW*1000.0
    return Cp_R*Rsp

def retrieveDBSpecies(h5, species):
    # strip off asterisk if they exist (output common from CEA)
    species = species.strip('*')
    
    #species = species.upper()
    print('\n-----------------------------------------------------')
    print('Retrieving species: %s' % species)
    print('-----------------------------------------------------')
    spList = h5.get('species')
    speciesObject = spList.get(species)
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
        print(item)
    f.close()
        
if __name__ == "__main__":

    speciesName = sys.argv[1]
    
    # this is the standard location after installation
    dbfile = os.path.expanduser('~/.proteusCFD/database/chemdb.hdf5')
    print('Reading: %s' % dbfile)
    
    h5 = h5py.File(dbfile, 'r')

    listSpecies(h5)
    
    spdb = retrieveDBSpecies(h5, speciesName)


    print('T(K)\t Cp(J/Kg.K)')
    for i in range(0,20):
        T_K = 287 + i*50.0
        cp = evaluateCEACp(spdb, T_K)
        print('%f\t %f' % (T_K, cp))
    
    h5.close()
