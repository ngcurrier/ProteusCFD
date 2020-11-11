#!/usr/bin/env python3
import os
import sys
import h5py
import numpy
import re
import math


# Explicit ranges for building the model
LOW_TEMP = 150
HI_TEMP = 6000
WARNINGS = True


R_UNIV = float(8.31446261815324) # J/K.mol

def getTemperatureRange(spdb):
    Tintervals = spdb.get('CEA_Tintervals')[()]
    lowT = spdb.get('CEA_T1')[()]
    hiT = spdb.get('CEA_T' + str(int(Tintervals)+1))[()]
    return [lowT, hiT]


# h - enthalpy (J/kg)
# s - entropy (J/kg.K)
# T_K - temperature (K)
def evaluateGibbsFreeEnergy(h, s, T_K):
    return h - T_K*s #(J/kg)

# h5 - hdf5 db object
# speciesList - list of species names (ie. ['H2O', 'H2'])
# massFractions - massfractions list in order of species, must sum to one
#                 NOTE: for best results pass largest massfraction species last
# T_K - temperature in Kelvin to evaluate
def evaluateMixtureS(h5, speciesList, massFractions, T_K):
        
    spdb = []
    for specie in speciesList:
        spdb.append(retrieveDBSpecies(h5, specie))

    slist = []
    for specie in spdb:
        slist.append(evaluateCEAS(specie, T_K))

    #compute mixture
    s = 0.0
    for imf, ies in zip(massFractions, slist):
        s += imf*ies

    return s # J/kg.K

# h5 - hdf5 db object
# speciesList - list of species names (ie. ['H2O', 'H2'])
# massFractions - massfractions list in order of species, must sum to one
#                 NOTE: for best results pass largest massfraction species last
# T_K - temperature in Kelvin to evaluate
def evaluateMixtureH(h5, speciesList, massFractions, T_K):
        
    spdb = []
    for specie in speciesList:
        spdb.append(retrieveDBSpecies(h5, specie))

    Hlist = []
    for specie in spdb:
        Hlist.append(evaluateCEAH(specie, T_K))

    #compute mixture
    h = 0.0
    for imf, ih in zip(massFractions, Hlist):
        h += imf*ih

    return h # J/kg


# evaluates mixture heat of formation at standard state
# h5 - hdf5 db object
# speciesList - list of species names (ie. ['H2O', 'H2'])
# massFractions - massfractions list in order of species, must sum to one
#                 NOTE: for best results pass largest massfraction species last
def evaluateMixtureHf(h5, speciesList, massFractions):
    spdb = []
    for specie in speciesList:
        spdb.append(retrieveDBSpecies(h5, specie))

    Hflist = []
    for specie in spdb:
        Hflist.append(evaluateCEAHf(specie))
        
    #compute mixture
    hf = 0.0
    for imf, ih in zip(massFractions, Hflist):
        hf += imf*ih

    return hf # J/kg

# evaluates mixture heat of formation at standard state
# h5 - hdf5 db object
# speciesList - list of species names (ie. ['H2O', 'H2'])
# massFractions - massfractions list in order of species, must sum to one
#                 NOTE: for best results pass largest massfraction species last
def evaluateMixtureDeltaHf(h5, speciesList, massFractions):
    spdb = []
    for specie in speciesList:
        spdb.append(retrieveDBSpecies(h5, specie))

    dHflist = []
    for specie in spdb:
        dHflist.append(evaluateCEADeltaHf(specie))
        
    #compute mixture
    dhf = 0.0
    for imf, ih in zip(massFractions, dHflist):
        dhf += imf*ih

    return dhf # J/kg


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
    MW = float(spdb.get('MW')[()])
    Rsp = R_UNIV/MW*1000.0

    Cp_R = a + b*T + c*T*T + d*T*T*T + e*T*T*T*T
    Cp = Cp_R*Rsp
    
    return Cp # J/kg.K
    
def evaluateNASA7H(spdb, T_K):
    T = T_K
    MW = float(spdb.get('MW')[()])
    Rsp = R_UNIV/MW*1000.0

    H_RT = a + b*T/2.0 + c*T*T/3.0 + d*T*T*T/4.0 + e*T*T*T*T/5.0 + f/T
    H = H_RT*T*Rsp

    return H # J/kg
    
def evaluateCEAS(spdb, T_K):
    # look for interval in which T_K lives
    Tintervals = spdb.get('CEA_Tintervals')[()]
    MW = float(spdb.get('MW')[()])
    Rsp = R_UNIV/MW*1000.0  # J/kg.K

    found = False
    for iT in range(0, Tintervals):
        T1key = 'CEA_T'+str(iT+1)
        T2key = 'CEA_T'+str(iT+2)
        T1 = spdb.get(T1key)[()]
        T2 = spdb.get(T2key)[()]
        if T_K > T1 and T_K <= T2:
            found = iT+1
    if not found:
        if WARNINGS:
            print('evaluate CEA H could not find temperature range for %f' % T_K)

    [lowT, hiT] = getTemperatureRange(spdb)
    if T_K < lowT:
        if WARNINGS:
            print('Bottom range is %f' % lowT)
            print('Extrapolating downward')
        found = 1
    elif T_K > hiT:
        if WARNINGS:
            print('Top range is %f' % hiT)
            print('Extrapolating upward')
        found = Tintervals

    exp = spdb.get('CEA_T'+str(found)+'_exp')[()]
    coeff = spdb.get('CEA7_T'+str(found))[()]

    # This is from the CEA user's manual
    # note, the NASA form has an a8*log(T) term here, CEA does not use
    S_R = -0.5*coeff[0]*T_K**exp[0] - \
            coeff[1]*T_K**exp[1] + \
            coeff[2]*T_K**exp[2]*math.log(T_K) + \
            coeff[3]*T_K**exp[3] + \
            coeff[4]*T_K**exp[4]*0.5 + \
            coeff[5]*T_K**exp[5]/3.0 + \
            coeff[6]*T_K**exp[6]*0.25 + \
            coeff[9] # coefficient b1 is indx 8 of 9, b2 is 9 of 9

    return S_R*Rsp # J/kg.K

# evaluates the enthalpy shift from 298.15 to 0K
def evaluateCEADeltaHf(spdb):
    # retrieve the heat of formation from CEA database
    dhf1 = spdb.get('CEA_hf298_T1_offset')[()] # J/mol
    dhf2 = spdb.get('CEA_hf298_T2_offset')[()] # J/mol
    # dhf1 and dhf2 are always the same values but we read them to 
    # have the ability to do a reader check here
    MW = float(spdb.get('MW')[()])

    # NOTE: that CEA gives this value in J/mol whereas the Burcat data
    #       gives us this value as hf(298.15K)/R
    dhf1_metric = dhf1/MW*1000.0 # J/kg
    dhf2_metric = dhf2/MW*1000.0 # J/kg
    return dhf1_metric #J/kg

# evaluates the heat of formation from CEA for a single species (j/kg)
def evaluateCEAHf(spdb):
    # retrieve the heat of formation from CEA database
    hf = spdb.get('CEA_hf')[()] # J/mol
    MW = float(spdb.get('MW')[()])

    # NOTE: that CEA gives this value in J/mol whereas the Burcat data
    #       gives us this value as hf(298.15K)/R
    hf_metric = hf/MW*1000.0 # J/kg
    # The evaluation of the curvefit and the hf(298.15K) should give the same value
    # or very close by definition

    return hf_metric #J/kg

# evaluates the heat of formation from CEA for a single species (j/kg)
def evaluateCEAH(spdb, T_K):
    # look for interval in which T_K lives
    Tintervals = spdb.get('CEA_Tintervals')[()]
    MW = float(spdb.get('MW')[()])
    Rsp = R_UNIV/MW*1000.0  # J/kg.K

    found = False
    for iT in range(0, Tintervals):
        T1key = 'CEA_T'+str(iT+1)
        T2key = 'CEA_T'+str(iT+2)
        T1 = spdb.get(T1key)[()]
        T2 = spdb.get(T2key)[()]
        if T_K > T1 and T_K <= T2:
            found = iT+1
    if not found:
        if WARNINGS:
            print('evaluate CEA H could not find temperature range for %f' % T_K)

    [lowT, hiT] = getTemperatureRange(spdb)
    if T_K < lowT:
        if WARNINGS:
            print('Bottom range is %f' % lowT)
            print('Extrapolating downward')
        found = 1
    elif T_K > hiT:
        if WARNINGS:
            print('Top range is %f' % hiT)
            print('Extrapolating upward')
        found = Tintervals

    exp = spdb.get('CEA_T'+str(found)+'_exp')[()]
    coeff = spdb.get('CEA7_T'+str(found))[()]

    # This is from the CEA user's manual
    H_RT = -coeff[0]*T_K**exp[0] + \
            coeff[1]*T_K**exp[1]*math.log(T_K) + \
            coeff[2]*T_K**exp[2] + \
            coeff[3]*T_K**exp[3]*0.5 + \
            coeff[4]*T_K**exp[4]/3.0 + \
            coeff[5]*T_K**exp[5]*0.25 + \
            coeff[6]*T_K**exp[6]*0.2 + \
            coeff[8]/T_K # coefficient b1 is indx 8 of 9, b2 is 9 of 9

    # NOTE: that CEA gives this value in J/mol whereas the Burcat data
    #       gives us this value as hf(298.15K)/R
    hf = spdb.get('CEA_hf')[()] # J/mol
    hf_metric = hf/MW*1000.0 # J/kg
    # The evaluation of the curvefit and the hf(298.15K) should give the same value
    # or very close by definition
       
    return H_RT*Rsp*T_K  #J/kg
    

def evaluateCEACp(spdb, T_K):
    # look for interval in which T_K lives
    Tintervals = spdb.get('CEA_Tintervals')[()]
    MW = float(spdb.get('MW')[()])
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
        if WARNINGS:
            print('evaluate CEA CP could not find temperature range for %f' % T_K)

    [lowT, hiT] = getTemperatureRange(spdb)
    if T_K < lowT:
        if WARNINGS:
            print('Bottom range is %f' % lowT)
            print('Extrapolating downward')
        found = 1
    elif T_K > hiT:
        if WARNINGS:
            print('Top range is %f' % hiT)
            print('Extrapolating upward')
        found = Tintervals

    exp = spdb.get('CEA_T'+str(found)+'_exp')[()]
    coeff = spdb.get('CEA7_T'+str(found))[()]

    # This is from the CEA user's manual
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

    if len(sys.argv) != 3:
        print('USAGE: %s <mixture CSV filename> <MW for non-dim (-1) for computed>' % sys.argv[0])
        exit()
    
    fileName = sys.argv[1]
    MW = float(sys.argv[2])

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
    lowT = []
    hiT = []
    #first line is header
    lines = f.readlines()
    for iline in range(1,len(lines)):
        line = lines[iline]
        # first value is species, second is massfraction
        # expect the file to be a csv
        tokens = line.replace(' ', '').split(',')
        speciesList.append(tokens[0])
        spdbs.append(retrieveDBSpecies(h5, speciesList[-1], verbose=True))
        print(spdbs[-1])
        massFractions.append(float(tokens[1]))
        MWs.append(getSpeciesMW(h5, speciesList[-1]))
        [lT, hT] = getTemperatureRange(spdbs[-1])
        lowT.append(lT)
        hiT.append(hT)

    # sanity check on mass fractions
    sum = 0.0
    for i in range(0,len(massFractions)-1):
        sum = sum + massFractions[i]

    # check for closeness of summation to unity
    if abs((1.0 - massFractions[-1]) - sum) > 1.0e-12:
        print('WARNING: last massfraction does not sum to one, it will be modified')

    massFractions[-1] = 1.0 - sum
    mixtureMW = 0.0
    print('\nSpecies\t\tMassFractions \tMW \tlowTemp \thiTemp')
    print('-------------------------------------------------------------')
    for isp, imf, imw, ilt, iht in zip(speciesList, massFractions, MWs, lowT, hiT):
        mixtureMW += imf/imw
        print('%s\t\t%.4f\t%.3f\t%.4f\t%.4f' % (isp+(' ')*6, imf, imw, ilt, iht))
    mixtureMW = 1.0/mixtureMW
        
    if MW < 0.0:
        MW = mixtureMW
    else:
        print('Using user input MW, not computed value!')

    listSpecies(h5)
    fout = open('resultsCP.csv', 'w')


    for spdb,isp in zip(spdbs,speciesList):
        print('\n*************** %s ********************' % isp)
        print('T(K)\t\t Cp(J/kg.K)')
        for i in range(0,120):
            T_K = 250 + i*50.0
            if T_K < LOW_TEMP:
                continue
            if T_K > float(HI_TEMP) :
                break
            cp = evaluateCEACp(spdb, T_K)
            print('%f\t %f' % (T_K, cp))
    
    Rmix = R_UNIV/MW*1000.0
    
    fout.write('T(K),Cv(J/kg.K),Cp(J/kg.K),Cp/R,H(J/kg),H(J/mol),h/RT,S(J/kg.K),s/R,G(J/kg),mu(Pa.s),k(W/m.K)\n')
    print('\n*************** MIXTURE ************************')
    print('T(K)\t\t Cv(J/kg.K) \tCp(J/kg.K)\tCp/R \t\t H(kJ/kg) \t H(kJ/mol) \t h/RT \t\t S(kJ/kg.K) \t s/R\t\t G(kJ/kg)\t mu(Pa.s)\t k(W/m.K)')
    print('----------------------------------------------------------------------------------------------------------------------------------------')
    for i in range(0,120):
        T_K = 250 + i*50.0
        if T_K < LOW_TEMP:
            continue
        if T_K > float(HI_TEMP):
            break
        cp = evaluateMixtureCp(h5, speciesList, massFractions, T_K)
        h = evaluateMixtureH(h5, speciesList, massFractions, T_K)
        s = evaluateMixtureS(h5, speciesList, massFractions, T_K)
        g = evaluateGibbsFreeEnergy(h, s, T_K)
        cv = cp - Rmix
        mu = 0.0 # todo
        k = 0.0  # todo
        print('%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f' % (T_K, cv, cp, cp/Rmix, h/1000.0, h/1000.0*MW/1000.0, h/(Rmix*T_K), s/1000.0, s/Rmix, g/1000.0, mu, k))
        fout.write('%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n' % (T_K, cv, cp, cp/Rmix, h, h*MW/1000.0, h/(Rmix*T_K), s, s/Rmix, g, mu, k))

    T_K = 298.15
    cp = evaluateMixtureCp(h5, speciesList, massFractions, T_K)
    h = evaluateMixtureH(h5, speciesList, massFractions, T_K)
    s = evaluateMixtureS(h5, speciesList, massFractions, T_K)
    g = evaluateGibbsFreeEnergy(h, s, T_K)
    cv = cp - Rmix
    mu = 0.0
    k = 0.0
    print()
    print('%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t%f' % (T_K, cv, cp, cp/Rmix, h/1000.0, h/1000.0*MW/1000.0, h/(Rmix*T_K), s/1000.0, s/Rmix, g/1000.0, mu, k))

    hf = evaluateMixtureHf(h5, speciesList, massFractions)
    dhf = evaluateMixtureDeltaHf(h5, speciesList, massFractions)

    print('\n*** Mixture MW computed (g/mol): %f ***' % mixtureMW)
    print('*** Mixture Rsp (J/kg.K): %f ***' % (R_UNIV/mixtureMW*1000.0))
    print('*** MW used in non-dimensionalization (g/mol): %f ***' % MW)
    print('*** Rsp used in non-dimensionalization (J/kg.K): %f ***' % Rmix)
    print('*** Hf(298.15K) computed (J/kg): %f ***' % hf)
    print('*** Hf(298.15K) - Hf(0K) computed (J/kg): %f ***' % dhf)
        
    fout.close()
    f.close()
    h5.close()
