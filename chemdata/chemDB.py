#!/usr/bin/env python3
import os
import sys
import h5py
import numpy
import re

def main():
   """
   This script digs into the Burcat and GRIMECH thermo databases, and the 
   McBride-Gordon database for transport coefficients, and dumps the data
   into a single database file (an HDF5 file) for use.  It also digs into
   our own text database for additional species data that is needed that
   happens not to be contained in the standard databases.
   """
  
   dbfile = 'chemdb.hdf5'

   os.system("rm %s"%dbfile);

   h5 = h5py.File(dbfile,'w')

   print("Processing Burcat database...")
   process_burcat(h5,'BURCAT_FIXED.THR')

   print("Processing CEA database...")
   process_CEA(h5, 'CEAthermo.inp')
   
   print("Processing Gordon-McBride database...")
   process_mcbride(h5,'trans.inp')

   print("Processing other database...")
   process_database(h5,'database.species.db')

   #Most of these coefficients are accurate over a lower 
   print("Processing GRIMECH thermo database...")
   process_grimech_thermo(h5,'GRIMECH.THR')
   
   h5.close()

   #Install to the correct location in home directory
   os.system("mkdir %s" % "~/.proteusCFD")
   os.system("mkdir %s" % "~/.proteusCFD/database")
   os.system("cp %s %s" % (dbfile, "~/.proteusCFD/database/"+dbfile))
   
   return


def process_database(h5,fname):
   f = open(fname,"r")
   lines = f.read() 
   
   # remove all of the comments
   pat = r"#.*"
   pat = re.compile(pat)
   lines = pat.sub('',lines)

   # find all of the groupings.
   pat = r"\{.*?\}"
   pat = re.compile(pat,re.DOTALL)
   matches = pat.findall(lines)

   # parse the groupings
   pat = r"[ ]*.*=.*"
   pat = re.compile(pat);
   for m in matches:
      d = {}
      keys = pat.findall(m);
      for k in keys:
         k = k.split('=')
         d[k[0].strip()] = k[1].strip()

      print("Database ==> %s"%d["compound"])

     
      # now we have a key, value dictionary with the data in there, but we need to
      # turn the value strings into floating point data.
      for k in d:
         v = d[k];
         if ',' in v:
            v = v.strip()
            if (v[-1] == ','): v = v[:-1]
            v = v.split(',')
            arr = numpy.zeros(len(v))
            for i,x in enumerate(v): arr[i] = float(x)
            d[k] = arr
         else:
            if k != "compound":
               d[k] = float(v)
         
      write_species(h5,d);

   return   


def process_mcbride(h5,fname):
   f = open(fname,"r")
   lines = f.readlines() 
   del lines[0]
 
   # this is a special list of ions were we assume that the transport properties are the same
   # as the neutrally charged compound.   We do this because there is no data for any ions in
   # the McBride-Gordon dataset.
   ions = {"O2":'+-',
            "O":'+-',
           "N2":'-+',
            "N":'+',
            "NO":'+', 
            "CO":'+',
            "OH":'+-',
            "H2":'+-',
            "H2O":'+'
           }

   while len(lines):
      d = {}
      header = lines[0]
      if header.strip() == "end": break
      compound = header.split()[0]
      d["compound"] = compound
      nv = int(header[35])
      nc = int(header[37])
      visc = lines[1:(1+nv)] 
      therm = lines[(1+nv):(1+nv+nc)] 
      del lines[0:(nv+nc+1)]

      if nv > 0:
        coeffs = numpy.zeros((nv,6))
        for i,l in enumerate(visc):
          T_lo = read_coeff(l,3,6)
          T_hi = read_coeff(l,11,7) 
          coeffs[i][0] = T_lo
          coeffs[i][1] = T_hi
          coeffs[i][2] = read_coeff(l,20)
          coeffs[i][3] = read_coeff(l,35)
          coeffs[i][4] = read_coeff(l,50)
          coeffs[i][5] = read_coeff(l,65)
        d["mu"] = coeffs

      if nc > 0:
        coeffs = numpy.zeros((nc,6))
        for i,l in enumerate(therm):
          T_lo = read_coeff(l,3,6)
          T_hi = read_coeff(l,11,7) 
          coeffs[i][0] = T_lo
          coeffs[i][1] = T_hi
          coeffs[i][2] = read_coeff(l,20)
          coeffs[i][3] = read_coeff(l,35)
          coeffs[i][4] = read_coeff(l,50)
          coeffs[i][5] = read_coeff(l,65)
        d["k"] = coeffs

      print("McBride/Gordon ==> %s"%d["compound"])
      write_species(h5,d)     

      # take care of the ions.
      if d["compound"] in ions:
         base = d["compound"]
         for i in ions[base]:
            d["compound"] = base+i
            print("McBride/Gordon ==> %s"%d["compound"])
            write_species(h5,d)     
        

      
   return

def process_CEA(h5,fname):
   f = open(fname,"r")
   
   # start reading after massive header
   lines = f.readlines()[42:] 
   nlines = len(lines)
   iline = 0
   
   # Records look like:
   #   Ar                Ref-Elm. Moore,1971. Gordon,1999.                             
   # 3 g 3/98 AR  1.00    0.00    0.00    0.00    0.00 0   39.9480000          0.000
   #    200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6197.428
   # 0.000000000D+00 0.000000000D+00 2.500000000D+00 0.000000000D+00 0.000000000D+00
   # 0.000000000D+00 0.000000000D+00                -7.453750000D+02 4.379674910D+00
   #   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6197.428
   # 2.010538475D+01-5.992661070D-02 2.500069401D+00-3.992141160D-08 1.205272140D-11
   #-1.819015576D-15 1.078576636D-19                -7.449939610D+02 4.379180110D+00
   #   6000.000  20000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6197.428
   #-9.951265080D+08 6.458887260D+05-1.675894697D+02 2.319933363D-02-1.721080911D-06
   # 6.531938460D-11-9.740147729D-16                -5.078300340D+06 1.465298484D+03


   while iline < nlines:
      d = {}
      # read species name 
      line = lines[iline]
      iline = iline + 1
      if line.find('END ') != -1:
         if iline >= nlines:
            break
         line = lines[iline]
         iline = iline + 1
      d["compound"] = line[0:17].strip()
      print("CEA database ==> %s" % d["compound"])
      # read #T intervals, MW, and heat of formation (298.15K)
      line = lines[iline]
      line=line.replace('D', 'E')
      iline = iline + 1
      Tintervals = int(line[0:2])
      d["CEA_Tintervals"] = Tintervals
      d["MW"] = float(line[52:64])
      d["hf"] = float(line[65:80])
      # read temperature ranges
      line = lines[iline]
      iline = iline + 1
      line=line.replace('D', 'E')
      items = line.split()
      d["CEA_T1"] = float(items[0])
      d["CEA_T2"] = float(items[1])
      d["hf298_T1_offset"] = float(items[-1])
      d["CEA_T1_exp"] = []
      # we have to slice out the correct region due to the way
      # the DB often writes numbers with no space between them
      expslice = line[24:58]
      items = expslice.split()
      d["CEA_T1_exp"].append(float(items[0]))
      d["CEA_T1_exp"].append(float(items[1]))
      d["CEA_T1_exp"].append(float(items[2]))
      d["CEA_T1_exp"].append(float(items[3]))
      d["CEA_T1_exp"].append(float(items[4]))
      d["CEA_T1_exp"].append(float(items[5]))
      d["CEA_T1_exp"].append(float(items[6]))

      if Tintervals > 0:
         # read coefficients
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         poly1 = []
         poly1.append(float(line[0:16]))
         poly1.append(float(line[16:32]))
         poly1.append(float(line[32:48]))
         poly1.append(float(line[48:64]))
         poly1.append(float(line[64:]))
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         poly1.append(float(line[0:16]))
         poly1.append(float(line[16:32]))
         if(line[32:48].strip() == ''):
            poly1.append(0.0)
         else:
            poly1.append(float(line[32:48]))
         # read integration constant b1
         poly1.append(float(line[48:64]))
         # read integration constant b2
         poly1.append(float(line[64:]))
         # Note: CEA uses a flexible form Cp curve fit where the coefficients are defined on the fly in the record
         # This is why we don't use the NASA7 key to identify them
         d["CEA7_T1"] = poly1
      # if the next (2nd) interval exists, read it
      if Tintervals > 1:
         # read temperature ranges
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         items = line.split()
         # todo: check T2 for sanity since it is listed twice
         d["CEA_T2"] = float(items[0])
         d["CEA_T3"] = float(items[1])
         d["hf298_T2_offset"] = float(items[-1])
         d["CEA_T2_exp"] = []
         # we have to slice out the correct region due to the way
         # the DB often writes numbers with no space between them
         expslice = line[24:58]
         items = expslice.split()
         d["CEA_T2_exp"] = []
         d["CEA_T2_exp"].append(float(items[0]))
         d["CEA_T2_exp"].append(float(items[1]))
         d["CEA_T2_exp"].append(float(items[2]))
         d["CEA_T2_exp"].append(float(items[3]))
         d["CEA_T2_exp"].append(float(items[4]))
         d["CEA_T2_exp"].append(float(items[5]))
         d["CEA_T2_exp"].append(float(items[6]))
         # read coefficients
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         poly1 = []
         poly1.append(float(line[0:16]))
         poly1.append(float(line[16:32]))
         poly1.append(float(line[32:48]))
         poly1.append(float(line[48:64]))
         poly1.append(float(line[64:]))
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         poly1.append(float(line[0:16]))
         poly1.append(float(line[16:32]))
         if(line[32:48].strip() == ''):
            poly1.append(0.0)
         else:
            poly1.append(float(line[32:48]))
         # read integration constant b1
         poly1.append(float(line[48:64]))
         # read integration constant b2
         poly1.append(float(line[64:]))
         # Note: CEA uses a flexible form Cp curve fit where the coefficients are defined on the fly in the record
         # This is why we don't use the NASA7 key to identify them
         d["CEA7_T2"] = poly1
      
      # if the next (3rd) interval exists, read it
      if Tintervals > 2:
         # read temperature ranges
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         items = line.split()
         # todo: check T3 for sanity since it is listed twice
         d["CEA_T3"] = float(items[0])
         d["CEA_T4"] = float(items[1])
         d["hf298_T3_offset"] = float(items[-1])
         d["CEA_T3_exp"] = []
         # we have to slice out the correct region due to the way
         # the DB often writes numbers with no space between them
         expslice = line[24:58]
         items = expslice.split()
         d["CEA_T3_exp"].append(float(items[0]))
         d["CEA_T3_exp"].append(float(items[1]))
         d["CEA_T3_exp"].append(float(items[2]))
         d["CEA_T3_exp"].append(float(items[3]))
         d["CEA_T3_exp"].append(float(items[4]))
         d["CEA_T3_exp"].append(float(items[5]))
         d["CEA_T3_exp"].append(float(items[6]))
         # read coefficients
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         poly1 = []
         poly1.append(float(line[0:16]))
         poly1.append(float(line[16:32]))
         poly1.append(float(line[32:48]))
         poly1.append(float(line[48:64]))
         poly1.append(float(line[64:]))
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         poly1.append(float(line[0:16]))
         poly1.append(float(line[16:32]))
         if(line[32:48].strip() == ''):
            poly1.append(0.0)
         else:
            poly1.append(float(line[32:48]))
         # read integration constant b1
         poly1.append(float(line[48:64]))
         # read integration constant b2
         poly1.append(float(line[64:]))
         # Note: CEA uses a flexible form Cp curve fit where the coefficients are defined on the fly in the record
         # This is why we don't use the NASA7 key to identify them
         d["CEA7_T3"] = poly1

      # if the next (4th) interval exists, read it
      if Tintervals > 3:
         # read temperature ranges
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         items = line.split()
         # todo: check T3 for sanity since it is listed twice
         d["CEA_T4"] = float(items[0])
         d["CEA_T5"] = float(items[1])
         d["hf298_T4_offset"] = float(items[-1])
         d["CEA_T4_exp"] = []
         # we have to slice out the correct region due to the way
         # the DB often writes numbers with no space between them
         expslice = line[24:58]
         items = expslice.split()
         d["CEA_T4_exp"].append(float(items[0]))
         d["CEA_T4_exp"].append(float(items[1]))
         d["CEA_T4_exp"].append(float(items[2]))
         d["CEA_T4_exp"].append(float(items[3]))
         d["CEA_T4_exp"].append(float(items[4]))
         d["CEA_T4_exp"].append(float(items[5]))
         d["CEA_T4_exp"].append(float(items[6]))
         # read coefficients
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         poly1 = []
         poly1.append(float(line[0:16]))
         poly1.append(float(line[16:32]))
         poly1.append(float(line[32:48]))
         poly1.append(float(line[48:64]))
         poly1.append(float(line[64:]))
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         poly1.append(float(line[0:16]))
         poly1.append(float(line[16:32]))
         if(line[32:48].strip() == ''):
            poly1.append(0.0)
         else:
            poly1.append(float(line[32:48]))
         # read integration constant b1
         poly1.append(float(line[48:64]))
         # read integration constant b2
         poly1.append(float(line[64:]))
         # Note: CEA uses a flexible form Cp curve fit where the coefficients are defined on the fly in the record
         # This is why we don't use the NASA7 key to identify them
         d["CEA7_T4"] = poly1

      # if the next (5th) interval exists, read it
      if Tintervals > 4:
         # read temperature ranges
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         items = line.split()
         # todo: check T3 for sanity since it is listed twice
         d["CEA_T5"] = float(items[0])
         d["CEA_T6"] = float(items[1])
         d["hf298_T5_offset"] = float(items[-1])
         d["CEA_T5_exp"] = []
         # we have to slice out the correct region due to the way
         # the DB often writes numbers with no space between them
         expslice = line[24:58]
         items = expslice.split()
         d["CEA_T5_exp"].append(float(items[0]))
         d["CEA_T5_exp"].append(float(items[1]))
         d["CEA_T5_exp"].append(float(items[2]))
         d["CEA_T5_exp"].append(float(items[3]))
         d["CEA_T5_exp"].append(float(items[4]))
         d["CEA_T5_exp"].append(float(items[5]))
         d["CEA_T5_exp"].append(float(items[6]))
         # read coefficients
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         poly1 = []
         poly1.append(float(line[0:16]))
         poly1.append(float(line[16:32]))
         poly1.append(float(line[32:48]))
         poly1.append(float(line[48:64]))
         poly1.append(float(line[64:]))
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         poly1.append(float(line[0:16]))
         poly1.append(float(line[16:32]))
         if(line[32:48].strip() == ''):
            poly1.append(0.0)
         else:
            poly1.append(float(line[32:48]))
         # read integration constant b1
         poly1.append(float(line[48:64]))
         # read integration constant b2
         poly1.append(float(line[64:]))
         # Note: CEA uses a flexible form Cp curve fit where the coefficients are defined on the fly in the record
         # This is why we don't use the NASA7 key to identify them
         d["CEA7_T5"] = poly1

      # if the next (6th) interval exists, read it
      if Tintervals > 5:
         # read temperature ranges
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         items = line.split()
         # todo: check T3 for sanity since it is listed twice
         d["CEA_T6"] = float(items[0])
         d["CEA_T7"] = float(items[1])
         d["hf298_T6_offset"] = float(items[-1])
         d["CEA_T6_exp"] = []
         # we have to slice out the correct region due to the way
         # the DB often writes numbers with no space between them
         expslice = line[24:58]
         items = expslice.split()
         d["CEA_T6_exp"].append(float(items[0]))
         d["CEA_T6_exp"].append(float(items[1]))
         d["CEA_T6_exp"].append(float(items[2]))
         d["CEA_T6_exp"].append(float(items[3]))
         d["CEA_T6_exp"].append(float(items[4]))
         d["CEA_T6_exp"].append(float(items[5]))
         d["CEA_T6_exp"].append(float(items[6]))
         # read coefficients
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         poly1 = []
         poly1.append(float(line[0:16]))
         poly1.append(float(line[16:32]))
         poly1.append(float(line[32:48]))
         poly1.append(float(line[48:64]))
         poly1.append(float(line[64:]))
         line = lines[iline]
         iline = iline + 1
         line=line.replace('D', 'E')
         poly1.append(float(line[0:16]))
         poly1.append(float(line[16:32]))
         if(line[32:48].strip() == ''):
            poly1.append(0.0)
         else:
            poly1.append(float(line[32:48]))
         # read integration constant b1
         poly1.append(float(line[48:64]))
         # read integration constant b2
         poly1.append(float(line[64:]))
         # Note: CEA uses a flexible form Cp curve fit where the coefficients are defined on the fly in the record
         # This is why we don't use the NASA7 key to identify them
         d["CEA7_T6"] = poly1

      write_species(h5,d)

   return


def process_burcat(h5,fname):
   f = open(fname,"r")
   
   # start reading after massive header
   lines = f.readlines()[88:] 

   lines = ''.join(lines)
 

   pat = r"[0-9]+\-[0-9]*\-[0-9].+?\n[ \t\n]*?\n"
   pat = re.compile(pat,re.DOTALL)

   matches = pat.findall(lines)

   fout = open('out.txt','w')
   for i in matches:
      fout.write(i+'\n')
   fout.close()
      
 
   for r in matches:
      r = r.strip().split('\n')
      record = r[-4:]

      if (len(record) != 4):
         print("Error: record not 4")
         continue

      # todo.......some of the burcat records are composite ones.  Right now, it just
      # takes the last one and ignores the rest.

      d = {}

      #read the compound name
      d["compound"] = record[0].split()[0] 

      print("Burcat database ==> %s"%d["compound"])

      recorderr = 0
      for i in range(4):
         record[i] = record[i].expandtabs().rstrip()

         if len(record[i]) != 80 or record[i][79] != str(i+1): 
           print("|" + record[i] + "|")
           print("length of record: %d" %len(record[i]))
           if (len(record) >= 80): print("Character: " + record[i][79])
           print("Error: record ID not %d"%(i+1))
           recorderr = 1

      if recorderr: continue


      # read the molecular weight
      try:
         d["MW"] = float(record[0][69:78])
      except:    
         print("="*66)
         print(record[0])
         print(record[1])
         print(record[2])
         print(record[3])
         print("="*66)
         print(record[0][70:78])
         print(record[0][60:80])
         raise

      # read the low and high temperatures of applicability for the curve fit.
      d["NASA7_burcat_T1"]  = float(record[0][47:55])
      d["NASA7_burcat_T2"]  = 1000.0
      d["NASA7_burcat_T3"]  = float(record[0][56:65])
    
      # read the high temperature polynomial coefficients (1000K-T2)
      polyhi = numpy.zeros(7) 
      polyhi[0] = read_coeff(record[1],0)
      polyhi[1] = read_coeff(record[1],15)
      polyhi[2] = read_coeff(record[1],30)
      polyhi[3] = read_coeff(record[1],45)
      polyhi[4] = read_coeff(record[1],60)
      polyhi[5] = read_coeff(record[2],0)
      polyhi[6] = read_coeff(record[2],15)

      # read the low temperature polynomial coefficients (T1-1000K)
      polylo = numpy.zeros(7) 
      polylo[0] = read_coeff(record[2],30)
      polylo[1] = read_coeff(record[2],45)
      polylo[2] = read_coeff(record[2],60)
      polylo[3] = read_coeff(record[3],0)
      polylo[4] = read_coeff(record[3],15)
      polylo[5] = read_coeff(record[3],30)
      polylo[6] = read_coeff(record[3],45)

      d["NASA7_burcat1"] = polylo
      d["NASA7_burcat2"] = polyhi

      # read the 15th coefficient
      d["NASA7_coeff15"] = read_coeff(record[3],60)


      write_species(h5,d)
       

   return




def process_grimech_thermo(h5,fname):
   """
   This routine processes the GRIMECH thermo file, replacing the burcat coefficients
   with those found in here.
   """
   f = open(fname,"r")
   lines = f.readlines()[5:] 

 
   while len(lines) >= 4:
      record = lines[:4]
      del lines[:4]

      if (len(record) != 4):
         print("Error: record not 4")
         continue

      d = {}

      #read the compound name
      d["compound"] = record[0].split()[0] 

      print("GRIMECH 3.0 database ==> %s"%d["compound"])

      recorderr = 0
      for i in range(4):
         record[i] = record[i].expandtabs().rstrip()

         if len(record[i]) != 80 or record[i][79] != str(i+1): 
           print("|" + record[i] + "|")
           print("length of record: %d" %len(record[i]))
           if (len(record) >= 80): print("Character: " + record[i][79])
           print("Error: record ID not %d"%(i+1))
           recorderr = 1

      if recorderr: continue



      # read the low and high temperatures of applicability for the curve fit.
      d["NASA7_grimech_T1"] = float(record[0][47:55])
      d["NASA7_grimech_T2"] = 1000.0 
      d["NASA7_grimech_T3"] = float(record[0][56:65])
      
      # read the high temperature polynomial coefficients (1000K-T2)
      polyhi = numpy.zeros(7) 
      polyhi[0] = read_coeff(record[1],0)
      polyhi[1] = read_coeff(record[1],15)
      polyhi[2] = read_coeff(record[1],30)
      polyhi[3] = read_coeff(record[1],45)
      polyhi[4] = read_coeff(record[1],60)
      polyhi[5] = read_coeff(record[2],0)
      polyhi[6] = read_coeff(record[2],15)

      # read the low temperature polynomial coefficients (T1-1000K)
      polylo = numpy.zeros(7) 
      polylo[0] = read_coeff(record[2],30)
      polylo[1] = read_coeff(record[2],45)
      polylo[2] = read_coeff(record[2],60)
      polylo[3] = read_coeff(record[3],0)
      polylo[4] = read_coeff(record[3],15)
      polylo[5] = read_coeff(record[3],30)
      polylo[6] = read_coeff(record[3],45)

      d["NASA7_grimech1"] = polylo
      d["NASA7_grimech2"] = polyhi

      # read the 15th coefficient
      #d["coeff15"] = read_coeff(record[3],60)


      write_species(h5,d,overwrite=True)
       

   return
      
      
         
def read_coeff(line,pos,nchar=15):
   """
   This fix is needed because there are some flaws in the burcat file; there is sometimes
   a space where the sign of the exponent should be.
   """
   coeff = line[pos:pos+nchar]

   if nchar == 15:
     if coeff[11:13] == "E ": coeff = coeff[0:12] + "+" + coeff[13:15] 
     if coeff[10:12] == " E": coeff = coeff[0:10] + "0" + coeff[11:15] 
   if "N/A" in coeff: return 0.0
   try:
     f = float(coeff)
   except: 
     print("line: "+line)
     print("attempted to convert: "+coeff)
     raise
   return f 
   

def write_species(h5,d,overwrite=False):
   
   compound = d["compound"]
   
   grp = h5.require_group('/species')
   grp = grp.require_group(compound)
   
   for k in d: 
      if k != "compound":
         if k not in grp: 
            grp.create_dataset(k,data=d[k])
         else:
            if not overwrite:
               print("Warning: skipping dataset "+compound+"/"+k)
            else:
               del grp[k]
               grp.create_dataset(k,data=d[k])
   return





if __name__ == "__main__":
   main()
