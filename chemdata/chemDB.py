#!/usr/bin/env python
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

   print("Processing Gordon-McBride database...")
   process_mcbride(h5,'trans.inp')

   print("Processing other database...")
   process_database(h5,'database.species.db')

   #Most of these coefficients are accurate over a lower range
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



def process_burcat(h5,fname):
   f = open(fname,"r")
   lines = f.readlines()[88:] 

   lines = ''.join(lines)
 

   pat = r"[0-9]+\-[0-9]*\-[0-9].+?\n[ \t\n]*?\n"
   pat = re.compile(pat,re.DOTALL)

   matches = pat.findall(lines)

 
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
      for i in xrange(4):
         record[i] = record[i].expandtabs().rstrip()

         if len(record[i]) != 80 or record[i][79] != str(i+1): 
           print "|" + record[i] + "|"
           print "length of record: %d" %len(record[i])
           if (len(record) >= 80): print("Character: " + record[i][79])
           print("Error: record ID not %d"%(i+1))
           recorderr = 1

      if recorderr: continue


      # read the molecular weight
      try:
         d["MW"] = float(record[0][69:78])
      except:    
         print("="*66)
         print record[0]
         print record[1]
         print record[2]
         print record[3]
         print("="*66)
         print record[0][70:78]
         print record[0][60:80]
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
      for i in xrange(4):
         record[i] = record[i].expandtabs().rstrip()

         if len(record[i]) != 80 or record[i][79] != str(i+1): 
           print "|" + record[i] + "|"
           print "length of record: %d" %len(record[i])
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
