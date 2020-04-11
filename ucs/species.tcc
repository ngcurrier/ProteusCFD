#include "exceptions.h"
#include <iomanip>

template <class Type>
Species<Type>::Species()
{
  hasViscousProps = false;
  k_coeff_curves = 0;
  mu_coeff_curves = 0;
  charge = 0;

  //set these to air standard, this will be used when other curvefits don't extend low enough
  //From White Viscous Fluid Flow p.29 and p.32
  k_coeff_White[0] = 0.0241;  //ref_k  - W/m.K
  k_coeff_White[1] = 273;  //ref_temp
  k_coeff_White[2] = 194;  //S
  k_transition_White = 1000;   //this number is simply picked to cover lower range

  mu_coeff_White[0] = 1.716e-5;  //ref_mu N.s/m^2 (Pa.s)
  mu_coeff_White[1] = 273;  //ref_temp
  mu_coeff_White[2] = 111;  //S
  mu_transition_White = 1000;   //this number is simply picked to cover lower range

  return;
}

template <class Type>
Species<Type>::~Species()
{
}

template <class Type>
void Species<Type>::Init(std::string name, std::string database)
{
  //set the name of the species
  this->symbol = name;
  //look up info in the DB
  std::cout << "SPECIES: Reading chemical database " << database << " for species " << symbol << std::endl;
  GetDBInfo(database);
}

template <class Type>
Type Species<Type>::GetCp(Type T)
{
  Type a[7];
  GetThermoCoeff(T, a);
  //this returns Cp/R
  Type cp_R = a[0] + T*(a[1] + T*(a[2] + T*(a[3] + T*a[4])));
  //multiply by R to get Cp
  Type cp = cp_R*R;
  
  return cp; // (J/kg.K)
}

template <class Type>
Type Species<Type>::GetH(Type T, Bool shiftCurve)
{
  Type a[7];
  GetThermoCoeff(T, a);
  Type h_R = a[5] + T*(a[0] + T*(a[1]/2.0 + T*(a[2]/3.0 + T*(a[3]/4.0 + T*a[4]/5.0))));
  Type h = h_R*R;

  if(shiftCurve){
    //from the definition of these curves we need to shift h by h(298K) left and then
    //add back the heat of formation at 298K to get the actual enthalpy
    h -= href;
    //h += hf298;
  }
  
  return (h);  // J/kg
}

template <class Type>
Type Species<Type>::GetG(Type T)
{
  Type a[7];
  GetThermoCoeff(T, a);
  Type G_RT = a[0]*(1.0 - log(T)) - T*(a[1]/2.0 + T*(a[2]/6.0 + T*(a[3]/12.0 + a[4]/20.0*T))) 
    + a[5]/T - a[6];
  Type G = G_RT*R*T;
  return (G);
}

template <class Type>
Type Species<Type>::GetS(Type T)
{
  Type a[7];
  GetThermoCoeff(T, a);
  Type s_R = a[0]*log(T) + T*(a[1] + T*(a[2]/2.0 + T*(a[3]/3.0 + T*a[4]/4.0))) + a[6];
  Type s = s_R*R;

  return (s); // J/kg.K
}

template <class Type>
void Species<Type>::GetThermoCoeff(Type T, Type* a)
{
  Int i;
  if(real(T) < 200.0){
    std::stringstream ss;
    ss << "WARNING: Temperature outside of thermo valid range! Low range @ 200K";
    ss << " -- T = " << T << std::endl;
    ss << "Pinning coefficients at low range" << std::endl;
    std::cerr << ss.str();
    std::cout << ss.str();
    for(i = 0; i < 7; i++){
      a[i] = thermo_coeff[0][i];
    }
    T = 200.0; //pin low
  }
  else if(real(T) > 6000.0){
    std::stringstream ss;
    ss << "WARNING: Temperature outside of thermo valid range! High range @ 6000K";
    ss << " -- T = " << T << std::endl;
    ss << "Pinning coefficients at high range" << std::endl;
    std::cerr << ss.str();
    std::cout << ss.str();
    for(i = 0; i < 7; i++){
      a[i] = thermo_coeff[1][i];
    }
    T = 6000.0; //pin high
  }
  else{
    //use the upper curvefit
    if(real(T) > 1000.0){
      for(i = 0; i < 7; i++){
	a[i] = thermo_coeff[1][i];
      }
    }
    //use the lower curvefit
    else{
      for(i = 0; i < 7; i++){
	a[i] = thermo_coeff[0][i];
      }
    }
  }
}

template <class Type>
Int Species<Type>::GetDBInfo(std::string database)
{
  Int i = 0;
  Int err = 0;
  hid_t file_id = -1;
  std::string directory;
  std::string filename = database;

  file_id = HDF_OpenFile(filename, 0);
  if(file_id < 0){
    Abort << "Species::GetDBInfo() could not open file --" + filename;
    return file_id;
  }

  directory = "/species/";
  directory += this->symbol;

  Real temp;
  std::string variable;
  
  //read molecular weight
  variable = "MW";
  HDF_ReadScalar(file_id, directory, variable, &temp);
  this->MW = temp;
  //convert molecular weight to kg/mol
  //this saves us the headache of having factors of 1000 around
  //anytime we calculate productions, etc.
  this->MW /= 1000;
  
  //read in hf298
  variable = "NASA7_burcat_coeff15";
  HDF_ReadScalar(file_id, directory, variable, &temp);
  this->hf298 = temp;
  
  //read in curvefit coefficients
  Int n = 7;
  Real* ctemp = new Real[n];
  variable = "NASA7_burcat1";
  err = HDF_ReadArray(file_id, directory, variable, &ctemp, &n); 
  if(err){
    err = 0;
    variable = "NASA7_gupta1";
    err = HDF_ReadArray(file_id, directory, variable, &ctemp, &n);
    if(err){
      variable = "NASA7_grimech1";
      err = HDF_ReadArray(file_id, directory, variable, &ctemp, &n);
      if(err){
	std::stringstream ss;
	ss << "Could not read species coefficients " << this->symbol << std::endl;
	Abort << ss.str();
      }
    }
  }
  for(i = 0; i < n; i++) this->thermo_coeff[0][i] = ctemp[i];
  
  variable = "NASA7_burcat2";
  err = HDF_ReadArray(file_id, directory, variable, &ctemp, &n); 
  if(err){
    err = 0;
    variable = "NASA7_gupta2";
    err = HDF_ReadArray(file_id, directory, variable, &ctemp, &n); 
    if(err){
      variable = "NASA7_grimech2";
      err = HDF_ReadArray(file_id, directory, variable, &ctemp, &n);
      if(err){
	std::stringstream ss;
	ss << "Could not read species coefficients " << this->symbol << std::endl;
	Abort << ss.str();
      }
    }
  }
  for(i = 0; i < n; i++) this->thermo_coeff[1][i] = ctemp[i];
  
  variable = "NASA7_burcat2";
  err = HDF_ReadArray(file_id, directory, variable, &ctemp, &n); 
  if(err){
    err = 0;
    variable = "NASA7_gupta3";
    err = HDF_ReadArray(file_id, directory, variable, &ctemp, &n);
    if(err){
      variable = "NASA7_grimech2";
      err = HDF_ReadArray(file_id, directory, variable, &ctemp, &n);
      if(err){
	std::stringstream ss;
	ss << "Could not read species coefficients " << this->symbol << std::endl;
	Abort << ss.str();
      }
    }
  }
  for(i = 0; i < n; i++) this->thermo_coeff[2][i] = ctemp[i];
  delete [] ctemp;

  hasViscousProps = true;
  //read in thermal conductivity curvefit
  variable = "k";
  Int ncols = -1;
  Int nrows = -1;
  ctemp = NULL;
  //read to get array size
  err = HDF_ReadArray(file_id, directory, variable, &ctemp, &nrows, &ncols);
  if(err){
    std::cout << "WARNING: Thermal conductivity not found for species " << this->symbol << std::endl;
    std::cout << "WARNING: Using thermal conductivity from N2 instead" << std::endl;
    std::cerr << "WARNING: Thermal conductivity not found for species " << this->symbol << std::endl;
    std::cerr << "WARNING: Using thermal conductivity from N2 instead" << std::endl;
    std::string tempDirectory = "/species/N2";
    err = HDF_ReadArray(file_id, tempDirectory, variable, &ctemp, &nrows, &ncols);
    //allocate tempspace
    ctemp = new Real[nrows*ncols];
    //read data
    err = HDF_ReadArray(file_id, tempDirectory, variable, &ctemp, &nrows, &ncols); 
    for(i = 0; i < nrows; i++){
      for(Int j = 0; j < ncols; j++){
	k_coeff[i][j] = ctemp[i*ncols + j];
      }
    }
    k_coeff_curves = nrows;
    delete [] ctemp;
  }
  else{
    //allocate tempspace
    ctemp = new Real[nrows*ncols];
    //read data
    err = HDF_ReadArray(file_id, directory, variable, &ctemp, &nrows, &ncols); 
    for(i = 0; i < nrows; i++){
      for(Int j = 0; j < ncols; j++){
	k_coeff[i][j] = ctemp[i*ncols + j];
      }
    }
    k_coeff_curves = nrows;
    delete [] ctemp;
  }
  
  //read in visosity curvefit
  variable = "mu";
  ncols = -1;
  nrows = -1;
  ctemp = NULL;
  //read to get array size
  err = HDF_ReadArray(file_id, directory, variable, &ctemp, &nrows, &ncols);
  if(err){
    std::cout << "WARNING: Molecular viscosity not found for species " << this->symbol << std::endl;
    std::cout << "WARNING: Using molecular viscosity from N2 instead" << std::endl;
    std::cerr << "WARNING: Molecular viscosity not found for species " << this->symbol << std::endl;
    std::cerr << "WARNING: Using molecular viscosity from N2 instead" << std::endl;
    std::string tempDirectory = "/species/N2";
    //read to get array size
    err = HDF_ReadArray(file_id, tempDirectory, variable, &ctemp, &nrows, &ncols);
    //allocate tempspace
    ctemp = new Real[nrows*ncols];
    //read data
    err = HDF_ReadArray(file_id, tempDirectory, variable, &ctemp, &nrows, &ncols);
    for(i = 0; i < nrows; i++){
      for(Int j = 0; j < ncols; j++){
	mu_coeff[i][j] = ctemp[i*ncols + j];
      }
    }
    mu_coeff_curves = nrows;
    delete [] ctemp;
  }
  else{
    //allocate tempspace
    ctemp = new Real[nrows*ncols];
    //read data
    err = HDF_ReadArray(file_id, directory, variable, &ctemp, &nrows, &ncols);
    for(i = 0; i < nrows; i++){
      for(Int j = 0; j < ncols; j++){
	mu_coeff[i][j] = ctemp[i*ncols + j];
      }
    }
    mu_coeff_curves = nrows;
    delete [] ctemp;
  }
  
  //compute R
  R = UNIV_R / MW;
  
  //the value read is actually delta Hf(298)/R
  hf298 *= R;

  //compute href and dhf298
  dhf298 = GetH(298.15, false);
  href = GetH(298.15, false);
  
  //TODO: read in other important coeff. for diffusion, etc.
  
  HDF_CloseFile(file_id);
  
  return 0;
}

template <class Type>
Type Species<Type>::GetdHdT(Type T)
{
  Type a[7];
  GetThermoCoeff(T, a);
  Type dhdt_R = a[0] + T*(a[1] + T*(a[2] + T*(a[3] + T*a[4])));
  Type dhdt = dhdt_R*R;

  //we do not need to shift by H(298K) and heat of formation here b/c 
  //we are simply looking for the derivative

  return (dhdt); //(J/kg) /(K)
}

template <class Type>
void Species<Type>::Print()
{
  std::cout << symbol << std::endl;
  std::cout << "===========================================" << std::endl;
  std::cout << "Molecular weight (kg/mol): " << MW << std::endl;
  std::cout << "Specific gas constant - R (J/kg.K): " << R << std::endl;
  std::cout << "Delta hf(298K): " << dhf298 << std::endl;
  std::cout << "hf(298K): " << hf298 << std::endl;
  std::cout << "href(298K): " << href << "\n\n";
  std::cout << "T(K)\t\th (J/kg)\tCp (J/kg.K)\tk (W/m.K)\tmu (Pa.s)" << std::endl;
  std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
  std::cout.setf(std::ios::scientific);
  std::cout << 298.0 << "\t" << std::setw(8) << GetH(298.0) << "\t" << std::setw(8) << GetCp(298.0) << std::endl;
  if(hasViscousProps){
    for(Int i = 0; i <= 20; i++){
      Type T = 300.0 + i*100.0;
      std::cout << T << "\t" << std::setw(8) << GetH(T) << "\t" << std::setw(8) << GetCp(T) << 
	"\t" << std::setw(8) << GetThermalConductivity(T) << "\t" << std::setw(8) 
		<< GetViscosity(T) << std::endl; 
    }
  }
  else{
    for(Int i = 0; i <= 20; i++){
      Type T = 300.0 + i*100.0;
      std::cout << T << "\t" << std::setw(8) << GetH(T) << "\t" << std::setw(8) << GetCp(T) << std::endl; 
    }
  }
  std::cout << std::endl;
}

template <class Type>
Type Species<Type>::GetViscosity(Type T)
{
  if(real(T) <= real(mu_transition_White)){
    Type mu0 = mu_coeff_White[0];
    Type T0 = mu_coeff_White[1]; 
    Type S = mu_coeff_White[2];
    return mu0*(pow(T/T0, 1.5))*((T0 + S)/(T + S));
  }

  if(mu_coeff_curves == 0){
    Abort << "WARNING: Species::GetViscosity() has no data for visosity curvefits for species " +
      symbol;
  }

  //This is from NASA RP-1311
  //find correct range to use
  Int range = -1;
  for(Int i = 0; i < mu_coeff_curves; i++){
    if(real(T) >= real(mu_coeff[i][0]) && real(T) <= real(mu_coeff[i][1])){
      range = i;
    }
  }

  if(range == -1){
    std::stringstream ss;
    ss << T;
    Abort << "WARNING: Species::GetViscosity() out of temperature valid range for species " +
      symbol + " Temperature: " + ss.str();
  }
  Type A = mu_coeff[range][2];
  Type B = mu_coeff[range][3];
  Type C = mu_coeff[range][4];
  Type D = mu_coeff[range][5];

  Type logmu  = A*log(T) + B/T + C/(T*T) + D;

  //Result is in microPoises we need N.s/m^2 (Pa.s)
  Type convFact = 1.0e-7;
  
  return exp(logmu)*convFact; // Pa.s
}

template <class Type>
Type Species<Type>::GetThermalConductivity(Type T)
{
  if(real(T) <= real(k_transition_White)){
    Type k0 = k_coeff_White[0];
    Type T0 = k_coeff_White[1]; 
    Type S = k_coeff_White[2];
    return k0*(pow(T/T0, 1.5))*((T0 + S)/(T + S));
  }

  if(k_coeff_curves == 0){
    Abort << "WARNING: Species::GetThermalConductivity() has no data for conductivity curvefits for species " +
      symbol;
  }

  //This is from NASA RP-1311
  //find correct range to use
  Int range = -1;
  for(Int i = 0; i < k_coeff_curves; i++){
    if(real(T) >= real(k_coeff[i][0]) && real(T) <= real(k_coeff[i][1])){
      range = i;
    }
  }

  if(range == -1){
    std::stringstream ss;
    ss << T;
    Abort << "WARNING: Species::GetThermalConductivity() out of temperature valid range for species " +
      symbol + " Temperature: " + ss.str();
  }
  Type A = k_coeff[range][2];
  Type B = k_coeff[range][3];
  Type C = k_coeff[range][4];
  Type D = k_coeff[range][5];

  Type logk  = A*log(T) + B/T + C/(T*T) + D;

  //result is returned in microWatt/cm.K, we need W/m.K
  Type convFact = 0.0001;
  
  //NOTE: kinetic theory also suggests the following:
  // k = (15.0/4.0)*R_sp*mu*((4.0/15.0)*(Cp(T)/R_sp) + (1.0/3.0))
	
  return exp(logk)*convFact; // W/(m.K)
}
