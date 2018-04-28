#ifndef CHEM_REACTION_H__
#define CHEM_REACTION_H__

#include "chem.h"
#include "species.h"
#include "chem_constants.h"
#include "strings_util.h"
#include "exceptions.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>

enum reactionType
{
  ArrheniusCoeff,
  ModArrheniusCoeff,
  GuptaModArrheniusCoeff,
  Power,
  NUM_RXN_TYPES
};

// Arrhenius k = A exp(-Ea/RT)
// ModArrhenius k = A T^n exp(-Ea/RT)
// GuptaModArrhenius  k = A T^n exp(-Ea/T)
// Power k = A T^n

const std::string rxnTypeNames[] = 
  {
    "Arrhenius form",        
    "modified Arrhenius form",
    "Gupta modified Arrhenius form",
    "power form",
    //this must match the above enumeration --- be warned!!!!
  };

template <class Type>
class Reaction
{
 public:

  Reaction();
  ~Reaction();
  
  std::string name;             //descriptor of reaction -- plain text

  Int nespecies;                //number of elemental species involved
  bool thirdBodiesPresent;      //1 - third body reaction, 0 - non-third body reaction
  std::vector<Type> TBEff;      //third body efficiencies including species in reaction
  Species<Type>* species;       //pointers to species in reaction, list includes species that show up as third bodies as well (reaction does NOT own this pointer)
  std::vector<Int> globalIndx;  //list of indices for every reaction species which indx's into the global list
  
  std::vector<std::string> speciesSymbols;  //list of species involved in reaction
  std::vector<Type> Nup;                     //stoich. coeff. left side
  std::vector<Type> Nupp;                    //stoich. coeff. right side

  //Arrhenius reaction rate curve fit coefficients
  Type A;                            
  Type EA;
  Type n;  //used in modified Arrhenius

  //Arrhenius reaction rate curve fit coefficients (backward)
  Bool backwardRateGiven;
  Type Ab; 
  Type EAb;
  Type nb;

  Int GetNspecies();
  void Print();
  
  //will read reaction from file and will complete a reaction object to pass back
  Int ReadReactionFromFile(std::string fileName, Int ReactionId);
  Int ReadReactionFromFileChemkin(std::string fileName, Int ReadtionId);
  Int ParseReactionBlockLine(std::string& line);
  void ParseReaction(std::string& rxn);
  void ParseThirdBodies(std::string& line);

  void SetSpeciesPointer(Species<Type>* globalSpecies, Int globalCount);
  Int GetLocalSpeciesFromGlobal(Int speciesId);

  Type GetForwardReactionRate(Type T);     //forward reaction rate
  Type GetBackwardReactionRate(Type T);    //backward reaction rate
  Type GetEquilibriumReactionRate(Type T); //equilibrium reaction rate
  Type GetMassProductionRate(Type* rhoi, Int globalNspecies, Type T, Int speciesId); 

  Int rxnType;   //see enum. above for types
  Int rxnTypeBackward; //see enum. above for types
            
  Type Arrhenius(Type A, Type EA, Type T);
  Type ModArrhenius(Type A, Type EA, Type n, Type T);
  Type GuptaModArrhenius(Type A, Type EA, Type n, Type T);
  Type PowerRxn(Type A, Type n, Type T);

  //returns string of rxn formatted for human readability
  std::string GetFormattedReaction();

 private:

};

//include implementations
#include "reaction.tcc"

#endif
