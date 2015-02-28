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

enum reactionType
{
  ArrheniusCoeff,
  ModArrheniusCoeff,
  GuptaModArrheniusCoeff,
  Power,
  NUM_RXN_TYPES
};

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
  //public constructor.. read data from reaction file
  Reaction(Int reactionId);

  ~Reaction();
  
  std::string name;             //descriptor of reaction -- plain text

  Int nespecies;                //number of elemental species involved
  Int nspecies;                 //number of species involved
  Int ntbodies;                 //number of rxn species which also act as catalysts
  Int ncatalysts;               //number of species acting as catalysts ONLY
  Int catalysts;                //1 - if any species act as catalysts, 0 - turn catalyst action off
  Int thirdBodies;              //1 - rxn species act as catalysts,0 - pure catalysts only
  Type* TBEff;                  //third body efficiencies including catalysts and species in reaction
  Species<Type>* species;       //pointers to species in reaction
  Int* globalIndx;              //global index of species needed by this reaction 
                                //(in order - rxn species, then catalysts and thirdbodies)

  std::string* speciesSymbols;  //list of species involved in reaction
  std::string* catSymbols;
  Type* Nup;                     //stoich. coeff. left side
  Type* Nupp;                    //stoich. coeff. right side


  //Arrhenius reaction rate curve fit coefficients
  Type A;                            
  Type EA;
  Type n;  //used in modified Arrhenius

  //Arrhenius reaction rate curve fit coefficients (backward)
  Bool backwardRateGiven;
  Type Ab; 
  Type EAb;
  Type nb;
    
  void Print();
  
  //will read reaction from file and will complete a reaction object to pass back
  Int ReadReactionFromFile(std::string fileName, Int ReactionId);
  Int ParseLine(std::string& line);
  Int ParseReaction(std::string& rxn);
  Int ParseCatalysts(std::string& line);
  void SetSpeciesPointer(Species<Type>* globalSpecies, Int globalCount);

  Int GetLocalSpeciesFromGlobal(Int speciesId);

  Type GetForwardReactionRate(Type T);     //forward reaction rate
  Type GetBackwardReactionRate(Type T);    //backward reaction rate
  Type GetEquilibriumReactionRate(Type T); //equilibrium reaction rate
  Type GetMassProductionRate(Type* rhoi, Type T, Int speciesId); 

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
