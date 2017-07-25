/* Collect data utility functions
 * 
 * Author: gridge
 */

#include "DataDefs.h"

#include "TTree.h"

int getMOOPBinCategory(int moopValue) 
{
  if (moopValue == MoopMap::MOOP_GREEN) return 0;
  else if ((moopValue == MoopMap::MOOP_YELLOW) or (moopValue == MoopMap::MOOP_RED)) return 1;
  //if ((moopValue == MoopMap::MOOP_GREEN) or (moopValue == MoopMap::MOOP_YELLOW)) return 0;
  //else if (moopValue == MoopMap::MOOP_RED) return 1;
  return -1;
}

void connectToySettingsBranches(TTree *tree, toySettings_t &instance, bool writing)
{
  if (writing) {
    tree->Branch("outputFile", instance.outputFile); //need a pointer to correctly read/write 
    tree->Branch("moopMapFile", instance.moopMapFile);
    tree->Branch("nToys", &instance.nToys, "nToys/i");
    tree->Branch("population", &instance.population, "population/i");
    tree->Branch("fractionVirgins", &instance.fracVirgins, "fractionVirgins/F");
    tree->Branch("asymmetryVirgins", &instance.asymmetryVirgins, "asymmetryVirgins/F");
    tree->Branch("smearing", &instance.smearing, "smearing/i");
    tree->Branch("smearingParameter", &instance.smearingParameter, "smearingParameter/F");  
  } else {    
    tree->SetBranchAddress("outputFile", &instance.outputFile); //need a pointer to correctly read/write 
    tree->SetBranchAddress("moopMapFile", &instance.moopMapFile);
    tree->SetBranchAddress("nToys", &instance.nToys);
    tree->SetBranchAddress("population", &instance.population);
    tree->SetBranchAddress("fractionVirgins", &instance.fracVirgins);
    tree->SetBranchAddress("asymmetryVirgins", &instance.asymmetryVirgins);
    tree->SetBranchAddress("smearing", &instance.smearing);
    tree->SetBranchAddress("smearingParameter", &instance.smearingParameter);  
  }
}

void connectDataBranches(TTree *tree, data_t &instance, bool writing)
{
  if (writing) {
    tree->Branch("generation", &instance.generation, "generation/i");
    tree->Branch("trueRadius", &instance.trueRadius, "trueRadius/F");
    tree->Branch("truePhi", &instance.truePhi, "truePhi/F");
    tree->Branch("trueMOOP", &instance.trueMOOP, "trueMOOP/I");
    tree->Branch("trueMOOPCategory", &instance.trueMOOPCategory, "trueMOOPCategory/I");
    tree->Branch("posRadius", &instance.posRadius, "posRadius/F");
    tree->Branch("posPhi", &instance.posPhi, "posPhi/F");
    tree->Branch("nPrevBurns", &instance.nPrevBurns, "nPrevBurns/I");    
  } else {
    tree->SetBranchAddress("generation", &instance.generation);
    tree->SetBranchAddress("trueRadius", &instance.trueRadius);
    tree->SetBranchAddress("truePhi", &instance.truePhi);
    tree->SetBranchAddress("trueMOOP", &instance.trueMOOP);
    tree->SetBranchAddress("trueMOOPCategory", &instance.trueMOOPCategory);
    tree->SetBranchAddress("posRadius", &instance.posRadius);
    tree->SetBranchAddress("posPhi", &instance.posPhi);
    tree->SetBranchAddress("nPrevBurns", &instance.nPrevBurns);        
  }
}

void connectAnalysisOutputBranches(TTree *tree, analysisOutput_t &instance, bool writing)
{
  if (writing) {
    tree->Branch("a", &instance.a, "a/F");
    tree->Branch("aSigmaUp", &instance.aSigmaUp, "aSigmaUp/F");
    tree->Branch("aSigmaDn", &instance.aSigmaDn, "aSigmaDn/F");
    tree->Branch("aTrue", &instance.aTrue, "aTrue/F");
    tree->Branch("pm1y0", &instance.p10, "pm1y0/F");
    tree->Branch("pm1y0SigmaUp", &instance.p10SigmaUp, "pm1y0SigmaUp/F");
    tree->Branch("pm1y0SigmaDn", &instance.p10SigmaDn, "pm1y0SigmaDn/F");
    tree->Branch("pm1y1", &instance.p11, "pm1y1/F");
    tree->Branch("pm1y1SigmaUp", &instance.p11SigmaUp, "pm1y1SigmaUp/F");
    tree->Branch("pm1y1SigmaDn", &instance.p11SigmaDn, "pm1y1SigmaDn/F");
    tree->Branch("population", &instance.population, "population/i");
    tree->Branch("y0Population", &instance.y0Population, "y0Population/i");
    tree->Branch("y1Population", &instance.y1Population, "y1Population/i");
    tree->Branch("invalidMoopData", &instance.invalidMoopData, "invalidMoopData/i");
  } else {
    tree->SetBranchAddress("a", &instance.a);
    tree->SetBranchAddress("aSigmaUp", &instance.aSigmaUp);
    tree->SetBranchAddress("aSigmaDn", &instance.aSigmaDn);
    tree->SetBranchAddress("aTrue", &instance.aTrue);
    tree->SetBranchAddress("pm1y0", &instance.p10);
    tree->SetBranchAddress("pm1y0SigmaUp", &instance.p10SigmaUp);
    tree->SetBranchAddress("pm1y0SigmaDn", &instance.p10SigmaDn);
    tree->SetBranchAddress("pm1y1", &instance.p11);
    tree->SetBranchAddress("pm1y1SigmaUp", &instance.p11SigmaUp);
    tree->SetBranchAddress("pm1y1SigmaDn", &instance.p11SigmaDn);
    tree->SetBranchAddress("population", &instance.population);
    tree->SetBranchAddress("y0Population", &instance.y0Population);
    tree->SetBranchAddress("y1Population", &instance.y1Population);
    tree->SetBranchAddress("invalidMoopData", &instance.invalidMoopData);
  }
}
