/* Exdcutable for analyzing toy data.
 * 
 * Author: gridge
 */

#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TH1F.h"
#include "TApplication.h"


#include "MoopMap.h"
#include "MoopMapPic.h"
#include "MoopDataAnalyzer.h"

#include "DataDefs.h"

//uncomment to enable DEBUG
//#define DEBUG_BUILD
#include "Utilities.h"

/**********************
 * Globals
 **********************/
///Input file with data
std::string inputFile;

///Input toys settings
toySettings_t toySettings;

///MoopMap file
std::string moopMapFile;

///Cache file for expensive computations
std::string cacheFile;

///Output file for storing statistics
std::string outputFile;

///Number of iterations for bootstrap method
unsigned int bootstrapIterations;

///Container for a given individual data
data_t data;

///Analyzer algorithm
MoopDataAnalyzer *mca;

//Usage function
void usage() {
  std::cerr << "Usage: analyzeToyData [options] inputFile.root" << std::endl;
  std::cerr << "Generate toys of population position correlated with observabes to study" << std::endl;
  std::cerr << "List of options:" << std::endl;
  std::cerr << " --moop: Override moop map file to use (otherwise grab from meta-data)" << std::endl;
  std::cerr << " --output: Output file to store statistics" << std::endl;
  std::cerr << " --bootstrap: Bootstrap method with N iterations for asymmetry determination (slow!)" << std::endl;
  std::cerr << " --cache: Cache file for expensive calculations." << std::endl;
  std::cerr << " " << std::endl;
}

/**********************
 * Main function
 **********************/
int main(int argc, char **argv) {

  std::cout << "Welcome to analyzeToyData." << std::endl;

  //Set defaults
  mca = 0;
  bootstrapIterations = 1; //single iteration

  //Read command-line options
  static struct option long_options[] = {
    {"moop",           required_argument, 0,  'm' },
    {"output",         required_argument, 0,  'o' },
    {"bootstrap",      required_argument, 0,  'b' },
    {"cache",          required_argument, 0,  'c' },
    {0,         0,                 0,  0 }
  };

  while (1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    int c;
    
    c = getopt_long(argc, argv, "m:o:b:",
		    long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'm':
      moopMapFile = optarg;
      break;
    case 'o':
      outputFile = optarg;
      break;
    case 'b':
      bootstrapIterations = atoi(optarg);
      break;
    default:
      std::cerr << std::endl;
      std::cerr << "Unrecognized option: " << (char)c << std::endl;
      std::cerr << std::endl;
      usage();
      return 1;
    }
  }

  //get non-parameters argument
  if (optind != (argc - 1)) {
    std::cerr << std::endl;
    std::cerr << "ERROR: One and only one input file required." << std::endl;
    std::cerr << std::endl;
    usage();
    return 1;
  }
  inputFile = argv[optind];

  //Init ROOT app
  TApplication theApp("analyzeToyData", &argc, argv);

  //Open input file and get meta-data
  std::cout << "Opening input file" << std::endl;
  TFile *fInROOT = TFile::Open(inputFile.c_str());
  if (! fInROOT) {
    std::cerr << "ERROR opening input file: " << inputFile << std::endl;
    return 1;
  }
  TTree *inSettingsTree = dynamic_cast<TTree*>(fInROOT->Get("toySettings"));
  if (!inSettingsTree) {
    std::cerr << "Unable to open meta-data." << std::endl;
    return 1;
  }
  //get needed information from meta-data
  connectToySettingsBranches(inSettingsTree, toySettings, false);
  
  inSettingsTree->GetEntry(0); //retrieve settings

  if (moopMapFile.empty()) {
    moopMapFile = *toySettings.moopMapFile;
  }

  //Load MOOP map
  MoopMapPic *moopMap; //!< Moop map
  std::cout << "Loading MOOP map" << std::endl;
  moopMap = new MoopMapPic();  
  if (not moopMap->loadFromXML(moopMapFile)) {
    std::cerr << "ERROR loading MOOP map: " << moopMapFile << std::endl;
    delete moopMap;
    return 1;
  } 

  //Setup analyzer
  mca = new MoopDataAnalyzer(moopMap);  
  if (toySettings.smearing == 0) {
    mca->probEvalAlg = MoopDataAnalyzer::probFromLikelihood;
    //mca->probEvalAlg = MoopDataAnalyzer::probFromMajority;
    mca->moopEvalAlg = MoopDataAnalyzer::moopFromPosition;
  } else if (toySettings.smearing == 1) {
    //mca->probEvalAlg = MoopDataAnalyzer::probFromMajority;
    mca->probEvalAlg = MoopDataAnalyzer::probFromLikelihood;
    mca->moopEvalAlg = MoopDataAnalyzer::moopFromArea;
    mca->radiusForMoopAlg = toySettings.smearingParameter;
  } else if (toySettings.smearing == 2) {
    //mca->probEvalAlg = MoopDataAnalyzer::probFromLikelihood;
    //mca->probEvalAlg = MoopDataAnalyzer::probFromPopAsym; //use fast method for toys
    mca->probEvalAlg = MoopDataAnalyzer::probFromPopAsymSlow; //test:use accurate method
    mca->moopEvalAlg = MoopDataAnalyzer::moopFromIntersection;
  }
  

  //Print settings
  std::cout << std::endl;
  std::cout << "Settings:" << std::endl;
  std::cout << " input file: " << inputFile << std::endl;
  std::cout << " output file: " << outputFile << std::endl;
  std::cout << " MOOP map file: " << moopMapFile << std::endl;
  std::cout << " cache file: " << cacheFile << std::endl;
  std::cout << " Analyzer, location algorithm: " << mca->moopEvalAlg << std::endl;
  std::cout << " Analyzer, moop prob algorithm: " << mca->probEvalAlg << std::endl;
  std::cout << " Bootstrap iterations: " << bootstrapIterations << std::endl;
  std::cout << "---------" << std::endl;
  std::cout << std::endl;
  
  //Cache results that are very expensive to compute and, if provided, load them from cache
  if (not cacheFile.empty()) {
    //load results from cache, if present and refer to the same MOOP map we have
    //@TODO
  }
  //trigger cache refresh, if not present already, for moop map areas
  if (toySettings.smearing == 2) {
    std::cout << "Caching information on MOOP areas" << std::endl;
    moopMap->getMoopArea(10.0);
    moopMap->getMoopAreaNearIntersection(moopMap->getIntersection(MoopMap::street_A, 4*3), 10); //just a radon valid intersection (e.g. A@3:00)				       
  }

  //Prepare statistics
  TH1F *h_residuals = new TH1F("h_residuals", "Residuals;a - a^{TRUE};Entries", 100, 0.0, toySettings.asymmetryVirgins*10);
  TH1F *h_pulls = new TH1F("h_pulls", "Pulls;(a - a^{TRUE})/#sigma(a);Entries", 50, -5, 5);

  //Retrieve data tree
  TTree *inDataTree = dynamic_cast<TTree*>(fInROOT->Get("data"));
  connectDataBranches(inDataTree, data, false);

  //Loop over data
  Long64_t nEntries = inDataTree->GetEntries();
  Long64_t startEntry = 0;
  Long64_t endEntry = 0;
  unsigned int currentGeneration=0;
  unsigned int invalidMoopData=0;
  unsigned int errorsResultCalc=0;
  std::vector<analysisOutput_t> results;
  std::cout << "Looping over data" << std::endl;
  // Loop structure is a bit nasty to accomodate bootstrap method:
  // - Loop over generations
  // - Loop over bootstrap iterations
  // - Loop over events within generation

    //@TODO: current loop has generations inside bootstrap loop!!
    //Need to change next loop to only loop over relevant entries and stop when generation changes
    //Then resume to next generation -> generation "whi;e {" need to be taken outside

  while (true) { //loop over generations, interrupted by break
    mca->clear(); //full-cleaning of counters only
    startEntry = endEntry; //move to next generation
    for (unsigned int ibootstrap=0; ibootstrap < bootstrapIterations; ibootstrap++) {
      DBG("Starting bootstrap iteration #" << ibootstrap << " from entry " << startEntry);
      mca->clear(true); //soft-cleaning of counters only
      Long64_t iEntry = startEntry;
      while (true) {
	Long64_t nb = inDataTree->GetEntry(iEntry);
	if (iEntry == startEntry) {
	  //initialize
	  currentGeneration = data.generation;
	}
	if (currentGeneration != data.generation or nb<=0) {
	  //end of generation or end of all data
	  endEntry = iEntry;
	  break;
	}
	TVector2 pos;
	pos.SetMagPhi(data.posRadius, data.posPhi);
	//call analyzer algorithm
	try {
	  mca->fill(data.nPrevBurns, pos, true); //use prior-correction
	  //DBG("trueMOOP: " << data.trueMOOP << ", m-category = " << getMOOPBinCategory(data.trueMOOP));
	} catch (int e) {
	  if ((e == MoopDataAnalyzer::ERROR_INVALID_MOOP) or (e == MoopDataAnalyzer::ERROR_POS_NOT_FIDUCIAL)) {
	    //can happen in some cases, keep track for stat reasons
	    invalidMoopData++;
	  } else {
	    //non-tolerable error
	    std::cerr << "Unexpected exception thrown: " << e << ". Continuing." << std::endl;
	    continue;
	  }
	}
	++iEntry;
      } //end loop over tree

      //Analyse data collected for this generatio
      if (bootstrapIterations > 1) {
	DBG("bootstrap iteration #" << ibootstrap);
	//compute and retrieve results
	analysisOutput_t result = mca->getResults();
	DBG("Asymmetry = " << result.a << " +"<< result.aSigmaUp << " -" << result.aSigmaDn << "(true: " << toySettings.asymmetryVirgins << ")");
      }
    } //end bootstrap loop

    DBG("Finished processing generation #" << data.generation);
    analysisOutput_t result = mca->getResults();
    if (result.invalidMoopData > 0) errorsResultCalc++;
    result.aTrue = toySettings.asymmetryVirgins;
    DBG("Asymmetry = " << result.a << " +"<< result.aSigmaUp << " -" << result.aSigmaDn << "(true: " << result.aTrue << ")");
    results.push_back(result);
    float res = result.a - toySettings.asymmetryVirgins;
    h_residuals->Fill(res);
    float pull;
    if (res >= 0) pull = res / result.aSigmaDn;
    else pull = res / result.aSigmaUp;
    h_pulls->Fill(pull);

    //check if we're done
    if ((currentGeneration == 0) or (currentGeneration % 5 == 4)) {
      std::cout << "Processed " << currentGeneration+1 << " generations so far." << std::endl;
    }
    if (endEntry == nEntries) break;

    //======= HACK
    //DEBUG: 1 GENERATION ONLY
    //std::cerr << "HACK: Aborting  after first iteration." << std::endl;
    //break;
    //======= HACK
      
  } //end loop over generations

  std::cout << std::endl;
  std::cout << "Done processing data." << std::endl;
  std::cout << "N. generation: " << currentGeneration+1 << std::endl;
  std::cout << "Total population: " << nEntries << std::endl;
  std::cout << "Invalid MOOP errors: " << invalidMoopData << std::endl;
  std::cout << "Errors in result calculation: " << errorsResultCalc << std::endl;
  std::cout << std::endl;

  //prepare output file
  std::cout << "Opening output file" << std::endl;
  TFile *fOutROOT = TFile::Open(outputFile.c_str(), "RECREATE");
  if (! fOutROOT) {
    std::cerr << "ERROR opening output file" << std::endl;
    return 1;
  }
  //Write historgrams
  h_residuals->Write();
  h_pulls->Write();
  //Write results tree
  TTree *outTree = new TTree("results", "Analysis results");
  analysisOutput_t outResults;
  connectAnalysisOutputBranches(outTree, outResults, true); //connect branches for writing
  for (int idx=0; idx < results.size(); ++idx) {
    outResults = results[idx];
    outTree->Fill();
  }
  outTree->Write();

  fOutROOT->Close();
  fInROOT->Close();

  std::cout << "All done." << std::endl;

  return 0;
}



