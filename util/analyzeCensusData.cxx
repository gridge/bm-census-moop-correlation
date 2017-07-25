/* Exdcutable for generation of toy experiments on persons location on map.
 * 
 * Author: gridge
 */

#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TApplication.h"

#include "MoopMap.h"
#include "MoopMapPic.h"
#include "MoopDataAnalyzer.h"

#include "DataDefs.h"

//uncomment to enable DEBUG
//#define DEBUG_BUILD
#include "Utilities.h"

//Usage function
void usage() {
  std::cerr << "Usage: analyzeCensusData [options] moopMapFile.xml inputFile.tsv" << std::endl;
  std::cerr << "Analyze input TSV (filtered with macros/parseCensusCSV.R) Census data." << std::endl;
  std::cerr << "List of options:" << std::endl;
  std::cerr << " --output: Output file to store statistics" << std::endl;
  std::cerr << " " << std::endl;
}

/**********************
 * Main function
 **********************/
int main(int argc, char **argv) {

  std::cout << "Welcome to analyzeCensusData." << std::endl;

  std::string outputFile;

  //Read command-line options
  static struct option long_options[] = {
    {"help",           no_argument, 0,  'h' },
    {"output",         required_argument, 0,  'o' },
    {0,         0,                 0,  0 }
  };

  while (1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    int c;
    
    c = getopt_long(argc, argv, "ho:",
		    long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'o':
      outputFile  = optarg;
      break;
    case 'h':
      usage();
      return 0;
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
  if (optind != (argc - 2)) {
    std::cerr << std::endl;
    std::cerr << "ERROR: One and only one input file required." << std::endl;
    std::cerr << std::endl;
    usage();
    return 1;
  }
  std::string moopMapFile(argv[optind]);
  std::string inputFile(argv[optind+1]);


  //Print settings
  std::cout << std::endl;
  std::cout << "Settings:" << std::endl;
  std::cout << " input file: " << inputFile << std::endl;
  std::cout << " MOOP Map file: " << moopMapFile << std::endl;
  std::cout << " output file: " << outputFile << std::endl;
  std::cout << "---------" << std::endl;
  std::cout << std::endl;

  //Init ROOT app
  TApplication theApp("analyzeCensusData", &argc, argv);

  //Open input file and create TTree
  TTree *tData = new TTree("data", "data");
  tData->ReadFile(inputFile.c_str(), "id/I:address_letter/C:address_hour/I:virgin_bool/I:weightnerds/F");
  int id;
  tData->SetBranchAddress("id", &id);
  char address_letter[10];
  tData->SetBranchAddress("address_letter", &address_letter);
  int address_hour;
  tData->SetBranchAddress("address_hour", &address_hour);
  int virgin_bool;
  tData->SetBranchAddress("virgin_bool", &virgin_bool);
  float weightnerds;
  tData->SetBranchAddress("weightnerds", &weightnerds);

  //Open MOOP Map
  MoopMapPic *moopMap; //!< Moop map
  std::cout << "Loading MOOP map" << std::endl;
  moopMap = new MoopMapPic();  
  if (not moopMap->loadFromXML(moopMapFile)) {
    std::cerr << "ERROR loading MOOP map: " << moopMapFile << std::endl;
    delete moopMap;
    return 1;
  } 
  std::cout << "Caching information on MOOP areas" << std::endl;
  moopMap->getMoopArea(10.0);
  moopMap->getMoopAreaNearIntersection(moopMap->getIntersection(MoopMap::street_A, 4*3), 10); //just a radon valid intersection (e.g. A@3:00)				       

  //Setup analyzer
  std::cout << "Loading MOOP data analyzer" << std::endl;
  MoopDataAnalyzer *mca = new MoopDataAnalyzer(moopMap);  
  mca->probEvalAlg = MoopDataAnalyzer::probFromPopAsymSlow; 
  mca->moopEvalAlg = MoopDataAnalyzer::moopFromIntersection;

  //Prepare results structure
  analysisOutput_t outResults;  
  TTree *outTree = new TTree("results", "Analysis results");
  connectAnalysisOutputBranches(outTree, outResults, true); //connect branches for writing  

  //Analyze data
  std::cout << "Analyzing data" << std::endl;
  Long64_t nEntries = tData->GetEntries();
  for (int iEntry=0; iEntry < nEntries; ++iEntry) {
    if ((iEntry % 1000 == 0) && (iEntry > 0))
      std::cout << "Processed " << iEntry << " entries." << std::endl;
    if (tData->GetEntry(iEntry) < 0) break;
    //prepare data
    TVector2 pos;
    MoopMap::roundStreet_t roundStreet = moopMap->getRoundStreet(address_letter);
    int hours, mins;
    mins = address_hour % 100;
    hours = address_hour / 100; //integer division
    MoopMap::radialStreet_t radialStreet = moopMap->getClosestRadialStreet(hours, mins);
    if (moopMap->isIntersectionForbidden(moopMap->getIntersection(roundStreet, radialStreet)))
      continue; //entry not valid
    pos = moopMap->getIntersectionPosition(roundStreet, radialStreet);
    int nBurns;
    if (virgin_bool == 1) nBurns = 0;
    else nBurns = 1;
    mca->fill(nBurns, pos);
  }

  std::cout << "Finished analyzing data." << std::endl;
  outResults = mca->getResults();  

  //display results
  std::cout << std::endl << "---- Results: " << std::endl;
  std::cout << "Original entries: " << nEntries << std::endl;
  std::cout << "Valid population: " << outResults.population << ", (y==0) fraction: " << (float)outResults.y0Population / outResults.population << std::endl;
  std::cout << "p10 = " << outResults.p10 << " +" << outResults.p10SigmaUp << " -" << outResults.p10SigmaDn << std::endl;
  std::cout << "p11 = " << outResults.p11 << " +" << outResults.p11SigmaUp << " -" << outResults.p11SigmaDn << std::endl;
  std::cout << "a = " << outResults.a << " +" << outResults.aSigmaUp << " -" << outResults.aSigmaDn << std::endl;
  std::cout << "---" << std::endl << std::endl;

  //write output file
  if (not outputFile.empty()) {
    std::cout << "Opening output file" << std::endl;
    TFile *fOutROOT = TFile::Open(outputFile.c_str(), "RECREATE");
    if (! fOutROOT) {
      std::cerr << "ERROR opening output file" << std::endl;
      return 1;
    }
    //write ttree
    //data->Write(); //for DEBUG only
    outTree->Fill();
    outTree->Write();
    fOutROOT->Close();
  }
  
  //Quit
  std::cout << "All Done." << std::endl;
  return 0;
}
