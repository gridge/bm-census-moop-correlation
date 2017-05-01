/* Exdcutable for generation of toy experiments on persons location on map.
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


#include "MoopMap.h"
#include "MoopMapPic.h"

#include "DataDefs.h"

//uncomment to enable DEBUG
//#define DEBUG_BUILD
#include "Utilities.h"

/**********************
 * Globals
 **********************/
/// Toy settings
toySettings_t settings;

/// Instance for output data storage
data_t outputData;

//Usage function
void usage() {
  std::cerr << "Usage: generateToys [options] outputFile.root" << std::endl;
  std::cerr << "Generate toys of population position correlated with observabes to study" << std::endl;
  std::cerr << "List of options:" << std::endl;
  std::cerr << " --moop: Moop map file to use" << std::endl;
  std::cerr << " --ntoys: Set total number of toy generations" << std::endl;
  std::cerr << " --population: Set total population for each toy" << std::endl;
  std::cerr << " --frac-virgin: Set fraction of population that are virgin-burners" << std::endl;
  std::cerr << " --asym-virgin: Set asymmetry for vigin-burner observable" << std::endl;
  std::cerr << " --smear: Set smearing method" << std::endl;  
  std::cerr << " --smear-value: Set smearing parameter value" << std::endl;  
  std::cerr << " " << std::endl;
}

/**********************
 * Main function
 **********************/
int main(int argc, char **argv) {

  std::cout << "Welcome to generateToys." << std::endl;

  //Set default parameters
  settings.nToys = 1;
  settings.population = 70000;
  settings.fracVirgins = 0.393; //fraction of nPrevBurns == 0
  settings.asymmetryVirgins = 0.0; //no asymmetry
  settings.smearing = 2; // closest intersection
  settings.smearingParameter = 0.0;
  settings.moopMapFile = new std::string("../data/MoopMap-2016.xml"); //2016 Moop map

  //Read command-line options
  static struct option long_options[] = {
    {"moop",           required_argument, 0,  'm' },
    {"ntoys",          required_argument, 0,  'n' },
    {"population",     required_argument, 0,  'p' },
    {"frac-virgin",    required_argument, 0,  'f' },
    {"asym-virgin",    required_argument, 0,  'a' },
    {"smear",          required_argument, 0,  's' },
    {"smear-value",    required_argument, 0,  'v' },
    {0,         0,                 0,  0 }
  };

  while (1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    int c;
    
    c = getopt_long(argc, argv, "m:n:p:f:a:s:v:",
		    long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'n':
      settings.nToys = atoi(optarg);
      break;
    case 'm':
      *settings.moopMapFile = optarg;
      break;
    case 'p':
      settings.population = atoi(optarg);
      break;
    case 'f':
      settings.fracVirgins = atof(optarg);
      break;
    case 'a':
      settings.asymmetryVirgins = atof(optarg);
      break;
    case 's':
      settings.smearing = atoi(optarg);
      break;
    case 'v':
      settings.smearingParameter = atof(optarg);
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
    std::cerr << "ERROR: One and only one output file required." << std::endl;
    std::cerr << std::endl;
    usage();
    return 1;
  }
  settings.outputFile = new std::string(argv[optind]);

  //Print settings
  std::cout << std::endl;
  std::cout << "Settings:" << std::endl;
  std::cout << " output file: " << *settings.outputFile << std::endl;
  std::cout << " MOOP map file: " << *settings.moopMapFile << std::endl;
  std::cout << " number of toys: " << settings.nToys << std::endl;
  std::cout << " single toy population: " << settings.population << std::endl;
  std::cout << " fraction of virgin population: " << settings.fracVirgins << std::endl;
  std::cout << " virgin asymmetry: " << settings.asymmetryVirgins << std::endl;
  std::cout << " smearing type: " << settings.smearing << std::endl;
  std::cout << " smearing parameter: " << settings.smearingParameter << std::endl;
  std::cout << "---------" << std::endl;
  std::cout << std::endl;

  //prepare output file
  std::cout << "Opening output file" << std::endl;
  TFile *fOutROOT = TFile::Open(settings.outputFile->c_str(), "NEW");
  if (! fOutROOT) {
    std::cerr << "ERROR opening output file (maybe the file exists already?)" << std::endl;
    return 1;
  }

  //Store settings
  TTree *outSettingsTree = new TTree("toySettings", "Toy Settings");
  connectToySettingsBranches(outSettingsTree, settings);
  outSettingsTree->Fill();
  outSettingsTree->Write();

  //Load MOOP map (@TODO load from xml with settings + map file instead!)
  MoopMapPic *moopMapPic; //!< Moop map
  std::cout << "Loading moop map" << std::endl;
  moopMapPic = new MoopMapPic();  
  if (not moopMapPic->loadFromXML(*settings.moopMapFile)) {
    std::cerr << "ERROR loading MOOP map: " << *settings.moopMapFile << std::endl;
    delete moopMapPic;
    return 1;
  }
  float mapMinRadius, mapMaxRadius, mapMinPhi, mapMaxPhi;
  moopMapPic->getFiducialRegion(mapMinRadius, mapMaxRadius, mapMinPhi, mapMaxPhi);

  //Cache info that are useful for toy generation (see also doc/stat.pdf)
  //WARNING: to change category assignment, needs to edit also inside the toy loop!
  std::cout << "Calculating total MOOP areas (may take a minute..)" << std::endl;
  DBG("Caching general calculation about map");
  std::map<int, float> moopAreas =  moopMapPic->getMoopArea(10.0);

  ////// HACK
  //MoopMap::streetIntersection_t filterByInt = moopMapPic->getIntersection(MoopMap::street_I, moopMapPic->getClosestRadialStreet(8, 30));
  //  moopAreas = moopMapPic->getMoopAreaNearIntersection(filterByInt, 10);
  ////// END HACK


  std::cout << "MOOP fiducial areas:" << std::endl;
  float totalArea = 0.0;
  float totalValidArea = 0.0;
  for (auto i : moopAreas) {
    totalArea += i.second;
    if (i.first != MoopMap::MOOP_NOTVALID) totalValidArea += i.second;
  }
  for (auto i : moopAreas) {
    TString moopStr, area, f1, f2;
    area.Form("%.2f*10^6 ft^2", i.second/TMath::Power(10,6));
    f1.Form("%.1f%%", i.second / totalArea * 100);
    f2.Form("%.1f%%", i.second / totalValidArea * 100);
    moopStr.Form("%9s:", MoopMap::MOOP_STR[i.first].c_str());
    std::cout << moopStr << area << " (" << f1;
    if (i.first != MoopMap::MOOP_NOTVALID) 
      std::cout << ", " << f2 << " of camping";
    std::cout << ")" << std::endl;
  }
  TString a1, a2;
  a1.Form("%.0f", totalArea/TMath::Power(10,6));
  a2.Form("%.0f", totalValidArea/TMath::Power(10,6));
  std::cout << "Total fiducial area = " << a1 << "*10^6 ft^2 ( " << a2 << " of camping )" << std::endl;
  float p_m0 = 0.0;
  for (int i = MoopMap::MOOP_GREEN; i <= MoopMap::MOOP_RED; ++i) 
    if (getMOOPBinCategory(i) == 0)
      p_m0 += moopAreas[i];
  p_m0 = p_m0 / (moopAreas[MoopMap::MOOP_GREEN]+moopAreas[MoopMap::MOOP_YELLOW]+moopAreas[MoopMap::MOOP_RED]);
  float p_m1 = 1 - p_m0;
  float p_y0 = settings.fracVirgins;
  float p_y1 = 1 - p_y0;
  float A = (2 - settings.asymmetryVirgins) / (2 + settings.asymmetryVirgins);
  float p_y1_m0 = p_y1 / p_m0 * ( 1 - (p_m1*A / (p_y0 + A*p_y1)) );
  float p_y1_m1 = p_y1*A / (p_y0 + A*p_y1);

  //Prepare data tree to hold toys output
  TTree *outTree = new TTree("data", "Toy data");
  connectDataBranches(outTree, outputData);
  
  //Start generation of toys
  std::cout << "Starting toy generation" << std::endl;
  gRandom = new TRandom3();
  gRandom->SetSeed(0); 
  #ifdef DEBUG_BUILD
  //set to some other number for reproducibility tests
  gRandom->SetSeed(2811);
  #endif
  for (unsigned int toy = 0; toy < settings.nToys; ++toy) {
    // -- fill general info
    if (toy > 0 && toy % 5 == 0) {
      std::cout << toy << " toys generated so far. Saving." << std::endl;
      outTree->AutoSave();
    }  
    outputData.generation = toy;
    // -- Now loop over individuals
    for (unsigned int individual = 0; individual < settings.population; ++individual) {
      // -- generate true position within fiducial region
      //@TODO: generate in polar coordinates for speed (p = k*dx*dy = k*r*dphi*dr)
      TVector2 truePosition;
      do {
	do {
	  float x, y;
	  x = (gRandom->Rndm()-0.5)*2*mapMaxRadius;
	  y = (gRandom->Rndm()-0.5)*2*mapMaxRadius;
	  truePosition.SetX(x);
	  truePosition.SetY(y);
	  DBG("Generated true position x,y = " << x << ", " << y << " -- r,phi=" << truePosition.Mod() << ", " << truePosition.Phi());
	} while (not moopMapPic->isFiducial(truePosition));      
	DBG("True position accepted as within fiducial region.");
	outputData.trueMOOP = moopMapPic->getValue(truePosition);
      } while (not moopMapPic->isCampingAllowed(outputData.trueMOOP));

     /// HACK
	//} while (not (moopMapPic->isCampingAllowed(outputData.trueMOOP) and (moopMapPic->getClosestIntersection(truePosition) == filterByInt)));
     /// END HACK
    
      outputData.trueRadius = truePosition.Mod();
      outputData.truePhi = truePosition.Phi();
      DBG("True position accepted as with valid MOOP value: " << outputData.trueMOOP);
      // -- smear true position
      DBG("Applying smearing");
      TVector2 recoPosition = truePosition;
      if (settings.smearing == 1) {
	float smear_x(0.0);
	float smear_y(0.0);    
	//Gaussian smearing
	smear_x = gRandom->Gaus(0, settings.smearingParameter);
	smear_y = gRandom->Gaus(0, settings.smearingParameter);
	DBG("Smearing values x,y = " << smear_x << ", " << smear_y);
	recoPosition.SetX(truePosition.X() + smear_x);
	recoPosition.SetY(truePosition.Y() + smear_y);
      } else if (settings.smearing == 2) {
	//Glue to closest valid intersection
	MoopMap::streetIntersection_t closestInt;
	closestInt = moopMapPic->getClosestIntersection(truePosition);
	recoPosition = moopMapPic->getIntersectionPosition(closestInt);
      }
      DBG("New position x,y = " << recoPosition.X() << ", " << recoPosition.Y() << " -- r,phi=" << recoPosition.Mod() << ", " << recoPosition.Phi());
      outputData.posRadius = recoPosition.Mod();
      outputData.posPhi = recoPosition.Phi();

      // -- generate observables with given true asymmetry
      // First get inclusive probability of true moop value (p_{m}) and probability of being at first burn (p_{y}(0))
      //  Note: in this model we only consider m=0 and m=1
      //        currently mapped by the function getMOOPBinCategory()
      // Now get conditional probability of p(y=1 | m) for the value of m of interest.
      // see doc/stat.pdf for further information
      float p_y1_m;
      if (getMOOPBinCategory(outputData.trueMOOP) == 0) {
	p_y1_m = p_y1_m0;
      } else { //mapped_m == 1
	p_y1_m = p_y1_m1;
      }
      float catRnd = gRandom->Rndm();
      if (catRnd < p_y1_m) {
	outputData.nPrevBurns = 1; //y=1
      } else {
	outputData.nPrevBurns = 0; //y=0
      }
      outputData.trueMOOPCategory = getMOOPBinCategory(outputData.trueMOOP);
      DBG("m-category: " << getMOOPBinCategory(outputData.trueMOOP) << ", p(y=1|m)=" << p_y1_m << ", rnd=" << catRnd << " => y-category: " << outputData.nPrevBurns);
      
      // -- Fill tree and keep track of number of toys generated
      outTree->Fill();
    } //loop over individuals
  } //loop over toys
  //finalize the writing of tree to disk
  outTree->Write();

  std::cout << "All toys generated and saved to file: " << *settings.outputFile << std::endl;

  //Clean-up
  fOutROOT->Close();

  return 0;
}
