/* Collect data structures for toy and real data.
 * 
 * Author: gridge
 */

#ifndef __DATA_STRUCTS_H__
#define __DATA_STRUCTS_H__

#include <string>

#include <TTree.h>

#include "MoopMap.h"

/**********************
 * data file format
 **********************/
struct data_t {

  ///Toy generation progressive index
  unsigned int generation; 

  ///True position generated: radius (MOOP coordinates)
  float trueRadius;

  ///True position generated: phi (MOOP coordinates)
  float truePhi;

  ///True value of MOOP at generated position
  int trueMOOP;

  ///True value of MOOP category at generated position
  int trueMOOPCategory;

  ///Position recorded: radius (MOOP coordinates)
  float posRadius;
  
  ///Position recorded: phi (MOOP coordinates)
  float posPhi;

  ///observable generated: number of previous burns (0 = virgin, >=1 = not virgin)
  int nPrevBurns;
  
};

/// Structure to keep settings (will be saved into file)
struct toySettings_t {
  /// Output ROOT file with generated toys
  std::string *outputFile;

  /// MOOP map used
  std::string *moopMapFile;

  /// Number of toy generations
  unsigned int nToys;

  /// Population for each toy
  unsigned int population; 

  /// Fraction of virgin burners
  float fracVirgins;

  /* MOOP asymmetry in virgin observable
   * In the schema of
   *  m = 0(no MOOP), 1 (MOOP)
   *  y = 0(virgin), 1 (experienced)
   * asymmetry is defined as:
   * a = 2 * (p(m=1 | y=0) - p(m=1 | y=1)) / (p(m=1 | y=0) + p(m=1 | y=1))
   * i.e. the difference of probabilities to produce MOOP for a virgin and an experienced person
   *      divided by their average.
   * See doc/stat.pdf for further information
   */
  float asymmetryVirgins; 

  /** Type of smearing applied 
   * 0: no smearing applied
   * 1: gaussian smearing (argument is radius in feets)
   * 2: glue to closest-intersection (argument introduce an "error" probability)
   */
  unsigned int smearing;

  /// Smearing parameter (optional, see above)
  float smearingParameter;

};

/// Output information from population analysis
struct analysisOutput_t {

  ///Category being analyzed
  std::string *category;

  /// Asymmetry
  float a;
  ///asymmetry uncertainty
  float aSigmaUp;
  ///asymmetry uncertainty
  float aSigmaDn;

  ///For toys, store true value of asymmetry as well (for convenience)
  float aTrue;

  ///Prob(m=1 | y=0)
  float p10;
  ///p10 uncertainty
  float p10SigmaUp;
  ///p10 uncertainty
  float p10SigmaDn;

  ///Prob(m=1 | y=1)
  float p11;
  ///p10 uncertainty
  float p11SigmaUp;
  ///p10 uncertainty
  float p11SigmaDn;

  ///Population (excludes invalid data)
  unsigned int population;
  ///Population category y=0
  unsigned int y0Population;
  ///Population category y=1
  unsigned int y1Population;

  ///counter for invalid moop errors
  unsigned int invalidMoopData;
  
};


/**********************
 * Helpers
 **********************/
/** Define MOOP binary category
 * Assumes already validated input!
 * Alternatives could be to cluster GREEN,YELLOW as "0" and RED as "1"
 */
int getMOOPBinCategory(int moopValue);

/// Connect branches for toySettings_t
void connectToySettingsBranches(TTree *tree, toySettings_t &instance, bool writing=true);

/// Connect branches for data_t
void connectDataBranches(TTree *tree, data_t &instance, bool writing=true);

/// Connect branches for analysisOutput_t
void connectAnalysisOutputBranches(TTree *tree, analysisOutput_t &instance, bool writing=true);

#endif
