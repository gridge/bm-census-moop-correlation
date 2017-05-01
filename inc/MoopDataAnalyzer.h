/* Class to analyze MOOP data.
 * 
 * Author: gridge
 * 
 */

#ifndef __MOOPDATAANALYZER_H__
#define __MOOPDATAANALYZER_H__

#include <vector>
#include <map>

#include "Math/Minimizer.h"
#include "Math/BrentMinimizer1D.h"
#include "TF1.h"

#include "MoopMap.h"
#include "DataDefs.h"

/** General interface for MOOP data analysis.  */
class MoopDataAnalyzer {
 public:
  /** Constructor.
   * @param moopMap pointer to the moop map
   */
  MoopDataAnalyzer(MoopMap *moopMap);

  ~MoopDataAnalyzer();

  /** Analyze individual and keep statistics.
   * @param y y-category (e.g. virgin)
   * @param pos position in MOOP coordinates
   * @param map pointer to MOOP map
   * @param usePriorAsymmetry use previous result for asymmetry as prior (used in bootstrap method)
   */
  virtual void fill(int y, TVector2 pos, bool usePriorAsymmetry=false);

  /** Clear analysis data and start over. 
   * @param soft only clear counters (for iterative bootstrap determination)
   */
  virtual void clear(bool soft=false);

  ///Compute and retrieve results
  virtual analysisOutput_t getResults();

 public:
  /// Config: select algorithm for MOOP determination from position 
  enum moopEvalAlg_t {
    moopFromPosition = 0,
    moopFromArea = 1,
    moopFromIntersection = 2,
    nMoopEvalAlgs
  } moopEvalAlg;

  /// Config: algorithm for probability evaluation
  enum probEvalAlg_t {
    probFromMajority = 0,
    probFromLikelihood = 1,
    probFromPopAsym = 2,
    probFromPopAsymSlow = 3,
    nProbEvalAlgs
  } probEvalAlg;

  ///Config: radius for moopFromArea algorithm
  float radiusForMoopAlg;

 protected:
  ///Output results structure
  analysisOutput_t results;

 protected:

  // Internal pointer to moop map
  MoopMap* map;

  /** Internal counters (could be weighted) for probability calculation. 
   * indices: [y, m]
   * y = y-category (e.g. virgin or not)
   * m = m-category (e.g. GREEN moop or YELLOW+R
   */
  float counters[2][2];

  ///Keep track of weights squared, if needed
  float countersWSq[2][2];

  ///Cache for prior correction. For each intersection store prior correction for each y-category value
  std::map<MoopMap::streetIntersection_t, std::vector<float>> m_cachedPriorCorr;

  ///Store population per-intersection for prior correction, if needed. Each vector entry is a y-category
  std::map<MoopMap::streetIntersection_t, std::vector<unsigned int>> m_populationByIntersection;

 private:
  /** Calculate bayesian efficiency
   * @param efficiency returns efficiency valie
   * @param effErrUp unceertainty: up
   * @param effErrUp unceertainty: down
   * @param pw Sum of (weighted) counts for "passed" (numerator)
   * @param tw Sum of (weighted) counts for "total" (denominator)
   * @param pw2 Sum of (weighted) sqaured counts for "passed" (numerator). If 0, assumes pw2=pw
   * @param tw2 Sum of (weighted) squared counts for "total" (denominator). If 0, assumes tw2=tw
   * @param poissonRatio If true, treat numerator/denominator as poisson variables
   */
  void calculateEfficiency(float &efficiency, float &effErrUp, float &effErrDn, float pw, float tw, float pw2=0, float tw2=0, bool poissonRatio=false);

  ///Minimizer for probability calculation
  ROOT::Math::Minimizer* m_minimizer;
  //ROOT::Math::BrentMinimizer1D *m_minimizer;

  /** Calculate p(m | y) combining likelihood for each intersection given population p(y|i) and parameters (p(y|i)=p(y|i) * first + second
   * @param prob return probability p(y|m,i)
   * @param probErrUp return probability uncertainty
   * @param probErrDn return probability uncertainty
   * @param counts_posterior_i pdf for each intersection: p(y|i)
   * @param prob_posterior_i_pars first and second intersection-dependent parameters to go from p(y|i) to p(m|y)
   */
  bool calculateProb(float &prob, float &probErrUp, float &probErrDn, const std::vector<TF1*> &counts_posterior_i, const std::vector<std::pair<float, float>> &prob_posterior_i_pars);

  /** Slow (but robust) method to calculate p(m | y) combining likelihood for each intersection given population p(y|i) and parameters (p(y|i)=p(y|i) * first + second
   * @param prob return probability p(y|m,i)
   * @param probErrUp return probability uncertainty
   * @param probErrDn return probability uncertainty
   * @param counts_posterior_i pdf for each intersection: p(y|i)
   * @param prob_posterior_i_pars first and second intersection-dependent parameters to go from p(y|i) to p(m|y)
   */
  bool calculateProbSlow(float &prob, float &probErrUp, float &probErrDn, const std::vector<TF1*> &counts_posterior_i, const std::vector<std::pair<float, float>> &prob_posterior_i_pars);

  ///Minimization function for probability calculation (-log(L))
  double pMLogLikelihood(const double *x);

  ///Input parameters for probability log likelihood
  std::vector<TF1*> m_counts_posterior_i;
  /// Parameters to derive probability from counts: prob = counts * first + second
  std::vector<std::pair<float, float>> m_prob_posterior_i_pars;

 public:
  // Error codes

  static const int ERROR_INVALID_ALG = 1;
  static const int ERROR_INVALID_MOOP = 1;

};


#endif
