/* Implementation of Interface to analyze MOOP data.
 * 
 * Author: gridge
 * 
 */

#include "MoopDataAnalyzer.h"

#include "TMath.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/BrentMinimizer1D.h"

#include "DataDefs.h"

//uncomment to enable DEBUG
//#define DEBUG_BUILD
#include "Utilities.h"


MoopDataAnalyzer::MoopDataAnalyzer(MoopMap* moopMap)
{
  moopEvalAlg = moopFromIntersection;
  probEvalAlg = probFromPopAsymSlow;
  radiusForMoopAlg = 10.0; //just a reasonable default
  map = moopMap;
  m_minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  m_minimizer->SetPrintLevel(0);
  //m_minimizer = new ROOT::Math::BrentMinimizer1D();
  clear();
}

MoopDataAnalyzer::~MoopDataAnalyzer()
{
  if (m_minimizer) delete m_minimizer;
}

void MoopDataAnalyzer::clear(bool soft)
{
  counters[0][0] = counters[0][1] = counters[1][0] = counters[1][1] = 0.0;
  countersWSq[0][0] = countersWSq[0][1] = countersWSq[1][0] = countersWSq[1][1] = 0.0;
  m_populationByIntersection.clear();
  for (auto itIntersection : map->getListOfIntersections()) {
    m_populationByIntersection[itIntersection].resize(2);
    m_populationByIntersection[itIntersection][0] = 0; //initialize
    m_populationByIntersection[itIntersection][1] = 0; //initialize    
  }
  results.population = results.y0Population = results.y1Population = 0;
  if (not soft) {
    results.a = results.aSigmaUp = results.aSigmaDn = 0.0; 
    results.p10 = results.p10SigmaUp = results.p10SigmaDn = 0.0;
    results.p11 = results.p11SigmaUp = results.p11SigmaDn = 0.0;    
    results.invalidMoopData = 0;
    results.aTrue = 0.0;
    m_cachedPriorCorr.clear();
    for (auto itIntersection : map->getListOfIntersections()) {
      m_cachedPriorCorr[itIntersection].resize(2);
      m_cachedPriorCorr[itIntersection][0] = 1.0;
      m_cachedPriorCorr[itIntersection][1] = 1.0;
    }  
  }
}

void MoopDataAnalyzer::fill(int y, TVector2 pos, bool usePriorAsymmetry)
{
  //DBG("MCA: Filling info for position (r,phi)=" << pos.Mod() << ", " << pos.Phi() << "), y-category = " << y);
  //get moop values from position
  std::map<int, float> moopAreas;
  MoopMap::streetIntersection_t intersection;
  if (moopEvalAlg == moopFromIntersection) {    
    intersection = map->getClosestIntersection(pos);
    moopAreas = map->getMoopAreaNearIntersection(intersection, 10.0);
    m_populationByIntersection[intersection][y] += 1; //increment population
  } else if (moopEvalAlg == moopFromPosition) {
    int moopVal = map->getValue(pos);
    moopAreas[moopVal] = 1; //arbitrary constant value
  } else if (moopEvalAlg == moopFromArea) {
    moopAreas = map->getMoopArea(pos, radiusForMoopAlg, 10.0); //10ft resolution hard-coded
  } else {
    throw ERROR_INVALID_ALG;
  }

  //now increment counters
  if (probEvalAlg == probFromMajority) {
    //find most likely moop value
    int mostLikelyMoop=-1;
    float mostLikelyArea=0.0;
    DBG("Scanning areas nearest to intersection");
    for (auto itA : moopAreas) {
      DBG(itA.first << ", " << itA.second);
      if ((itA.second > mostLikelyArea) && (map->isCampingAllowed(itA.second))) {
	mostLikelyArea = itA.second;
	mostLikelyMoop = itA.first;
      }
    }
    //fill corresponding counter
    if (mostLikelyMoop < 0) {
      DBG("Invalid MOOP category. Can happen in rare cases.");
      results.invalidMoopData++;
      throw ERROR_INVALID_MOOP;
    }
    int moopCat = getMOOPBinCategory(mostLikelyMoop);
    counters[y][moopCat]++; //note: no need for prior correction in this case
    DBG("Determined MOOP category based on majority algorithm: " << moopCat);
  } else if (probEvalAlg == probFromLikelihood) {
    //fill counters with weights
    float area_m0 = 0.0;
    float area_m1 = 0.0;
    for (int i = MoopMap::MOOP_GREEN; i <= MoopMap::MOOP_RED; ++i) 
      if (getMOOPBinCategory(i) == 0)
	area_m0 += moopAreas[i];
      else
	area_m1 += moopAreas[i];
    float weight_m1 = area_m1 / (area_m0 + area_m1);
    if ((area_m0 + area_m1) > 0) {
      //calculate prior correction, if any
      float priorCorr = 1.0;
      if (usePriorAsymmetry) {
	//no prior needed for moopEvalAlg == moopFromPosition	
	if (moopEvalAlg == moopFromIntersection) {	 
	  priorCorr = m_cachedPriorCorr[intersection][y];
	} else if (moopEvalAlg == moopFromArea) {       
	  std::cerr << "MoopDataAnalyzer: Prior correction for moopEvalAlg == moopFromArea not yet implemented" << std::endl;
	  throw ERROR_INVALID_ALG;
	}
      }
      float p_m1y = weight_m1 * priorCorr;
      //if (p_m1y > 1.0) p_m1y = 1.0; // cap weight at 1.0
      counters[y][1] += p_m1y; 
      counters[y][0] += (1 - p_m1y);
      countersWSq[y][1] += p_m1y*p_m1y;
      countersWSq[y][0] += (1-p_m1y)*(1-p_m1y);
      DBG("Determined MOOP category weights for m=0,1: " << 1-weight_m1 << ", " << weight_m1 << "(priorCorr = " << priorCorr << ")");
    } else {
      DBG("Invalid MOOP category. Can happen in rare cases.");
      results.invalidMoopData++;
      throw ERROR_INVALID_MOOP;
    }
  } else if ((probEvalAlg == probFromPopAsym) or (probEvalAlg == probFromPopAsymSlow)) {
    //nothing else to do here, it's enough to count the population per intersection
  } else {
    throw ERROR_INVALID_ALG;
  }
  //update population
  results.population++;
  if (y == 0) results.y0Population++;
  else if (y==1) results.y1Population++;
}

analysisOutput_t MoopDataAnalyzer::getResults()
{
  //Calculate asymmetry and erros. note: counters[y][m] !!
  //float c_0 = counters[0][0] + counters[0][1]; //tot population for category y=0
  //float c_1 = counters[1][0] + counters[1][1]; //tot population for category y=1
  float c_0 = results.y0Population;
  float c_1 = results.y1Population;

  if ((c_1 <= 0) or (c_0 <= 0)) {    
    return results;
  }

  if (probEvalAlg == probFromMajority) {
    //compute frequentist p(m=1|y=0) and p(m=1|y=1)    
    results.p10 = counters[0][1] / c_0;
    results.p11 = counters[1][1] / c_1;
    //easy to calculate error under Gaussian approximation (-> symmetric errors)
    float c11 = counters[1][1];
    float c01 = counters[0][1];
    //Assume numerator uncertainty is the only one that counts (i.e. numerator << denominator)
    results.p10SigmaUp = TMath::Sqrt(c01) / c_0;
    results.p10SigmaDn = results.p10SigmaUp;
    results.p11SigmaUp = TMath::Sqrt(c11) / c_1;
    results.p11SigmaDn = results.p11SigmaUp;
  } else if (probEvalAlg == probFromLikelihood) {
    //get p10 first
    calculateEfficiency(results.p10, results.p10SigmaDn, results.p10SigmaUp, counters[0][1], 0, c_0, countersWSq[0][0] + countersWSq[0][1]);
    //std::cout << "p10: tw=" << tw << ", tw2=" << tw2 << ", pw=" << pw <<" -> p10=" << results.p10 << std::endl;
    //then get p11
    calculateEfficiency(results.p11, results.p11SigmaDn, results.p11SigmaUp, counters[1][1], 0, c_1, countersWSq[1][0] + countersWSq[1][1]);
    //std::cout << "p11: tw=" << tw << ", tw2=" << tw2 << ", pw=" << pw <<" -> p11=" << results.p11 << std::endl;
  } else if ((probEvalAlg == probFromPopAsym) or (probEvalAlg == probFromPopAsymSlow)) {
    //get common things first
    std::map<int, float> moopAreas =  map->getMoopArea(10.0);
    float area_m0 = 0.0;
    float area_m1 = 0.0;
    for (int i = MoopMap::MOOP_GREEN; i <= MoopMap::MOOP_RED; ++i) 
      if (getMOOPBinCategory(i) == 0)
	area_m0 += moopAreas[i];
      else
	area_m1 += moopAreas[i];
    float pm0 = area_m0 / (area_m0 + area_m1);
    float pm1 = 1 - pm0;
    float py1 = c_1 / results.population;
    float py0 = 1 - py1;
    //now loop over intersections and store results for each
    std::vector<float> p10_i;
    std::vector<float> p10ErrUp_i;
    std::vector<float> p10ErrDn_i;
    std::vector<TF1*> p10_posterior_i;
    std::vector<std::pair<float, float>> p10_posterior_i_pars;
    std::vector<float> p11_i;
    std::vector<float> p11ErrUp_i;
    std::vector<float> p11ErrDn_i;
    std::vector<TF1*> p11_posterior_i;
    std::vector<std::pair<float, float>> p11_posterior_i_pars;
    std::vector<float> a_i;
    std::vector<float> aErrUp_i;
    std::vector<float> aErrDn_i;
    DBG("Results per interaction:");
    for (auto i : map->getListOfIntersections()) {
      float pm0_i = 0.0;
      float pm1_i = 0.0;
      float py1_i = 0.0;
      float py1_i_errUp = 0.0;
      float py1_i_errDn = 0.0;
      float py0_i = 0.0;
      float py0_i_errUp = 0.0;
      float py0_i_errDn = 0.0;
      //get intersecion-specific info
      moopAreas = map->getMoopAreaNearIntersection(i, 10.0);
      area_m0 = 0.0;
      area_m1 = 0.0;
      for (int i = MoopMap::MOOP_GREEN; i <= MoopMap::MOOP_RED; ++i) 
	if (getMOOPBinCategory(i) == 0)
	  area_m0 += moopAreas[i];
	else
	  area_m1 += moopAreas[i];
      pm1_i = area_m1 / (area_m0 + area_m1);
      pm0_i = 1 - pm1_i;
      auto popByInt = m_populationByIntersection[i];      
      calculateEfficiency(py0_i, py0_i_errUp, py0_i_errDn, popByInt[0], popByInt[0]+popByInt[1]);
      calculateEfficiency(py1_i, py1_i_errUp, py1_i_errDn, popByInt[1], popByInt[0]+popByInt[1]);

      //calculate probabilities, asymmetry and uncertainties
      float p10 = (py0_i / py0 - pm0_i / pm0) / (pm1_i / pm1 - pm0_i / pm0); 
      float p10ErrUp = (py0_i_errUp / py0) / (pm1_i / pm1 - pm0_i / pm0);
      float p10ErrDn = (py0_i_errDn / py0) / (pm1_i / pm1 - pm0_i / pm0);      
      p10_i.push_back(p10);      
      p10ErrUp_i.push_back(p10ErrUp);
      p10ErrDn_i.push_back(p10ErrDn);
      TF1* p10_pdf = new TF1((TString("p10_pdf_")+p10_i.size()).Data(), "TMath::BetaDist(x,[0],[1])",0.0, 1.0);
      p10_pdf->SetParameter(0, 1.0+popByInt[0]); //@TODO does not account for possible weights, OK for now
      p10_pdf->SetParameter(1, 1.0+popByInt[1]);
      p10_posterior_i.push_back(p10_pdf);
      p10_posterior_i_pars.push_back(std::make_pair<float,float>(1.0 / (py0 * ( pm1_i / pm1 - pm0_i / pm0 )), -1.0*(pm0_i / pm0) / (pm1_i / pm1 - pm0_i / pm0) ));
//      #ifdef DEBUG_BUILD
//      c->cd(1);
//      p10_pdf->Draw();
//      #endif

      float p11 = (py1_i / py1 - pm0_i / pm0) / (pm1_i / pm1 - pm0_i / pm0); 
      float p11ErrUp = (py1_i_errUp / py1) / (pm1_i / pm1 - pm0_i / pm0);
      float p11ErrDn = (py1_i_errDn / py1) / (pm1_i / pm1 - pm0_i / pm0);
      p11_i.push_back(p11);
      p11ErrUp_i.push_back(p11ErrUp);
      p11ErrDn_i.push_back(p11ErrDn);
      TF1* p11_pdf = new TF1((TString("p11_pdf_")+p11_i.size()).Data(), "TMath::BetaDist(x,[0],[1])",0.0, 1.0);
      p11_pdf->SetParameter(0, 1.0+popByInt[1]); //@TODO does not account for possible weights, OK for now
      p11_pdf->SetParameter(1, 1.0+popByInt[0]);
      p11_posterior_i.push_back(p11_pdf);
      p11_posterior_i_pars.push_back(std::make_pair<float,float>(1.0 / (py1 * ( pm1_i / pm1 - pm0_i / pm0 )), -1.0*(pm0_i / pm0) / (pm1_i / pm1 - pm0_i / pm0) ));
//      #ifdef DEBUG_BUILD
//      c->cd(2);
//      p11_pdf->Draw();
//      #endif

      float a = 2 * (p10 - p11) / (p10 + p11);
      float aErrUp = 4 * TMath::Sqrt( (TMath::Power(p11*p10ErrUp,2) + TMath::Power(p10*p11ErrUp,2)) ) / (TMath::Power(p11+p10,2));
      float aErrDn = 4 * TMath::Sqrt( (TMath::Power(p11*p10ErrDn,2) + TMath::Power(p10*p11ErrDn,2)) ) / (TMath::Power(p11+p10,2));
      a_i.push_back(a);
      aErrUp_i.push_back(aErrUp);
      aErrDn_i.push_back(aErrDn);
      DBG("pm0_i=" << pm0_i << " p10 = " << p10 << "(+-" << p10ErrUp << ","<<p10ErrDn << "), p11 = " << p11 << "(+-" << p11ErrUp << ","<<p11ErrDn << "), a = " << a << "(+-" << aErrUp << ","<< aErrDn << ")");      
      //      #ifdef DEBUG_BUILD
      //c->WaitPrimitive();
      //      #endif
    } //end loop over intersections

    // Perform combination of likelihoods for extracting probabilities
    //--- Using Minuit2
    if (probEvalAlg == probFromPopAsymSlow) {
      if (not calculateProbSlow(results.p11, results.p11SigmaUp, results.p11SigmaDn, p11_posterior_i, p11_posterior_i_pars)) {
	DBG("Error in minimization, keep track.");
	results.invalidMoopData++;
      }
      if (not calculateProbSlow(results.p10, results.p10SigmaUp, results.p10SigmaDn, p10_posterior_i, p10_posterior_i_pars)) {
	DBG("Error in minimization, keep track.");
	results.invalidMoopData++;
      }
    } else { //probFromPopAsym
      if (not calculateProb(results.p11, results.p11SigmaUp, results.p11SigmaDn, p11_posterior_i, p11_posterior_i_pars)) {
	DBG("Error in minimization, keep track.");
	results.invalidMoopData++;
      }
      if (not calculateProb(results.p10, results.p10SigmaUp, results.p10SigmaDn, p10_posterior_i, p10_posterior_i_pars)) {
	DBG("Error in minimization, keep track.");
	results.invalidMoopData++;
      }
    }

    //Check for failures (weak check, but that's ok for now)
    if ((results.p11SigmaUp < 0.00001) or (results.p11SigmaDn < 0.00001)) {
      //Do not believe it, error estimation failed 
      //Assign similar relative error as p10
      DBG("Patching failed error estimation on p11 with p10 relative errors");
      results.p11SigmaUp = results.p10SigmaUp / results.p10 * results.p11;
      results.p11SigmaDn = results.p10SigmaDn / results.p10 * results.p11;
    } else if ((results.p10SigmaUp < 0.00001) or (results.p10SigmaDn < 0.00001)) {
      //assign similar relative error
      DBG("Patching failed error estimation on p10 with p11 relative errors");
      results.p10SigmaUp = results.p11SigmaUp / results.p11 * results.p10;
      results.p10SigmaDn = results.p11SigmaDn / results.p11 * results.p10;
    }

    // Cleanup
    for (auto f : p10_posterior_i) delete f;
    for (auto f : p11_posterior_i) delete f;
  }
  //calculate asymmetry 
  results.a = 2 * (results.p10 - results.p11) / (results.p10 + results.p11);
  results.aSigmaUp = 4 * TMath::Sqrt( (TMath::Power(results.p11*results.p10SigmaUp,2) + TMath::Power(results.p10*results.p11SigmaUp,2)) ) / (TMath::Power(results.p11+results.p10,2));
  results.aSigmaDn = 4 * TMath::Sqrt( (TMath::Power(results.p11*results.p10SigmaDn,2) + TMath::Power(results.p10*results.p11SigmaDn,2)) ) / (TMath::Power(results.p11+results.p10,2));
  DBG("RESULTS: p10 = " << results.p10 << "(+-" << results.p10SigmaUp << ","<<results.p10SigmaDn << "), p11 = " << results.p11 << "(+-" << results.p11SigmaUp << ","<< results.p11SigmaDn << "), a = " << results.a << "(+-" << results.aSigmaUp << ","<< results.aSigmaDn << ")");

  //store in cache results for prior correction (in case is needed) for each possible intersection
  // Note: we HAVE to cache this, since it depends on population in sensitive area that is only known after first iteration
  if ((moopEvalAlg == moopFromIntersection) && (probEvalAlg == probFromLikelihood)) {
    m_cachedPriorCorr.clear();
    for (auto popByInt : m_populationByIntersection) {
      float pop = popByInt.second[0]+popByInt.second[1];
      if (pop == 0) continue; //not needed anyway
      float py0 = static_cast<float>(popByInt.second[0]) / pop;
      float py1 = 1 - py0;
      float A = (2-results.a) / (2+results.a);
      m_cachedPriorCorr[popByInt.first].resize(2);
      m_cachedPriorCorr[popByInt.first][0] = 1.0 / (py0 + A*py1);
      m_cachedPriorCorr[popByInt.first][1] = A / (py0 + A*py1);
    }
  } else if ((moopEvalAlg == moopFromArea) && (probEvalAlg == probFromLikelihood)) {
    //@TODO if needed
    throw ERROR_INVALID_ALG;
  }

  return results;
}

void MoopDataAnalyzer::calculateEfficiency(float &efficiency, float &effErrUp, float &effErrDn, float pw, float tw, float pw2, float tw2, bool poissonRatio)
{
  //use same treatment as in TGraphAsymErrors::Divide, with bayesian treatment and CL=68.3%
  //https://root.cern.ch/doc/master/TGraphAsymmErrors_8cxx_source.html#l00858
  //@TODO allow poisson-ratio case!
  double aa, bb;
  double norm;
  double alpha = 1.0;
  double beta = 1.0;
  double conf = 0.683; //1-sigma of normal distribution
  if (tw2<=0) tw2 = tw;
  if (pw2<=0) pw2 = pw;
  
  norm = tw/tw2;
  aa = pw * norm + alpha;
  bb = (tw - pw) * norm + beta;
  efficiency = TEfficiency::BetaMode(aa,bb);
  effErrDn = results.p10 - TEfficiency::BetaCentralInterval(conf,aa,bb,false);
  effErrUp = TEfficiency::BetaCentralInterval(conf,aa,bb,true) - results.p10;
  return;
}


double MoopDataAnalyzer::pMLogLikelihood(const double *x)
{
  double val = 0.0; 
  for (int ipdf = 0; ipdf < m_counts_posterior_i.size(); ++ipdf) {
    double xprime = (x[0] - m_prob_posterior_i_pars[ipdf].second) / m_prob_posterior_i_pars[ipdf].first;
    if ((xprime > 0) and (xprime < 1))
      val = val - TMath::Log(m_counts_posterior_i[ipdf]->Eval(xprime));
  } 
  return val;
}

bool MoopDataAnalyzer::calculateProb(float &prob, float &probErrUp, float &probErrDn, const std::vector<TF1*> &counts_posterior_i, const std::vector<std::pair<float, float>> &prob_posterior_i_pars)
{
  bool status = true;
  //--- Using Minuit2
  double eup, edn;
  ROOT::Math::Functor f_ll(this, &MoopDataAnalyzer::pMLogLikelihood, 1);
  //ROOT::Math::Functor1D f_ll(this, &MoopDataAnalyzer::pMLogLikelihood);
  m_minimizer->SetFunction(f_ll);
  //m_minimizer->SetTolerance(0.1);

  m_counts_posterior_i = counts_posterior_i;
  m_prob_posterior_i_pars = prob_posterior_i_pars;
  m_minimizer->SetVariable(0, "f_ll", 0.1, 0.01);
  if (not m_minimizer->Minimize()) {
    DBG("Failed first minimization.");
    //try to recover
    m_minimizer->SetVariable(0, "f_ll", 0.1, 0.0001);
    m_minimizer->SetVariableLimits(0, 0.0, 1.0);
    if (not m_minimizer->Minimize()) {
      std::cout << "Failed second minimization too. Revert to slow minimization." << std::endl;
      status = calculateProbSlow(prob, probErrUp, probErrDn, counts_posterior_i, prob_posterior_i_pars);
      return status;
    }
  }
  prob = m_minimizer->X()[0];
  status = status or m_minimizer->GetMinosError(0, eup, edn);
  probErrUp = eup;
  probErrDn = edn;
  return status;
}

bool MoopDataAnalyzer::calculateProbSlow(float &prob, float &probErrUp, float &probErrDn, const std::vector<TF1*> &counts_posterior_i, const std::vector<std::pair<float, float>> &prob_posterior_i_pars)
{
  //--- Using TF1
  //Define -log(likelihood) combining intersections
  TF1 * prob_pdf = new TF1("prob_pdf",[&](double*x, double *p){ 
      double val = 0.0; 
      for (int ipdf = 0; ipdf < counts_posterior_i.size(); ++ipdf) {
	double xprime = (x[0] - prob_posterior_i_pars[ipdf].second) / prob_posterior_i_pars[ipdf].first;
	if ((xprime > 0) and (xprime < 1))
	  val = val - TMath::Log(counts_posterior_i[ipdf]->Eval(xprime));
      } 
      return val;
    }, 0.0, 1.0, 0);
  prob_pdf->SetNpx(10000);
  Double_t q[2];
  Double_t minXLL=0.0;
  prob = prob_pdf->GetMinimumX();
  double minLL = prob_pdf->Eval(prob);
  q[0] = prob_pdf->GetX(minLL+0.5);
  q[1] = prob_pdf->GetX(minLL+0.5, prob);
  probErrDn = prob - q[0];
  probErrUp = q[1] - prob;
#ifdef DEBUG_BUILD
  TCanvas *c = new TCanvas();
  prob_pdf->Draw();
  c->WaitPrimitive();
  delete c;
#endif
  delete prob_pdf;

  if ((probErrUp/prob < 0.00001) or (probErrDn/prob < 0.00001)) {
    return false;
  }
  return true;
}
