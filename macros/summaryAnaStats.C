// Create summary table/graphs for various data analyses
void summaryAnaStats(TString inputFiles)
{
  TChain *resultsTree = new TChain("results");
  resultsTree->Add(inputFiles);

  float a, aSigmaUp, aSigmaDn;
  unsigned int population, y0Population;
  std::string *category = new std::string();
  resultsTree->SetBranchAddress("category", &category);
  resultsTree->SetBranchAddress("a", &a);
  resultsTree->SetBranchAddress("aSigmaUp", &aSigmaUp);
  resultsTree->SetBranchAddress("aSigmaDn", &aSigmaDn);
  resultsTree->SetBranchAddress("population", &population);
  resultsTree->SetBranchAddress("y0Population", &y0Population);

  Long64_t nEntries = resultsTree->GetEntries();
  //format: category, population, fraction y==0, asymmetry, a uncertainty, significance, A
  TString formatCols  = "%-20s, %10d, %4.1f, %+5.2f, +%4.2f -%4.2f, %4.1f, %5.2f";
  TString formatHeader= "%-20s, %10s, %4s, %5s, %11s, %4s, %5s";
  std::cout << TString::Format(formatHeader, "Category", "Population", "% y0", "a", "s(a)", "sig.", "A") << std::endl;  
  for (Long64_t i=0; i< nEntries; ++i) {
    resultsTree->GetEntry(i);
    float significance; //very rough
    if (a > 0) significance = a / aSigmaDn;
    else significance = TMath::Abs(a / aSigmaUp);
    std::cout << TString::Format(formatCols, 
				 category->c_str(), population, float(y0Population)/population*100,
				 a, aSigmaUp, aSigmaDn, 
				 significance, (2-a)/(2+a)) << std::endl;				 
  }
  std::cout << std::endl;
}
