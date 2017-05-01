void drawToyResults(TString inputFiles) {  

  std::vector<TString> inputResultsFiles;
  inputResultsFiles.push_back(TString("ana_pop_toys_2016_5k_GlueSmear_A0p0.root"));
  inputResultsFiles.push_back(TString("ana_pop_toys_2016_5k_GlueSmear_A0p1.root"));
  inputResultsFiles.push_back(TString("ana_pop_toys_2016_5k_GlueSmear_A0p5.root"));
  inputResultsFiles.push_back(TString("ana_pop_toys_2016_5k_GlueSmear_A1p0.root"));
  inputResultsFiles.push_back(TString("ana_pop_toys_2016_5k_GlueSmear_A1p5.root"));
  inputResultsFiles.push_back(TString("ana_pop_toys_2016_5k_GlueSmear_A2p0.root"));
  inputResultsFiles.push_back(TString("ana_pop_toys_2016_5k_GlueSmear_Am0p5.root"));
  inputResultsFiles.push_back(TString("ana_pop_toys_2016_5k_GlueSmear_Am1p0.root"));

  TGraphErrors *g_response = new TGraphErrors();
  TGraphErrors *g_pull = new TGraphErrors();

  int nPoints=0;
  for (auto inFStr : inputResultsFiles) {
    std::cout << "Processing " << inFStr << std::endl;
    
    TFile *fIn = TFile::Open(inFStr.Data());
    TTree *tree = (TTree*)fIn->Get("results");
    float a, aSigmaUp, aSigmaDn, aTrue;
    tree->SetBranchAddress("a", &a);
    tree->SetBranchAddress("aSigmaUp", &aSigmaUp);
    tree->SetBranchAddress("aSigmaDn", &aSigmaDn);
    tree->SetBranchAddress("aTrue", &aTrue);

    float mean(0.0), rms(0.0);
    float mean_pull(0.0), rms_pull(0.0);
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
      tree->GetEntry(i);
      mean += a;
      rms += a*a;
      if ((a-aTrue) > 0) {
	if (aSigmaDn > 0) {
	  mean_pull += (a-aTrue) / aSigmaDn;
	  rms_pull += (a-aTrue)*(a-aTrue)/aSigmaDn/aSigmaDn;
	} else {
	  std::cout << "Warning: invalid Dn error: " << aSigmaDn << " for asymmetry: " << a << "(aTrue = " << aTrue << ")" << std::endl;
	}
      } else {
	if (aSigmaUp > 0) {
	  mean_pull += (a-aTrue) / aSigmaUp;
	  rms_pull += (a-aTrue)*(a-aTrue)/aSigmaUp/aSigmaUp;	
	} else {
	  std::cout << "Warning: invalid Up error: " << aSigmaDn << " for asymmetry: " << a << "(aTrue = " << aTrue << ")" << std::endl;
	}
      }
    }
    mean = mean / nEntries;
    mean_pull = mean_pull / nEntries;
    rms = TMath::Sqrt(rms/nEntries - mean*mean);
    rms_pull = TMath::Sqrt(rms_pull/nEntries - mean_pull*mean_pull);
    g_response->SetPoint(nPoints, aTrue, mean);
    g_response->SetPointError(nPoints, 0, 2*rms); //convert to 2*sigma ~ 95%CL
    g_pull->SetPoint(nPoints, aTrue, mean_pull);
    g_pull->SetPointError(nPoints, 0, rms_pull); 
    std::cout << "aTrue = " << aTrue << ", mean +- rms = " << mean << " +- " << rms << std::endl;
    nPoints++;
  }

  //Draw final plot
  TCanvas *c = new TCanvas();
  g_response->Draw("AP");
  g_response->GetHistogram()->GetXaxis()->SetTitle("True asymmetry");
  g_response->GetHistogram()->GetYaxis()->SetTitle("<Measured asymmetry>");
  g_response->GetHistogram()->SetTitle("Average response and approximate 95% C.L. band");
  g_response->SetMarkerStyle(20);
  c->SetGridx();
  c->SetGridy();
  c->Update();

  TCanvas *cp = new TCanvas();
  g_pull->Draw("AP");
  g_pull->GetHistogram()->GetXaxis()->SetTitle("True asymmetry");
  g_pull->GetHistogram()->GetYaxis()->SetTitle("Pull mean #pm rms");
  g_pull->GetHistogram()->SetTitle("Pulls as function of true asymmetry");
  g_pull->SetMarkerStyle(20);
  cp->SetGridx();
  cp->SetGridy();
  cp->Update();

}
