/* Executable for checking digitized map content
 * 
 * Author: gridge
 */

#include <iostream>

#include "TROOT.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TH2I.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TASImage.h"
#include "TPaveText.h"

//uncomment to enable DEBUG
//#define DEBUG_BUILD

#include "Utilities.h"
#include "MoopMap.h"
#include "MoopMapPic.h"

//Extern variables defined in checkMap.h which are needed for ROOT dictionary generation
//Interactive map objects
extern float mainWndWidth; //!< main window default width
extern TCanvas* mainWnd; //!< main frame for display, could also consider TGMainFrame
extern TImage* imageMap;  //!< image map, could also conside TGImageMap
extern TPaveText *text; //!< Text box with coordinates
extern TBox *moopColor; //!< rectangle showing selected MOOP color
extern void onClickUpdateText(); //!< Signal slot for interactive map
extern MoopMapPic *moopMapPic; //!< Moop map

//Program valid actions
std::vector<std::string> availableActions = {"rnd", "mooparea", "int"};

//Usage function
void usage() {
  std::cout << "Usage: checkMap action MoopMapFile" << std::endl;
  std::cout << "MoopMapFile is an XML file containing the information about the MOOP map" << std::endl;
  std::cout << "Available actions below. " << std::endl;
  std::cout << " rnd Fill moop map throwing random numbers and display it" << std::endl;
  std::cout << " mooparea Display information about total and MOOP area calculated from the map" << std::endl;
  std::cout << " int Display interactive map from picture to check individual points" << std::endl;
}

//Main function
int main(int argc, char **argv) {

  std::cout << "Welcome to checkMap." << std::endl;
  if ((argc < 3) or (argc > 3)) {
    usage();
    return 1;
  }


  std::string action = argv[1];
  std::string fileName = argv[2];
  if (std::find(availableActions.begin(), availableActions.end(), action) == availableActions.end()) {
    std::cerr << "Unrecognized action: " << action << std::endl;
    usage();
    return 2;
  }
  std::cout << "Input map:" << fileName << " (action: " << action << ")" << std::endl;

  //Init ROOT app
  TApplication theApp("checkMap", &argc, argv);

  //load map
  MoopMap *moopMap(0);  
  moopMap = new MoopMapPic();
  moopMapPic = dynamic_cast<MoopMapPic*>(moopMap);
  moopMapPic->loadFromXML(fileName);
  float mapMinRadius, mapMaxRadius, mapMinPhi, mapMaxPhi;
  moopMapPic->getFiducialRegion(mapMinRadius, mapMaxRadius, mapMinPhi, mapMaxPhi);

  if (action == "rnd") {

    //fill historam with map through random sampling
    unsigned long ntrials = TMath::Power(10,6);
    std::cout << "Generating random data to populate map, using " << ntrials << " random data." << std::endl;
    gRandom = new TRandom3();
    TH2I *h_moopMap = new TH2I("h_moopMap", "MOOP map rendering", 600, -6000, 6000, 600, -6000, 6000); //20ft resolution
    for (unsigned long i=0; i < ntrials; i++) {
      if (i % 100000 == 0) std::cout << "Generated " << i << " data so far." << std::endl;
      //generate random position, but restrict to fiducial volume
      float x, y;
      TVector2 pos;
      do {	
	x = (gRandom->Rndm()-0.5)*2*mapMaxRadius;
	y = (gRandom->Rndm()-0.5)*2*mapMaxRadius;
	pos.SetX(x);
	pos.SetY(y);
	DBG("Generated true position x,y = " << x << ", " << y << " -- r,phi=" << pos.Mod() << ", " << pos.Phi());
      } while (not moopMapPic->isFiducial(pos));      
      int val = moopMap->getValue(pos);
      int xbin = h_moopMap->GetXaxis()->FindBin(x);
      int ybin = h_moopMap->GetXaxis()->FindBin(y);
      //move value of 0 to value of 4 to display it clearly
      if (val == 0) val = 4;
      h_moopMap->SetBinContent(xbin, ybin, val);
      DBG("Setting " << x << "(" << xbin << "), " << y << "(" << ybin << ") = " << val);
    }

    //draw results
    std::cout << "Displaying results." << std::endl;
    gStyle->SetOptStat(0);
    TCanvas *c_map = new TCanvas("c_map", "MOOP Map", 800, 800);
    c_map->SetFixedAspectRatio();
    Int_t colors[] =      {kWhite, kGreen+8, kYellow, kRed, kBlue};
    Double_t levels[] = {0.0,   0.5,     1.5,      2.5,    3.5,  4.0};
    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    h_moopMap->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
    h_moopMap->Draw("COLZ");
    h_moopMap->GetZaxis()->SetRangeUser(0, 4.0);
    h_moopMap->GetXaxis()->SetTitle("ft");
    h_moopMap->GetYaxis()->SetTitle("ft");
    h_moopMap->GetYaxis()->SetTitleOffset(1.5);
    c_map->Update();    
    std::cout << "All Done. Close window or double-click to exit." << std::endl;
    c_map->WaitPrimitive(); 
  }

  if (action == "mooparea") {
    std::cout << "Calculating total MOOP areas (may take a minute..)" << std::endl;
    std::map<int, float> moopAreas =  moopMapPic->getMoopArea(10.0);
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
    std::cout << "Total fiducial area = " << a1 << "*10^6 ft^2 ( " << a2 << " for camping )" << std::endl;    
  }


  //create canvas with interactive check of colors
  if (action == "int") {
    std::cout << "Opening interactive MOOP map." << std::endl;
    
    //Init interactive map components
    imageMap = TImage::Open(moopMapPic->inputFile.c_str());
    mainWndWidth = 600;
    float aspect_ratio = imageMap->GetWidth() / imageMap->GetHeight();
    float mainWndHeight=static_cast<int>(mainWndWidth/aspect_ratio);
    mainWnd = new TCanvas("InteractiveMoopMap", "Interactive MOOP Map", mainWndWidth, mainWndHeight); 
    mainWnd->SetWindowSize(mainWndWidth, mainWndHeight);
    mainWnd->SetFixedAspectRatio();     
    DBG("Set canvas window size (WxH) = (" << mainWndWidth << ", " << mainWndHeight << ")");
    imageMap->SetConstRatio(kFALSE);
    text = new TPaveText(0.01, 0.01, 0.20, 0.11, "");
    text->AddText("0 ft, 00:00");
    text->AddText("[N.A.]");
    text->SetBorderSize(0);
    text->SetFillColor(kWhite);
    moopColor = new TBox(0.20, 0.01, 0.25, 0.11);
    moopColor->SetLineColor(kBlack);
    moopColor->SetFillColor(kWhite);
    moopColor->SetFillStyle(1001);        

    //Draw image in canvas
    mainWnd->cd();
    imageMap->Draw("X");
    text->Draw();
    moopColor->Draw();
    mainWnd->AddExec("exUpdateText", "onClickUpdateText()");
    mainWnd->Modified();
    mainWnd->Update();
    
  }

  //run ROOT session
  theApp.Run();

  //clean-up
  delete moopMap;

  return 0;
}


/// Signal handler for interactive map
void onClickUpdateText()
{
  int event = mainWnd->GetEvent();
  if (event != kButton1Down) return;

  //get raw pixel x position
  int raw_px = gPad->GetEventX() - gPad->XtoAbsPixel(0);
  if (raw_px < 0) raw_px = gPad->XtoAbsPixel(0);  
  //get raw pixel y position: note that y position for ROOT is zero at the bottom (corrected later)
  int raw_py = gPad->GetEventY() - gPad->YtoAbsPixel(1);
  if (raw_py < 0) raw_py = gPad->YtoAbsPixel(1);

  //Get zoomed area in pixels (if no zoom returns image in pixels)
  unsigned int zoom_x, zoom_y, zoom_w, zoom_h;
  TObject *select = gPad->GetSelected();
  if (not select) {
    //assume no zoom
    DBG("No object selected. Assume no zoom.");
    zoom_x = zoom_y = 0;
    zoom_w = imageMap->GetWidth();
    zoom_h = imageMap->GetHeight();
  } else {
    TASImage *img = dynamic_cast<TASImage*>(select);
    if (!img) {
      DBG("Failed dynamic cast. Assume no zoom");
      zoom_x = zoom_y = 0;
      zoom_w = imageMap->GetWidth();
      zoom_h = imageMap->GetHeight();
    } else {
      //get bottom-left corner of selected area
      //note that y position for ROOT is zero at the bottom and refers to bottom-left corner.
      //Apply correction so that y=0 is at the upper edge of the image and zoom_y refer to the upper edge
      // of the zoomed area (convenient below)
      img->GetZoomPosition(zoom_x, zoom_y, zoom_w, zoom_h);
      //invert height - y
      zoom_y = imageMap->GetHeight() - zoom_y - zoom_h;
    }
  }

  //correct for image/canvas scaling and zoom (zoom coordinates given in pixel!)
  //zoom_x,y now refer to the upper-left edge of the zoomed area
  int px = static_cast<float>(raw_px) / gPad->GetWw() * zoom_w + zoom_x;
  int py = static_cast<float>(raw_py) / gPad->GetWh() * zoom_h + zoom_y;

  DBG( "Click button called at gPad X,Y = (" << gPad->GetEventX() << ", " << gPad->GetEventY() << ")");
  DBG( "Zoom ROOT x=" << zoom_x << ", y=" << zoom_y << ", w=" << zoom_w << ", h=" << zoom_h );
  DBG( "Ww = " << gPad->GetWw() << ", Wh = " << gPad->GetWh());
  DBG( "image width = " << imageMap->GetWidth() << ", height = " << imageMap->GetHeight());
  DBG( "Derived pixel position: " << px << ", " << py);
  DBG( "gPad->XtoAbsPixel(0) = " << gPad->XtoAbsPixel(0) << ", gPad->YtoAbsPixel(1) = " << gPad->YtoAbsPixel(1) );

  //find moop color for given px, py
  TVector2 mousePosition = moopMapPic->getMoopCoordinatesFromPixel(px, py);
  DBG("MOOP coordinates: " << mousePosition.X() << ", " << mousePosition.Y() << "(" << mousePosition.Mod() << ", " << mousePosition.Phi() << ")");
  //find out MOOP value
  int moopValue = moopMapPic->getValue(mousePosition);
  Color_t newColor;
  if (moopValue == MoopMap::MOOP_GREEN) newColor = kGreen;
  else if (moopValue == MoopMap::MOOP_YELLOW) newColor = kYellow;
  else if (moopValue == MoopMap::MOOP_RED) newColor = kRed;
  else if (moopValue == MoopMap::MOOP_NOTVALID) newColor = kWhite;
  moopColor->SetFillColor(newColor);
  text->Clear();
  int hrs,min;
  MoopMap::getHoursFromPhi(mousePosition.Phi(), hrs, min);
  TString strCoord;
  strCoord.Form("%.0f ft, %02d:%02d", mousePosition.Mod(), hrs, min);
  text->AddText(strCoord);
  TString str_closestIntersection;
  str_closestIntersection = "[";
  str_closestIntersection += moopMapPic->getStrIntersection(moopMapPic->getClosestIntersection(mousePosition)).c_str();
  str_closestIntersection += "]";
  text->AddText(str_closestIntersection);
  bool isFiducial = moopMapPic->isFiducial(mousePosition);
  strCoord.Form("IsFiducial: %d", isFiducial);
  text->AddText(strCoord);
  mainWnd->Modified();
  mainWnd->Update();
  DBG("R="<< mousePosition.Mod() << ", phi=" << mousePosition.Phi() << ", MOOP = " << moopValue);
  DBG("Closest intersection: " << str_closestIntersection);
}

