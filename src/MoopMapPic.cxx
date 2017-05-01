/* Implementation of MoopMapPic class
 *
 * Author: gridge
 */

#include "MoopMapPic.h"

#include "TMath.h"

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

//ROOT dictionary
ClassImp(MoopMapPic)

//uncomment for DEBUG of this specific class
//#define DEBUG_BUILD
#include "Utilities.h"

//Initialization of static data member
int MoopMapPic::nInstances = 0;
const std::vector<unsigned char> MoopMapPic::MOOP_GREEN_RGB = {0, 166, 80};
const std::vector<unsigned char> MoopMapPic::MOOP_YELLOW_RGB = {255, 255, 0};
const std::vector<unsigned char> MoopMapPic::MOOP_RED_RGB = {255, 0, 0};
const std::vector<unsigned char> MoopMapPic::MOOP_BLACK_RGB = {255, 0, 0};
const std::vector<unsigned char> MoopMapPic::MOOP_WHITE_RGB = {255, 255, 255};
 

MoopMapPic::MoopMapPic()
{
  //Init SILLY library
  if (MoopMapPic::nInstances == 0) {
    MoopMapPic::nInstances++;
    SILLY::SILLYInit();
  }
  m_sizePixel = 0.0;

  m_pos_center.Set(0.0,0.0);
  m_pos_scale = 0.0; 
  m_tolerance_RGB = 50;
}

MoopMapPic::~MoopMapPic()
{
  if (m_pictureLocation) delete m_pictureLocation;
  if (m_moopPicture) delete m_moopPicture;

  if (nInstances == 1) {
    nInstances--;
    SILLY::SILLYCleanup();
  }
}

int MoopMapPic::getValue(TVector2 position)
{
  DBG("Requested MOOP value for position x="<<position.X() << ", y=" << position.Y());
  int moopValue = MOOP_NOTVALID;

  //check if we are in a valida region
  size_t idx = ToPicBufferIdx(position);  
  if (idx < 0) {
    DBG("Requested position outside image range.");
    return MOOP_NOTVALID;
  }  

  //check if we have a match spot-on (one pixel)
  moopValue = getSinglePixelColor(position);
  DBG("Single-Pixel MOOP = " << moopValue);
  if (moopValue != MOOP_UNKNOWN) {
    //got it!
    return moopValue;
  }

  //else we can't tell from this pixel alone, use average of window
  //@TODO: make tolerance configurable, for now hard-code 20 feets
  moopValue = getProminentPixelColor(position, 20.0);
  DBG("Prominent-Pixel MOOP = " << moopValue);
  //do not allow MOOP_UNKNOWN however (internal of this implementation)
  if (moopValue == MOOP_UNKNOWN) moopValue = MOOP_NOTVALID;

  return moopValue;
}

bool MoopMapPic::load(std::string pic, TVector2 positionCenter, TVector2 positionOtherRef, TVector2 refPoint)
{
  DBG("Opening image " << pic);
  DBG("Position center = (" << positionCenter.X() << ", " << positionCenter.Y() << ")");
  DBG("Position other = (" << positionOtherRef.X() << ", " << positionOtherRef.Y() << ")");
  DBG("reference Point = (" << refPoint.X() << ", " << refPoint.Y() << ")");
  m_pictureLocation = new SILLY::FileDataSource(pic.c_str());
  m_moopPicture = new SILLY::Image(*m_pictureLocation);
  if (not m_moopPicture->loadImageHeader()) {
    DBG("ERROR opening header for image: " << pic);
    return false;;
  }
  m_moopPicture->loadImageData(SILLY::PF_RGB); //load 24bit pixels (R, G, B)
  m_sizePixel = 3; //3-bytes per-pixel in PF_RGB
  if (not m_moopPicture->isValid()) {
    throw ERROR_LOAD_MAP;
  }

  m_pos_center = positionCenter;
  m_pos_scale =  (positionOtherRef - positionCenter).Mod() / refPoint.Mod();

  DBG("m_pos_center = (" << m_pos_center.X() << ", " << m_pos_center.Y() << "), m_pos_scale = " << m_pos_scale);
  DBG("Map image and reference frame loaded");
  return true;
}

bool MoopMapPic::loadFromXML(std::string mapXMLFile)
{
  xmlDocPtr doc;
  xmlNodePtr cur;
  DBG("Opening XML file " << mapXMLFile);
  doc = xmlParseFile(mapXMLFile.c_str());
  
  if (doc == NULL ) {
    std::cerr << "Error in loading document: " << mapXMLFile << std::endl;
    return false;
  }
  
  cur = xmlDocGetRootElement(doc);
  
  if (cur == NULL) {
    std::cerr << "Error in loading document (empty document):" << mapXMLFile << std::endl;
    xmlFreeDoc(doc);
    return false;
  }
  
  if (xmlStrcmp(cur->name, (const xmlChar *) "MoopMap")) {
    std::cerr << "Input map of the wrong type, root node != MoopMap" << std::endl;
    xmlFreeDoc(doc);
    return false;
  }
  
  //now parse all information
  TVector2 centerPixelPosition;
  TVector2 refPixelPosition;
  TVector2 refPosition;
  bool refFrameLoaded=false;
  bool roundStreetLocationsLoaded=false;

  cur = cur->xmlChildrenNode;
  DBG("Parsing XML tree");
  while (cur != NULL) {
    if (cur->type != XML_ELEMENT_NODE) {
      cur = cur->next;
      continue;
    }
    DBG(" Element: '" << cur->name << "', type: " << cur->type);
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"version"))) {
      xmlChar *key;
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      xmlVersion = reinterpret_cast<char*>(key);
      xmlFree(key);
      DBG("Retrieved xmlVersion: " << xmlVersion);
    } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"year"))) {
      xmlChar *key;
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      mapYear = reinterpret_cast<char*>(key);
      xmlFree(key);
      DBG("Retrieved mapYear: " << mapYear);      
    } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"Data"))) {
      xmlNodePtr childNode = cur->xmlChildrenNode;
      while (childNode != NULL) {
	if (childNode->type == XML_ELEMENT_NODE) {
	  if (!xmlStrcmp(childNode->name, (const xmlChar *)"type")) {
	    xmlChar *key;
	    key = xmlNodeListGetString(doc, childNode->xmlChildrenNode, 1);
	    inputType=reinterpret_cast<char*>(key);
	    DBG("Retrieved data.type: " << inputType);
	    if (xmlStrcmp(key, (const xmlChar*)"picture")) {
	      //if not "picture", we don't know how to handle it for now
	      std::cerr << "MOOP file with unknown data.type: " << key << std::endl;
	      xmlFree(key);
	      xmlFreeDoc(doc);
	      return false;
	    }
	    xmlFree(key);
	  } else if (!xmlStrcmp(childNode->name, (const xmlChar*)"file")) {
	      xmlChar *key;
	      key = xmlNodeListGetString(doc, childNode->xmlChildrenNode, 1);
	      inputFile = reinterpret_cast<char*>(key);
	      DBG("Retrieved data.inputFile map: " << inputFile);
	      xmlFree(key);
	  }
	} //node == XML_ELEMENT_NODE
	childNode = childNode->next;
      }
    } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"ReferenceFrame"))){
      xmlNodePtr childNode = cur->xmlChildrenNode;
      while (childNode != NULL) {
	DBG("node " << childNode->name << ", type=" << childNode->type);
	if (childNode->type == XML_ELEMENT_NODE) {
	  if (!xmlStrcmp(childNode->name, (const xmlChar *)"pixelPositionCenterMap")) {
	    int x,y;
	    bool s=true;
	    s = s and getXMLProperty(childNode, "x", x);
	    s = s and getXMLProperty(childNode, "y", y);
	    if (!s) return false; //error loading
	    centerPixelPosition.Set((float)x, (float)y);
	    DBG("Retrieved ReferenceFrame.pixelPositionCenterMap = " << centerPixelPosition.X() << ", " << centerPixelPosition.Y());
	  } else if (!xmlStrcmp(childNode->name, (const xmlChar *)"pixelPositionOnMap")) {
	    int x,y,hh,mm;
	    float r;
	    bool s = true;
	    s = s and getXMLProperty(childNode, "x", x);
	    s = s and getXMLProperty(childNode, "y", y);
	    s = s and getXMLProperty(childNode, "radius", r);
	    s = s and getXMLProperty(childNode, "hours", hh);
	    s = s and getXMLProperty(childNode, "minutes", mm);
	    if (!s) return false; //error loading
	    refPixelPosition.Set((float)x,(float)y);
	    refPosition.SetMagPhi(r, getPhiFromHours(hh, mm));
	    DBG("Retrieved ReferenceFrame.pixelPositionOnMap. (x,y) = (" << refPixelPosition.X() << ", " << refPixelPosition.Y()
		<< "), MOOP coordinates (r,phi = (" << refPosition.Mod() << ", " << refPosition.Phi() << ")");
	  } else {
	    std::cerr << "WARNING: Ignoring un-recognized XML tag inside ReferenceFrame: " << childNode->name << std::endl;
	  }
	} //type == XML_ELEMENT_NODE
	childNode = childNode->next;
      } //loop over childs      
      refFrameLoaded=true;
    } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"FiducialRegion"))) {
      xmlNodePtr childNode = cur->xmlChildrenNode;
      while (childNode != NULL) {	
	DBG("node " << childNode->name << ", type=" << childNode->type);
	if (childNode->type == XML_ELEMENT_NODE) {
	  if (!xmlStrcmp(childNode->name, (const xmlChar *)"radius")) {	    
	    float minR,maxR;
	    bool s=true;
	    s = s and getXMLProperty(childNode, "min", minR);
	    s = s and getXMLProperty(childNode, "max", maxR);
	    if (!s) return false; //error loading
	    DBG("Retrieved FiducialRegion.Radius = " << minR << " - " << maxR);
	    setFiducialRadii(minR, maxR);
	  } else if (!xmlStrcmp(childNode->name, (const xmlChar *)"angle")) {
	    int hhMin, hhMax, mmMin, mmMax;
	    bool s = true;
	    s = s and getXMLProperty(childNode, "hoursMin", hhMin);
	    s = s and getXMLProperty(childNode, "hoursMax", hhMax);
	    s = s and getXMLProperty(childNode, "minutesMin", mmMin);
	    s = s and getXMLProperty(childNode, "minutesMax", mmMax);
	    if (!s) return false; //error loading
	    DBG("Retrieved FiducialRegion.angle: " << hhMin << ":" << mmMin << " - " << hhMax << ":" << mmMax);
	    setFiducialPhi(getPhiFromHours(hhMin, mmMin), getPhiFromHours(hhMax, mmMax));
	  } else if (!xmlStrcmp(childNode->name, (const xmlChar *)"excludeCircle")) {
	    float centerRadius, circleRadius;
	    int centerHH, centerMM;
	    bool s = true;
	    s = s and getXMLProperty(childNode, "centerRadius", centerRadius);
	    s = s and getXMLProperty(childNode, "centerHours", centerHH);
	    s = s and getXMLProperty(childNode, "centerMinutes", centerMM);
	    s = s and getXMLProperty(childNode, "circleRadius", circleRadius);
	    if (!s) return false; //error loading
	    DBG("Retrieved Excluded circular FiducialRegion: (r,HH:MM)=" << centerRadius << "," << centerHH << ":" << centerMM << "), r = " << centerRadius);
	    setFiducialExclusionCircle(centerRadius, getPhiFromHours(centerHH, centerMM), circleRadius);
	  } else {
	    std::cerr << "WARNING: Ignoring un-recognized XML tag inside FiducialRegion: " << childNode->name << std::endl;
	  }
	} //type == XML_ELEMENT_NODE
	childNode = childNode->next;
      } //loop over childs  
    } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"RoundStreetsLocation"))) {
      xmlNodePtr childNode = cur->xmlChildrenNode;
      while (childNode != NULL) {	
	DBG("node " << childNode->name << ", type=" << childNode->type);
	if (childNode->type == XML_ELEMENT_NODE) {
	  if (!xmlStrcmp(childNode->name, (const xmlChar *)"street")) {	    
	    std::string name;
	    float radius;
	    bool s=true;
	    s = s and getXMLProperty(childNode, "name", name);
	    s = s and getXMLProperty(childNode, "radius", radius);
	    if (!s) return false; //error loading
	    MoopMap::roundStreet_t street = getRoundStreet(name);
	    if (street == street_NRoundStreets) {
	      std::cerr << "Invalid street name: " << name << std::endl;
	      return false; //do not ignore these errors!
	    }
	    m_streetRadius[street] = radius;	    
	    DBG("Retrieved Street location: " << name << " = " << radius);
	  } else {
	    std::cerr << "WARNING: Ignoring un-recognized XML tag inside RoundStreetsLocation: " << childNode->name << std::endl;
	  }
	} //type == XML_ELEMENT_NODE
	childNode = childNode->next;
      } //loop over childs  
    } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"ForbiddenIntersections"))) {
      xmlNodePtr childNode = cur->xmlChildrenNode;
      while (childNode != NULL) {	
	DBG("node " << childNode->name << ", type=" << childNode->type);
	if (childNode->type == XML_ELEMENT_NODE) {
	  if (!xmlStrcmp(childNode->name, (const xmlChar *)"intersection")) {	    
	    std::string round;
	    int radialHH, radialMM;
	    bool s=true;
	    s = s and getXMLProperty(childNode, "round", round);
	    s = s and getXMLProperty(childNode, "radialHours", radialHH);
	    s = s and getXMLProperty(childNode, "radialMinutes", radialMM);
	    if (!s) return false; //error loading
	    MoopMap::roundStreet_t roundSt = getRoundStreet(round);
	    if (roundSt == street_NRoundStreets) {
	      std::cerr << "Invalid street name: " << round << std::endl;
	      return false; //do not ignore these errors!
	    }
	    addForbiddenIntersection(getIntersection(roundSt, getClosestRadialStreet(radialHH, radialMM)));
	    DBG("Added forbidden intersections: " << round << " @ " << radialHH << ":" << radialMM);
	  } else if (!xmlStrcmp(childNode->name, (const xmlChar *)"intersectionRadialRange")) {
	    std::string minRound, maxRound;
	    int radialHH, radialMM;
	    bool s=true;
	    s = s and getXMLProperty(childNode, "minRound", minRound);
	    s = s and getXMLProperty(childNode, "maxRound", maxRound);
	    s = s and getXMLProperty(childNode, "radialHours", radialHH);
	    s = s and getXMLProperty(childNode, "radialMinutes", radialMM);
	    if (!s) return false; //error loading
	    MoopMap::roundStreet_t minRoundSt = getRoundStreet(minRound);
	    MoopMap::roundStreet_t maxRoundSt = getRoundStreet(maxRound);
	    if (minRoundSt == street_NRoundStreets) {
	      std::cerr << "Invalid street name: " << minRound << std::endl;
	      return false; //do not ignore these errors!
	    }
	    if (maxRoundSt == street_NRoundStreets) {
	      std::cerr << "Invalid street name: " << maxRound << std::endl;
	      return false; //do not ignore these errors!
	    }
	    if (minRoundSt > maxRoundSt) {
	      //swap them
	      MoopMap::roundStreet_t tmp;
	      tmp = minRoundSt;
	      minRoundSt = maxRoundSt;
	      maxRoundSt = tmp;
	    }
	    MoopMap::roundStreet_t street = minRoundSt;
	    maxRoundSt = getNextRoundStreet(maxRoundSt); //to allow inclusive range
	    do {
	      addForbiddenIntersection(getIntersection(street, getClosestRadialStreet(radialHH, radialMM)));
	      DBG("Added forbidden intersections: " << roundStreetNames[street] << " @ " << radialHH << ":" << radialMM);
	      street = getNextRoundStreet(street);
	    } while (street != maxRoundSt);
	  } else {
	    std::cerr << "WARNING: Ignoring un-recognized XML tag inside ForbiddenIntersections: " << childNode->name << std::endl;
	  }
	} //type == XML_ELEMENT_NODE
	childNode = childNode->next;
      } //loop over childs        
    } else {
      std::cerr << "WARNING: Ignoring unrecognized XML tag:'" << cur->name << "' in MOOP map file: " << mapXMLFile << std::endl;
    }
    cur = cur->next;
  }  
  
  xmlFreeDoc(doc);

  //Now load map
  if (refFrameLoaded and not inputFile.empty()) {
    load(inputFile, centerPixelPosition, refPixelPosition, refPosition);
  } else {
    std::cerr << "Not enough information to load map! Need input file and reference frame at least." << std::endl;
    return false;
  }

  return true;

}

void MoopMapPic::setRGBTolerance(float newTolerance)
{
  m_tolerance_RGB = newTolerance;
}

TVector2 MoopMapPic::getMoopCoordinatesFromPixel(size_t pixel_x, size_t pixel_y)
{
  DBG("Converting (" << pixel_x << ", " << pixel_y << ").");
  DBG("m_pos_scale = " << m_pos_scale << ", m_pos_center = (" << m_pos_center.X() << ", " << m_pos_center.Y() << ")");
  TVector2 pixelCoord(pixel_x, pixel_y);
  TVector2 moopCoords = (pixelCoord - m_pos_center) / m_pos_scale;
  return moopCoords;
}

std::pair<size_t, size_t> MoopMapPic::ToPicCoordinates(TVector2 position)
{
  std::pair<size_t, size_t> imageCoord;
  TVector2 tmpImageCoord = m_pos_center + position*m_pos_scale;
  imageCoord.first = static_cast<size_t>(tmpImageCoord.X());
  imageCoord.second = static_cast<size_t>(tmpImageCoord.Y());
  return imageCoord;
}

size_t MoopMapPic::ToPicBufferIdx(TVector2 position)
{
  std::pair<size_t, size_t> imageCoord = ToPicCoordinates(position);
  return ToPicBufferIdx(imageCoord.first, imageCoord.second);
}

size_t MoopMapPic::ToPicBufferIdx(size_t pixel_x, size_t pixel_y)
{
  size_t idx = (m_moopPicture->getWidth()*pixel_y) + pixel_x;
  idx = idx * m_sizePixel; //multiply by size of pixel data
  if (idx >= m_moopPicture->getPixelsDataSize()) idx = -1; //invalid index
  return idx;
}

float MoopMapPic::getColorDistance(std::vector<unsigned char> color1, std::vector<unsigned char> color2)
{
  float distance = 0.0;
  for (unsigned int i=0; i < 3; ++i) {
    distance += TMath::Power(color1[i] - color2[i],2);
  }  
  return TMath::Sqrt(distance);
}

int MoopMapPic::getSinglePixelColor(TVector2 position)
{
  std::pair<size_t, size_t> imageCoord;
  imageCoord = ToPicCoordinates(position);
  return getSinglePixelColor(imageCoord.first, imageCoord.second);
}

int MoopMapPic::getSinglePixelColor(size_t pixel_x, size_t pixel_y)
{
  size_t pixelIdx = ToPicBufferIdx(pixel_x, pixel_y);
  std::vector<unsigned char> RGB = {0, 0, 0};
  const unsigned char *p_pixel = 0;
  p_pixel = m_moopPicture->getPixelsDataPtr() + pixelIdx;
  RGB[0] = *p_pixel;
  RGB[1] = *(p_pixel+1);
  RGB[2] = *(p_pixel+2);
  //Hierarchical test: RED, YELLOW, GREEN
  if (getColorDistance(MOOP_RED_RGB, RGB) < m_tolerance_RGB) {
    return MOOP_RED;
  } else if (getColorDistance(MOOP_YELLOW_RGB, RGB) < m_tolerance_RGB) {
    return MOOP_YELLOW;
  } else if (getColorDistance(MOOP_GREEN_RGB, RGB) < m_tolerance_RGB) {
    return MOOP_GREEN;
  } else if (getColorDistance(MOOP_WHITE_RGB, RGB) < m_tolerance_RGB) {
    return MOOP_NOTVALID;
  }
  //else returns that we can't tell
  return MOOP_UNKNOWN;
}

int MoopMapPic::getProminentPixelColor(TVector2 position, float size)
{
  //store in a map of MOOP_* and pixel counts
  std::map<int, int> counts;
  size_t halfEdgeSizePixels = static_cast<size_t>(size * m_pos_scale / 2);
  std::pair<size_t, size_t> pixelPosition = ToPicCoordinates(position);
  
  //raster the square
  for (int x = pixelPosition.first - halfEdgeSizePixels; x < pixelPosition.first + halfEdgeSizePixels; ++x) {
      for (int y = pixelPosition.second - halfEdgeSizePixels; y < pixelPosition.second + halfEdgeSizePixels; ++y) {
	  counts[getSinglePixelColor(x, y)]++;
    } //raster on y
  } //raset on x

  //now return most prominent one, ignoring UNKNOWN
  // i.e. will return UNKNOWN only if no other pixel exist
  int maxValue = 0;
  int maxKey = MOOP_UNKNOWN;
  for (auto c : counts) {
    if (c.first == MOOP_UNKNOWN) continue; //ignore unknown
    if (c.second > maxValue) {
      maxValue = c.second;
      maxKey = c.first;
    }
  }

  return maxKey;
}

bool MoopMapPic::getXMLProperty(xmlNodePtr node, std::string name, float &x)
{
  xmlChar *prop;
  prop = xmlGetProp(node, (xmlChar*)name.c_str());
  if (!prop) {
    std::cerr << "Malformed property " << name << std::endl;
    return false;
  }
  try {
    x = atof(reinterpret_cast<char*>(prop));  
  } catch (...) {
    std::cerr << "Error in converting property " << name << " to float: " << x << std::endl;
    xmlFree(prop);
    return false;
  }
  xmlFree(prop);
  return true;
}

bool MoopMapPic::getXMLProperty(xmlNodePtr node, std::string name, int &x)
{
  xmlChar *prop;
  prop = xmlGetProp(node, (xmlChar*)name.c_str());
  if (!prop) {
    std::cerr << "Malformed property " << name << std::endl;
    return false;
  }
  try {
    x = atoi(reinterpret_cast<char*>(prop));  
  } catch (...) {
    std::cerr << "Error in converting property " << name << " to int: " << x << std::endl;
    xmlFree(prop);
    return false;
  }
  xmlFree(prop);
  return true;
}

bool MoopMapPic::getXMLProperty(xmlNodePtr node, std::string name, std::string &x)
{
  xmlChar *prop;
  prop = xmlGetProp(node, (xmlChar*)name.c_str());
  if (!prop) {
    std::cerr << "Malformed property " << name << std::endl;
    return false;
  }
  try {
    x = reinterpret_cast<char*>(prop);  
  } catch (...) {
    std::cerr << "Error in converting property " << name << " to string: " << x << std::endl;
    xmlFree(prop);
    return false;
  }
  xmlFree(prop);
  return true;
}
