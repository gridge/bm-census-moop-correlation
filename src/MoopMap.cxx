/* Class to interface with moop map data.
 *
 * Author: gridge
 */

#include "MoopMap.h"

#include <algorithm>

#include "TMath.h"
#include "TString.h"

//uncomment for DEBUG of this specific class
//#define DEBUG_BUILD
#include "Utilities.h"

/**************************
 * Main functions implementation in base class
 **************************/

//ROOT dictionary
ClassImp(MoopMap)

// Declare constants
const int MoopMap::MOOP_NOTVALID=0; //!< Not valid
const int MoopMap::MOOP_GREEN=1; //!< Green
const int MoopMap::MOOP_YELLOW=2; //!< Yellow
const int MoopMap::MOOP_RED=3; //!< Red

const int MoopMap::ERROR_COORDINATE_NOT_VALID = 1;
const int MoopMap::ERROR_STREET_UNKNOWN = 2;
const int MoopMap::ERROR_LOAD_MAP = 3;
const int MoopMap::ERROR_INVALID_ASSUMPTION = 4;

const std::vector<std::string> MoopMap::MOOP_STR = {"NOT_VALID", "GREEN", "YELLOW", "RED"};

const std::vector<std::string> MoopMap::roundStreetNames = {"Esp", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "Invalid"};

const float MoopMap::dPhi_hour = 2*TMath::Pi() / 12; //15min granularity

MoopMap::MoopMap() 
{
  //init streets data (default values)
  initStreetData();

  //disable fiducial region
  setFiducialRadii(-1, 6000);
  setFiducialPhi(-1, -1); //disabled
  setFiducialExclusionCircle(0, 0, 0); //disabled
}

MoopMap::~MoopMap()
{
}

/**************************
 * Map fiduciality
 **************************/

void MoopMap::setFiducialRadii(float radiusMin, float radiusMax)
{
  m_fidMaxRadius = radiusMax;
  m_fidMinRadius = radiusMin;
}

void MoopMap::setFiducialPhi(float phiMin, float phiMax)
{
  m_fidMinPhi = phiMin;
  m_fidMaxPhi = phiMax;
  m_fidMaxPhiPrime = TVector2::Phi_0_2pi(m_fidMaxPhi-m_fidMinPhi);
}

void MoopMap::setFiducialExclusionCircle(float radiusCenter, float phiCenter, float radiusCircle)
{
  m_fidExclCirclePosition.SetMagPhi(radiusCenter, phiCenter);
  m_fidExclCircleRadius = radiusCircle;
}

bool MoopMap::isFiducial(float radius, float phi)
{
  if (radius <= m_fidMinRadius) return false;
  if (radius >= m_fidMaxRadius) return false;
  //for phi, we need to assume clock-wise!
  //Use convenience axis rotated so that x-axis match phiMin
  if (m_fidMinPhi >= 0 and m_fidMaxPhi > 0) {
    float phiPrime = TVector2::Phi_0_2pi(phi-m_fidMinPhi);
    if (phiPrime >= m_fidMaxPhiPrime) return false;
  }
  //finally check exclusion circle
  if (m_fidExclCircleRadius > 0) {
    TVector2 pos; 
    pos.SetMagPhi(radius, phi);
    if ((m_fidExclCirclePosition - pos).Mod() < m_fidExclCircleRadius) 
      return false; //inside exclusion circle
  }
  //all checks passed
  return true;
}

bool MoopMap::isFiducial(TVector2 position)
{
  return isFiducial(position.Mod(), position.Phi());
}

void MoopMap::getFiducialRegion(float& radiusMin, float& radiusMax, float& phiMin, float& phiMax)
{
  radiusMin = m_fidMinRadius;
  radiusMax = m_fidMaxRadius;
  phiMin = m_fidMinPhi;
  phiMax = m_fidMaxPhi;
}

bool MoopMap::isCampingAllowed(TVector2 position) 
{
  //check simply if a valid MOOP value is found,
  //and is within the fiducial position
  // if so then assume this is an allowed camping position
  if (not isFiducial(position)) return false;
  int moopVal = getValue(position);
  return isCampingAllowed(moopVal);
}

bool MoopMap::isCampingAllowed(int moopVal)
{
  if ((moopVal == MOOP_GREEN) or 
      (moopVal == MOOP_YELLOW) or
      (moopVal == MOOP_RED)) return true;
  return false;
}

std::map<int, float> MoopMap::getMoopArea(float epsilon)
{
  //scan map and gather information on area for each category
  //Only considers fiducial area, if specified
  if (not m_moopAreaCached.empty()) return m_moopAreaCached;  
  m_moopAreaCached[MOOP_NOTVALID] = 0.0;
  m_moopAreaCached[MOOP_GREEN] = 0.0;
  m_moopAreaCached[MOOP_YELLOW] = 0.0;
  m_moopAreaCached[MOOP_RED] = 0.0;
  float start_y = - m_fidMaxRadius + epsilon/2;
  float start_x = -m_fidMaxRadius + epsilon/2;
  for (float y= start_y; y < m_fidMaxRadius; y += epsilon) {
    for (float x = start_x; x < m_fidMaxRadius; x += epsilon) {      
      TVector2 pos(x,y);
      if (not isFiducial(pos)) continue;
      m_moopAreaCached[getValue(pos)] += epsilon*epsilon;
    }
  }
  return m_moopAreaCached;
}

std::map<int, float> MoopMap::getMoopAreaNearIntersection(streetIntersection_t intersection, float epsilon)
{
  auto result = m_cacheAreasNearIntersection.find(intersection);
  if (result == m_cacheAreasNearIntersection.end()) {
    if (isIntersectionForbidden(intersection)) {
      //intersection not valid
      return std::map<int, float>();
    }
    DBG("Invalid cache. Recalculating areas for each intersection (may take a minute)");
    //cache not valud, re-calculate
    float start_y = - m_fidMaxRadius + epsilon/2;
    float start_x = -m_fidMaxRadius + epsilon/2;
    for (float y= start_y; y < m_fidMaxRadius; y += epsilon) {
      for (float x = start_x; x < m_fidMaxRadius; x += epsilon) {      
	TVector2 pos(x,y);
	if (not isFiducial(pos)) continue;
	streetIntersection_t currentInt;
	currentInt = getClosestIntersection(pos);
	int moopVal = getValue(pos);
	if (m_cacheAreasNearIntersection.find(currentInt) == m_cacheAreasNearIntersection.end()) {
	  //initialize area
	  m_cacheAreasNearIntersection[currentInt][MOOP_NOTVALID] = 0.0;
	  m_cacheAreasNearIntersection[currentInt][MOOP_GREEN] = 0.0;
	  m_cacheAreasNearIntersection[currentInt][MOOP_YELLOW] = 0.0;
	  m_cacheAreasNearIntersection[currentInt][MOOP_RED] = 0.0;    
	}
	m_cacheAreasNearIntersection[currentInt][moopVal] += epsilon*epsilon;
      }
    }    
    //retrieve result needed by this query
    result = m_cacheAreasNearIntersection.find(intersection);
  }
  return result->second;
}

std::map<int, float> MoopMap::getMoopArea(TVector2 pos, float radius, float epsilon)
{
  //scan map around the given position
  //Only considers fiducial area, if specified
  std::map<int, float> areas;
  areas[MOOP_NOTVALID] = 0.0;
  areas[MOOP_GREEN] = 0.0;
  areas[MOOP_YELLOW] = 0.0;
  areas[MOOP_RED] = 0.0;
  float start_y = pos.Y() - radius;
  float start_x = pos.X() - radius;;
  for (float y= start_y; y <= pos.Y() + radius; y += epsilon) {
    for (float x = start_x; x <= pos.X() + radius; x += epsilon) {      
      TVector2 currentPos(x,y);
      if (not isFiducial(currentPos)) continue;
      areas[getValue(currentPos)] += epsilon*epsilon;
    }
  }
  return areas;
}


/**************************
 * Functions for Roads
 **************************/

float MoopMap::getPhiFromHours(int hours, int minutes) 
{
  if ((hours < 0) or (hours > 12)) {
    throw ERROR_COORDINATE_NOT_VALID;
  }
  if ((minutes < 0) or (minutes >= 60)) {
    throw ERROR_COORDINATE_NOT_VALID;
  }
  //clock-based phi
  float phi_prime = (hours * 2*TMath::Pi()/12 + minutes * 2*TMath::Pi()/12/60);
  //transform in current coordinates (phi=0 @ 3:00)
  return TVector2::Phi_0_2pi(phi_prime - TMath::Pi()/2);
}

void MoopMap::getHoursFromPhi(float phi, int& hours, int& minutes)
{
  phi = TVector2::Phi_0_2pi(phi);
  hours = static_cast<int>(phi / dPhi_hour);
  minutes = 60 * (phi / dPhi_hour - hours);
  //adjust for ration of coordinate system
  hours = hours + 3;
  if (hours > 12) hours = hours - 12;
}

MoopMap::radialStreet_t MoopMap::getClosestRadialStreet(int hours, int minutes) 
{
  //index is simply a 15-minutes granularity
  //cast to float to allow rounding, then cast back to int
  return 4*hours + round((float)minutes / 15.0);
}

MoopMap::radialStreet_t MoopMap::getClosestRadialStreet(float phi)
{
  static float dPhi_15mins = 2*TMath::Pi() / 12 / 4;
  radialStreet_t tmpStreet = round(phi / dPhi_15mins);
  //now need to add offset (phi=0 -> index = 3*4)
  tmpStreet = tmpStreet + 3*4; //3:00 * 4 (15mins intervals)
  if (tmpStreet > 12*4) {
    tmpStreet = tmpStreet - 12*4;
  }
  return tmpStreet;
}

void MoopMap::getHoursFromRadialStreet(radialStreet_t index, int& hours, int &minutes)
{
  hours = index / 4;
  minutes = (index % 4)*15;  
}

float MoopMap::getPhiFromRadialStreet(radialStreet_t index)
{
  int hours,mins;
  getHoursFromRadialStreet(index, hours, mins);
  return getPhiFromHours(hours, mins);
}

TVector2 MoopMap::getPointAlongStreet(roundStreet_t street, int hours, int minutes) 
{
  float radius = m_streetRadius[street];
  float phi = getPhiFromHours(hours, minutes);
  TVector2 v;
  v.SetMagPhi(radius, phi);

  return v;
}

MoopMap::roundStreet_t MoopMap::getNextRoundStreet(roundStreet_t street, bool outwardDirection)
{
  int idx = static_cast<int>(street);
  if (outwardDirection) idx = idx + 1;
  else idx = idx - 1;
  if ((idx >= static_cast<int>(street_NRoundStreets)) or (idx < 0)) {
    //outside boundaries
    return street_NRoundStreets;
  }
  return static_cast<roundStreet_t>(idx);
}

MoopMap::streetIntersection_t MoopMap::getIntersection(roundStreet_t streetRound, radialStreet_t streetRadial)
{
  return std::make_pair(streetRound, streetRadial);
}

TVector2 MoopMap::getIntersectionPosition(roundStreet_t streetRound, radialStreet_t streetRadial) 
{
  return getIntersectionPosition(getIntersection(streetRound, streetRadial));
}

TVector2 MoopMap::getIntersectionPosition(streetIntersection_t intersection)
{
  int hours, mins;
  getHoursFromRadialStreet(intersection.second, hours, mins);
  return getPointAlongStreet(intersection.first, hours, mins);
}

void MoopMap::addForbiddenIntersection(streetIntersection_t intersection)
{
  if (not isIntersectionForbidden(intersection)) //check if not already forbidden
    m_forbiddenIntersections.push_back(intersection);
}

void MoopMap::clearForbiddenIntersections()
{
  m_forbiddenIntersections.clear();
}

bool MoopMap::isIntersectionForbidden(streetIntersection_t intersection)
{
  return (std::find(m_forbiddenIntersections.begin(), m_forbiddenIntersections.end(), intersection) != m_forbiddenIntersections.end());
}


MoopMap::streetIntersection_t MoopMap::getClosestIntersection(TVector2 pos, bool mustBeValid)
{
  //Start finding the two closest round streets
  // Takes advantage of radial ordering of streets to make fast search
  //and stores at the same time the two closest points
  int closestIdx=-1;
  int secondClosestIdx=-1;
  int idx_min = 0;
  int idx_max = street_NRoundStreets;
  int idx_probe = idx_max / 2;
  
  do {
    float dist = pos.Mod() - m_streetRadius[idx_probe];
    if (dist > 0) {
      //we need to go to larger radii
      idx_min = idx_probe;
      idx_probe = (idx_max-idx_min)/2 + idx_min;
    } else if (dist < 0) {
      //we need to go to smaller radii
      idx_max = idx_probe;
      idx_probe = (idx_max-idx_min)/2 + idx_min;      
    } else {
      //rare occasion! The two positions exactly equal
      closestIdx = idx_probe; 
      //if closest is the last element, we know the answer
      if (closestIdx == street_NRoundStreets-1) secondClosestIdx = street_NRoundStreets-2;
      //if closest is the first element, we know the answer
      else if (closestIdx == 0) secondClosestIdx = 1;
      else {
	float dist_down = TMath::Abs(pos.Mod() - m_streetRadius[closestIdx-1]);
	float dist_up = TMath::Abs(pos.Mod() - m_streetRadius[closestIdx+1]);
	if (dist_down < dist_up) secondClosestIdx = closestIdx - 1;
	else secondClosestIdx = closestIdx + 1;
      }
      break;
    }
  } while (not ((idx_probe == idx_max) or (idx_probe == idx_min)));
  if (closestIdx < 0) {
    float min_1 = TMath::Abs(pos.Mod() - m_streetRadius[idx_min]);
    float min_2 = TMath::Abs(pos.Mod() - m_streetRadius[idx_max]);
    if (min_1 < min_2) {
      closestIdx = idx_min;
      secondClosestIdx = idx_max;
    } else {
      closestIdx = idx_max;
      secondClosestIdx = idx_min;
    }    
   
  }
  roundStreet_t closestRoundSt = static_cast<roundStreet_t>(closestIdx);
  roundStreet_t secondClosestRoundSt = static_cast<roundStreet_t>(secondClosestIdx);
  DBG("Closest and second closest round streets in radius are: " << closestRoundSt << ", " << secondClosestRoundSt);

  //Now find the two closest radial streets
  // Takes advantage of the simple relantionship between phi and radial streets
  radialStreet_t closestRadialSt = getClosestRadialStreet(pos.Phi());
  radialStreet_t secondClosestRadialSt;      
  static float dPhi_15mins = 2*TMath::Pi() / 12 / 4;
  static int maxRadialStreets = 12*4;
  float phi_closest = getPhiFromRadialStreet(closestRadialSt);
  float dphi = pos.Phi() - phi_closest;
  if (dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  if (dphi > 0) {
    if (closestRadialSt < maxRadialStreets-1)
      secondClosestRadialSt = closestRadialSt + 1;
    else secondClosestRadialSt = 0;
  }    else {
    //dphi < 0
    if (closestRadialSt > 0)
      secondClosestRadialSt = closestRadialSt - 1;
    else secondClosestRadialSt = maxRadialStreets-1;
  }
  DBG("Closest and second closest radial street indeces: " << closestRadialSt << ", " << secondClosestRadialSt);

  //now make a decision
  streetIntersection_t closestInt = getIntersection(closestRoundSt, closestRadialSt);
  if ((not mustBeValid) or (not isIntersectionForbidden(closestInt))) {
    DBG("Closest point allowed:" << getStrIntersection(closestInt));
    return closestInt;
  }
  //otherwise we need to figure out which of the other three intersections is closer
  //We assume that the forbidden intersections are never two-in-a-row (true until now)
  // which ensure thsi is enough to always be correct.
  std::vector<float> dist = {-1, -1, -1};
  std::vector<streetIntersection_t> points;
  points.reserve(3);
  points[0]=std::make_pair(closestRoundSt, secondClosestRadialSt);
  if (not isIntersectionForbidden(points[0]))
    dist[0]=(getIntersectionPosition(points[0]) - pos).Mod();
  points[1]=std::make_pair(secondClosestRoundSt, closestRadialSt);
  if (not isIntersectionForbidden(points[1]))
    dist[1]=(getIntersectionPosition(points[1]) - pos).Mod();
  points[2]=std::make_pair(secondClosestRoundSt, secondClosestRadialSt);
  if (not isIntersectionForbidden(points[2]))
    dist[2]=(getIntersectionPosition(points[2]) - pos).Mod();
  int idxMin=-1;
  float distMin=-1;
  DBG("Closest point forbidden. Checking neighboroods:");
  for (int i=0; i < 3; ++i) {
    DBG("Point " << i << ", distance" << dist[i]);
    if (dist[i] >= 0) {
      if ((distMin < 0) or (dist[i] < distMin)) {
	distMin = dist[i];
	idxMin = i;
      }	
    }
  }
  if (idxMin < 0) {
    DBG("Invalida assumptionexception in calculating closest intersection for (r,phi)=(" << pos.Mod() << ", " << pos.Phi() << ")");
    //return invalid intersection
    // This may happen for instance in some areas near the central camp
    // @TODO: could instead use a sllower and more robust algorithm, but there's no need right now
    return getIntersection(street_NRoundStreets, 0);
  }
  return points[idxMin];

}

std::string MoopMap::getStrIntersection(streetIntersection_t intersection)
{
  TString retInt;
  int hours,mins;
  getHoursFromRadialStreet(intersection.second, hours, mins);
  
  retInt.Form("%s @ %02d:%02d", roundStreetNames[(unsigned int)intersection.first].c_str(), hours, mins);
  
  return std::string(retInt.Data());
}

MoopMap::roundStreet_t MoopMap::getRoundStreet(std::string name)
{
  for (unsigned int idx=street_Esplanade; idx < street_NRoundStreets; idx++) {
    if (roundStreetNames[idx] == name) {
      //found it
      return static_cast<roundStreet_t>(idx);
    }
  }
  //check for known variants
  if (name == "Esplanade") {
    return street_Esplanade;
  }
  return street_NRoundStreets; //not found
}

std::vector<MoopMap::streetIntersection_t> MoopMap::getListOfIntersections(bool onlyValidIntersections)
{
  std::vector<MoopMap::streetIntersection_t> listIntersections;
  //loop over round streets
  for (unsigned int idxRStreet=street_Esplanade; idxRStreet < street_NRoundStreets; ++idxRStreet) {
    roundStreet_t iRound = static_cast<roundStreet_t>(idxRStreet);
    //now loop over radial streets (every 15 min.)
    for (float idxPhi=0; idxPhi < 12*4; ++idxPhi) {
      radialStreet_t iRadial = static_cast<radialStreet_t>(idxPhi);
      streetIntersection_t intersection = getIntersection(iRound, iRadial);
      //check if intersection is fiducial
      TVector2 pos = getIntersectionPosition(intersection);
      if (not isFiducial(pos)) continue;
      //check if intersection is forbidden
      if (isIntersectionForbidden(intersection)) continue;
      //add to the list
      listIntersections.push_back(intersection);
    }
  }
  return listIntersections;
}

/**************************
 * Private members
 **************************/

//Init internal data
void MoopMap::initStreetData() 
{
  //Hard-coded radii of known streets
  //Will be overridden by XML loaded map
  m_streetRadius.reserve(street_NRoundStreets);
  m_streetRadius[street_Esplanade] = 2510.0;
  m_streetRadius[street_A] = 2950.0;
  m_streetRadius[street_B] = 3190.0;
  m_streetRadius[street_C] = 3430.0;
  m_streetRadius[street_D] = 3670.0;
  m_streetRadius[street_E] = 3910.0;
  m_streetRadius[street_F] = 4150.0;
  m_streetRadius[street_G] = 4390.0;
  m_streetRadius[street_H] = 4630.0;
  m_streetRadius[street_I] = 4870.0;
  m_streetRadius[street_J] = 5110.0;
  m_streetRadius[street_K] = 5350.0;
  m_streetRadius[street_L] = 5590.0;
}

