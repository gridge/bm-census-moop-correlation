/* Class to interface with moop map data.
 * 
 * Author: gridge
 * 
 */

#ifndef _MOOPMAP_H_
#define _MOOPMAP_H_

#include <string>
#include <map>

#include "TVector2.h"

#include "Utilities.h"

/* Class to interface with moop map data.
 * Position within a map is stored in a TVector2
 * The main coordinate system used for the MOOP map is polar
 * r = radius [ft]
 * phi = polar angle [rad]
 * The convention is that the phi = 0 at 3:00 hrs and it is measured clockwise in radians. Domain: [0,2*pi)
 * The radius is consistent with the MOOP map reference system, with zero being at the Man and it is measured in feets.
 * Cartesian coordinates can also be used (units: [ft]), with x-axis pointing left-to-right in the MOOP image and y-axis
 * pointing up-to-bottom in the MOOP map images.
 * Recognized street names are identified by:
 *  - Esplanade
 *  - A to L
 */
class MoopMap {
 public:

  /* @defgroup Virtual Pure virtual functions that need implementation in a derived class.
   * @{ */

  /// Get MOOP value for a given point
  virtual int getValue(TVector2 position) = 0;

  /* @} */

  /** Constructor */
  MoopMap();

  ///Default desctructor
  ~MoopMap();  

  /// Check if given point allows camping, including fiducial position
  virtual bool isCampingAllowed(TVector2 position);

  /// Check if given point allows camping (no fiducial position check possible!)
  virtual bool isCampingAllowed(int moopVal);

  /** \defgroup Roads Manage roads and intersections.
   * Streets at constant "radius" are named "round streets"
   * they are defined in an enum from the lowest radius and increasing.
   * Street at constant "phi" are named "radial streets" and called by the hours. 
   * The granularity of these streets is fixed at 15min.
   * What is most important are the valid intersections. By default every 
   * intersection between the two type of roads is considered, however a set of 
   * intersections can be excluded from a given list to re-create the real situation
   * observed in the map.
   * @{
   */

  ///Enum defining valid round street names (at constant phi). Needs to start with 0 (zero)
  enum roundStreet_t {
    street_Esplanade=0,
    street_A,
    street_B,
    street_C,
    street_D,
    street_E,
    street_F,
    street_G,
    street_H,
    street_I,
    street_J,
    street_K,
    street_L,
    street_NRoundStreets
  };

  /// Vector of street names with the same order as in roundStreet_t
  static const std::vector<std::string> roundStreetNames;

  /** Hours:minutes representation of streets (at constant radius).      
   * The integer value is simply translated into a 15minutes granularity
   * counter. E.g. 2:30 -> 2*4 + 30/15 = 8 + 2 = 10
   * Use the utility function getStreetHourIdx(int hours, int minutes)
   * and getHoursFromIndex(int index, int& hours, int& minutes)
   * to easily convert.
   */
  typedef int radialStreet_t;

  /// Type for intersection specification
  typedef std::pair<roundStreet_t, radialStreet_t> streetIntersection_t;

  /// Convert hours to radiants
  static float getPhiFromHours(int hours, int minutes=0);

  /// Convert phi into hours and minutes
  static void getHoursFromPhi(float phi, int& hours, int& minutes);

  /** Convert hours and minutes to closest radial street index
   * @param hours
   * @param minutes
   * @return street index corresponding to the closest 15min interval to hours:minutes 
   */  
  static radialStreet_t getClosestRadialStreet(int hours, int minutes);

  /** Convert phi to closest radial street index
   * @param phi angle in MOOP coordinates
   * @return street index corresponding to the closest 15min interval to hours:minutes 
   */  
  static radialStreet_t getClosestRadialStreet(float phi);

  /** Convert radial street index back to hours and minutes
   * @param index hours index
   * @param hours returns hours correposnding to index
   * @param minutes returns minutes correposnding to index
   */
  static void getHoursFromRadialStreet(radialStreet_t index, int& hours, int &minutes);

  ///Get phi corresponding to a given radial sreet
  static float getPhiFromRadialStreet(radialStreet_t index);

  /// Get coordinates of a given point on map on a street Name
  virtual TVector2 getPointAlongStreet(roundStreet_t street, int hours, int minutes=0);

  /** Get next/previous round street in increasing radius
   * @param street input round street
   * @param outwardDirection true is increasing radius, false is decreasing radius
   * @return next or previous radial street, street_NRoundStreets if outisde bounds
   */
  virtual roundStreet_t getNextRoundStreet(roundStreet_t street, bool outwardDirection=true);

  /// Identify intersection (just builds pair)
  virtual streetIntersection_t getIntersection(roundStreet_t streetRound, radialStreet_t streetRadial);

  /// Get coordinates of a given street intersection
  virtual TVector2 getIntersectionPosition(roundStreet_t streetRound, radialStreet_t streetRadial);

  /// Get coordinates of a given street intersection

  virtual TVector2 getIntersectionPosition(streetIntersection_t intersection);

  /// Add to the list of forbidden intersections
  virtual void addForbiddenIntersection(streetIntersection_t intersection);

  /// Clear list of forbidden intersections
  virtual void clearForbiddenIntersections();

  /// Check if intersection is forbidden
  virtual bool isIntersectionForbidden(streetIntersection_t intersection);

  ///Retrieve list of forbidden intersections

  /** Find closest valid intersection for a given point
   * Note: this is optimized for the current setup, assumes that taking the closest two streets of each type will
   *       always have a valid intersections
   * @param pos position in MOOP coordinates
   * @param mustBeValid requires the intersection to be marked as valid (not forbidden)
   * @return closest intersection found
   */
  virtual streetIntersection_t getClosestIntersection(TVector2 pos, bool mustBeValid=true);

  /// Convert intersection in a compact string
  virtual std::string getStrIntersection(streetIntersection_t intersection);

  /// Get Round street from string its name
  virtual roundStreet_t getRoundStreet(std::string name);

  /// Get list of intersections inside fiducial region
  virtual std::vector<streetIntersection_t> getListOfIntersections(bool onlyValidIntersections=true);

  /** @} */

  /** Get total MOOP area by category
   * Units are ft^2.
   * @param epsilon step for integral evaluation (in ft). Increase to faster results.
   * @return a map of moop category and area
   */
  virtual std::map<int, float> getMoopArea(float epsilon=1.0);

  /** Get MOOP area by category around a given position
   * Units are ft^2. Only fiducial areas are considered.
   * @param pos center position
   * @param radius edge of square to consider around the position
   * @param epsilon step for integral evaluation (in ft). Increase to faster results.
   * @return a map of moop category and area
   */
  std::map<int, float> getMoopArea(TVector2 pos, float radius, float epsilon=1.0);

  /** Get MOOP areas for points that are closest to the given intersection
   * Units are ft^2. Only valid intersections are considered.
   * Results are evaluated at first request and cached.
   * @param intersection intersection to consider
   * @param epsilon step for integral evaluation (in ft). Increase to faster results.
   * @return a map of moop category and area
   */
  virtual std::map<int, float> getMoopAreaNearIntersection(streetIntersection_t intersection, float epsilon=1.0);



  // Set simple fiducial regions
  /// Fiducial radius. 
  virtual void setFiducialRadii(float radiusMin, float radiusMax);
  /// Fiducial phi. IMPORTANT: Assumes clock-wise direction from phiMin to phiMax
  virtual void setFiducialPhi(float phiMin, float phiMax);
  //Set exclusion circle for e.g. center camp
  virtual void setFiducialExclusionCircle(float radiusCenter, float phiCenter, float radiusCircle);
  /// Test fiducial region
  virtual bool isFiducial(float radius, float phi);
  /// Test fiducial region (MOOP coordinates)
  virtual bool isFiducial(TVector2 position);
  /// Get fiducial region
  virtual void getFiducialRegion(float& radiusMin, float& radiusMax, float& phiMin, float& phiMax);
  
 public:  

  /* @defgroup ErrorCodes Error codes
   * @{
   */

  static const int ERROR_COORDINATE_NOT_VALID ;
  static const int ERROR_STREET_UNKNOWN;
  static const int ERROR_LOAD_MAP;
  static const int ERROR_INVALID_ASSUMPTION;

  /* @} */

  /// MOOP map values
  static const int MOOP_NOTVALID; //!< Not valid
  static const int MOOP_GREEN; //!< Green
  static const int MOOP_YELLOW; //!< Yellow
  static const int MOOP_RED; //!< Red

  static const std::vector<std::string> MOOP_STR;

 protected:  
  /* @addGroup Roads Data to identify street positions
     @{
   */

  /// Map of street radii, in increasing order
  std::vector<float> m_streetRadius;

  /* @} */

  /* @defgroup fiducialRegion Set a fiducial region on the map to limit data to consider
     @{
   */
  
  ///Simple radial fiducial region: minimu radius
  float m_fidMinRadius;
  ///Simple radial fiducial region: maximum radius
  float m_fidMaxRadius;  
  ///Simple angular fiducial region: minimum phi
  float m_fidMinPhi;  
  ///Simple angular fiducial region: maximum phi
  float m_fidMaxPhi;  
  ///Fiducial region circle exclusion (e.g. for center camp): center radius
  TVector2 m_fidExclCirclePosition;
  ///Fiducial region circle exclusion (e.g. for center camp): circle radius
  float m_fidExclCircleRadius;

  /* @} */

 private:
  /// Constant granularity of radial streets
  static const float dPhi_hour;

  ///Build internal street representation
  void initStreetData();

  /// @addGroup fiducialRegion Convenience variable to store rotated-frame for fast fiducial check
  float m_fidMaxPhiPrime;

  /// @addGroup Roads Store forbidden intersections
  std::vector<streetIntersection_t> m_forbiddenIntersections;

  /// Cache for Moop total area
  std::map<int, float> m_moopAreaCached;

  /// Cache for areas near intersection
  std::map<streetIntersection_t, std::map<int, float>> m_cacheAreasNearIntersection;

  //Root dictionary
  ClassDef(MoopMap, 1)

};

#endif

