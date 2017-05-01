/* Class to interface with moop map data.
 * Derives from moopmap, and takes moop and camping data from a JPG of the moop map.
 * The JPG is edited for clarity:
 *  - only camping sites are left with a color (used to identify is a point is a camping point)
 * In order to load a picture, both input file and reference points need to be set. 
 *
 * Author: gridge
 */

#ifndef _MOOPMAPPIC_H_
#define _MOOPMAPPIC_H_

#include "MoopMap.h"

#include "SILLY/SILLY.h"
#include <libxml/parser.h>

/// Specific implementation of MoopMap with input JPG.
class MoopMapPic : public MoopMap {
 public:

  MoopMapPic();

  ~MoopMapPic();

  /// Get MOOP value for a given point
  virtual int getValue(TVector2 position);

  /** Load map from file, given two reference points to translate coordinates
   * @pic input picture (edited for clarity, leave with color only camping sites)
   * @param positionCenter (x,y) position of center (man) in the picture, in pixels
   * @param positionOtherRef (x,y) position of a second refernce point given by refPoint
   * @param refPoint local coordinates of second reference point
   * @return true on success
   */
  bool load(std::string pic, TVector2 positionCenter, TVector2 positionOtherRef, TVector2 refPoint);

  /** Load map from XML file.
   * This is the recommended approach since the XML file will contain most of the useful settings already
   * without the need to add them after manually in each executable.
   * @param XML input file.
   @ @return true on success
   */
  bool loadFromXML(std::string mapXMLFile);

  /** Set RGB tolerance.
     Tolerance is defined calculating an euclidean distance in the "RGB" space.
     @param tolerance new tolerance allowed
   */
  void setRGBTolerance(float newTolerance);

  /** Get MOOP coordinates from pixel x,y.
   * @param pixel_x position of pixel along row
   * @param pixel_y position of pixel along column
   * @return MOOP cordinates
   */
  TVector2 getMoopCoordinatesFromPixel(size_t pixel_x, size_t pixel_y);

 public:
  /** \defgroup Settings Publicly accessible members describing the loaded MAP
   * @{
   */
  
  /// Version of XML file used in loading
  std::string xmlVersion;

  /// MOOP map year
  std::string mapYear;

  /// File used to load MOOP data
  std::string inputFile;

  /// Type of input file used
  std::string inputType;  

  /** @} */

 protected:
  ///Internal storage for MOOP picture
  //Store image location and data
  SILLY::FileDataSource *m_pictureLocation;
  SILLY::Image *m_moopPicture;
  ///Size of each pixel loaded
  size_t m_sizePixel;

  ///Static member to keep track of initialization/destruction of SILLY resources
  static int nInstances;

  /** @defgroup CoorUtilities Methods to transform from/to pixel coordinates.
   * @{ */

  /** Convert MOOP coordinates into pixel (x,y) of loaded image
   * @param position position in MOOP map coordinates
   * @return position in pixels (first=x, second=y) for image manipulation
   */
  std::pair<size_t, size_t> ToPicCoordinates(TVector2 position);

  /** Convert MOOP coordinates into pixel index of loaded image
   * @position position in MOOP map coordinates
   * @return position as 1-dimensional index in image buffer
   */
  size_t ToPicBufferIdx(TVector2 position);

  /** Convert pixel position (x,y) into pixel index of loaded image
   * @param pixel_x position of pixel along row
   * @param pixel_y position of pixel along column
   * @return position in pixels for image manipulation
   */
  size_t ToPicBufferIdx(size_t pixel_x, size_t pixel_y);

  /// Reference point 1: center
  TVector2 m_pos_center;

  /** Scale factor from image pixel distance and MOOP map radial coordinate.
     Defined as distance_image_pixels / distance_internalCoords
   */
  float m_pos_scale;

  /** @} */

  /** @defgroup ColorUtilities Define color manipulation and retrival methods.
   * @{ */

  /// Tolerance for color-matching (sqrt(delta_R**2+delta_G**2+delta_B**2)
  float m_tolerance_RGB;

  /** Get distance in color-sapce.
   * Using euclidean distance. */
  float getColorDistance(std::vector<unsigned char> color1, std::vector<unsigned char> color2);

  ///Define extra MOOP type to keep track of areas where we can't tell (internal usage only)
  static const int MOOP_UNKNOWN=-1; //!< Can't tell

  //Reference RGB values for MOOP colors
  static const std::vector<unsigned char> MOOP_GREEN_RGB;
  static const std::vector<unsigned char> MOOP_YELLOW_RGB;
  static const std::vector<unsigned char> MOOP_RED_RGB;
  //Add reference values for white background (no camp) and black (writings)
  static const std::vector<unsigned char> MOOP_BLACK_RGB;
  static const std::vector<unsigned char> MOOP_WHITE_RGB;

  /** Utility function to get MOOP color of single pixel.     
   * @param position position in MOOP coordinates
   * @return one of the MOOP_* enums
   */
  int getSinglePixelColor(TVector2 position);

  /** get MOOP color of single pixel.
     Input are pixel coordinates.
     @param x pixel position in the row
     @param y pixel position in the column
     @return one of the MOOP_* enums
   */
  int getSinglePixelColor(size_t pixel_x, size_t pixel_y);

  /** Utility function to get MOOP most prominent color in a given squared area.
   * @param position position (MOOP coordinates) of center of the area
   * @param size edge size (in feets) of rectangle centered on position
   */
  int getProminentPixelColor(TVector2 position, float size);

  /** @} */

  /** Utility to get xml tag property
   * @param node XML node pointer
   * @param name property name
   * @param x property value
   * @return true if successful
   */
  bool getXMLProperty(xmlNodePtr node, std::string name, float &x);

  /** Utility to get xml tag property
   * @param node XML node pointer
   * @param name property name
   * @param x property value
   * @return true if successful
   */
  bool getXMLProperty(xmlNodePtr node, std::string name, int &x);

  /** Utility to get xml tag property
   * @param node XML node pointer
   * @param name property name
   * @param x property value
   * @return true if successful
   */
  bool getXMLProperty(xmlNodePtr node, std::string name, std::string &x);

  //Root dictionary generation
  ClassDef(MoopMapPic, 1)

};

#endif
