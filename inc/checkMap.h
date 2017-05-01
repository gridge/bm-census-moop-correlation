/* Dictionary header for executable that checks digitized map content
 * Created separated header to allow ROOT dictionary creation
 * Only things that are needed to be ROOT visible are declared here.
 * and are re-declared in checkMap.cxx as external
 * 
 * Author: gridge
 */

#ifndef __CHECKMAP_H__
#define __CHECKMAP_H__

#include "TCanvas.h"
#include "TImage.h"
#include "TPaveText.h"
#include "TBox.h"
#include "MoopMapPic.h"

float mainWndWidth; //!< main window default width
float zoomRatio; //!< main window zoom ratio for initial image (for coord conversion)
TCanvas* mainWnd; //!< main frame for display, could also consider TGMainFrame
TImage* imageMap;  //!< image map, could also conside TGImageMap
TPaveText *text; //!< Text box with coordinates
TBox *moopColor; //!< rectangle showing selected MOOP color
void onClickUpdateText(); //!< Signal slot for interactive map
MoopMapPic *moopMapPic(0); //!< Moop map

#endif //__CHECKMAP_H__


