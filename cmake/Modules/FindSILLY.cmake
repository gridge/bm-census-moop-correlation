# Custom written to find libSILLY in the system
# author: gridge
#
#First search dependencies

#Search for lib and includes
find_path(SILLY_INCLUDE_DIR SILLY/SILLY.h)
find_library(SILLY_LIBRARY NAMES SILLY libSILLY)

#Now put everything together
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SILLY DEFAULT_MSG SILLY_LIBRARY SILLY_INCLUDE_DIR)
mark_as_advanced(SILLY_INCLUDE_DIR SILLY_LIBRARY )

set(SILLY_LIBRARIES ${SILLY_LIBRARY} )
set(SILLY_INCLUDE_DIRS ${SILLY_INCLUDE_DIR} )

#Specify to use SILLY Options file to enable inline definitions
set(SILLY_CXX_FLAGS -DUSE_SILLYOPTIONS_H)