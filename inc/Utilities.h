/* General utility functions
 *
 * Author: gridge
 */

#include <iostream>

///Debug macro, define DEBUG_BUILD to enable where necessary and recompile
//Clean previous define if it exists to allow per-file debug preferences
#ifdef DBG 
#undef DBG
#endif
//now define DBG if requested
#ifdef DEBUG_BUILD
#pragma GCC warning "Enabled DEBUG feature."
#define DBG(x) do { std::cerr << "DEBUG: " << x << std::endl; } while (0)
#else
#define DBG(x)
#endif

#ifndef __UTILITIES_H__
#define __UTILITIES_H__

#endif
