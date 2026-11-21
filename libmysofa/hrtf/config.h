#ifndef _CONFIG_H
#define _CONFIG_H

/* Minimal replacement for the CMake-generated config.h of the vendored
 * libmysofa (see ../README.md for the upstream revision). */
#define MYSOFA_PACKAGE_VERSION "1.3.5"
#define CPACK_PACKAGE_VERSION_MAJOR 1
#define CPACK_PACKAGE_VERSION_MINOR 3
#define CPACK_PACKAGE_VERSION_PATCH 5

/* Only used if mysofa_load() falls back to the default data set; the runtime
 * loader always supplies data from the file, so this path is unreachable. */
#define CMAKE_INSTALL_PREFIX "/usr"

#endif
