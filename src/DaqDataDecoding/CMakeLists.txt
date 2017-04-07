PROJECT                 (DaqDataDecoding)

# Type checks
INCLUDE                 (CheckTypeSize)
CHECK_TYPE_SIZE         ("char     "    SIZEOF_CHAR)
CHECK_TYPE_SIZE         ("short    "    SIZEOF_SHORT)
CHECK_TYPE_SIZE         ("int      "    SIZEOF_INT)
CHECK_TYPE_SIZE         ("long     "    SIZEOF_LONG)
CHECK_TYPE_SIZE         ("long long"    SIZEOF_LONG_LONG)

##########################
# Packages configuration #
##########################

SET                     (RFIO_DIR ${PROJECT_SOURCE_DIR}/.cmake)
FIND_PACKAGE            (RFIO)

SET                     (DATE_DIR ${PROJECT_SOURCE_DIR}/.cmake)
FIND_PACKAGE            (DATE)


SET                     (DDD_SRC ${PROJECT_SOURCE_DIR}/src)
SET                     (DDD_EXE ${PROJECT_SOURCE_DIR}/examples/how-to)

#######################
# Read version number #
#######################

FILE                    (READ ${DDD_SRC}/version.h VERSION_H)
STRING                  (REGEX REPLACE "(.*DaqDataDecoding_VERSION +)([0-9]+).*" "\\2" VERSION_MAJOR ${VERSION_H})
STRING                  (REGEX REPLACE "(.*DaqDataDecoding_VERSION_MINOR +)([0-9]+).*" "\\2" VERSION_MINOR ${VERSION_H})
MESSAGE                 (STATUS "DaqDataDecoding version (from src/version.h): ${VERSION_MAJOR}.${VERSION_MINOR}")

#############################
# Create configuration file #
#############################

CONFIGURE_FILE          (${DDD_SRC}/config.h.cmake ${DDD_SRC}/config.h)

#######################
# Compilation/Linking #
#######################

INCLUDE_DIRECTORIES     (${DDD_SRC} ${DATE_INCLUDE_DIR})

FILE                    (GLOB DDD_SRC ${DDD_SRC}/*.cc)

ADD_LIBRARY             (DaqDataDecoding STATIC ${DDD_SRC})

SET                     (COMMON_LIBS expat dl popt DaqDataDecoding ${DATE_LIBRARIES} ${RFIO_LIBRARIES})

ADD_EXECUTABLE          (ddd ${DDD_EXE}/ddd.cc)
TARGET_LINK_LIBRARIES   (ddd ${COMMON_LIBS})

ADD_EXECUTABLE          (maps ${DDD_EXE}/maps.cc)
TARGET_LINK_LIBRARIES   (maps ${COMMON_LIBS})

############################
# Documentation generation #
############################

INCLUDE                 (FindDoxygen)

CONFIGURE_FILE          (${PROJECT_SOURCE_DIR}/Doxyfile.cmake ${PROJECT_SOURCE_DIR}/Doxyfile)
ADD_CUSTOM_TARGET       (doc
                         COMMAND            ${DOXYGEN}
                         WORKING_DIRECTORY  ${PROJECT_SOURCE_DIR})

################
# Installation #
################

FILE                    (GLOB DDD_INC ${DDD_SRC}/*.h)
INSTALL                 (FILES ${DDD_INC} DESTINATION include/DaqDataDecoding)
#INSTALL                 (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ddd ${CMAKE_CURRENT_BINARY_DIR}/maps DESTINATION bin)

INSTALL                 (TARGETS ddd maps DaqDataDecoding
                        RUNTIME DESTINATION bin
                        LIBRARY DESTINATION lib
                        ARCHIVE DESTINATION lib)

FILE                    (GLOB CMAKE_MODULES ${PROJECT_SOURCE_DIR}/.cmake/*.cmake)
INSTALL                 (FILES ${CMAKE_MODULES} DESTINATION share/CMake)

#INSTALL                 (DIRECTORY ${PROJECT_SOURCE_DIR}/doc
#                         DESTINATION doc/DaqDataDecoding)
#INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/doc DESTINATION doc)
