# Create a library called "libShipData" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/shipdata
${CMAKE_SOURCE_DIR}/genfit/core/include
${CMAKE_SOURCE_DIR}/genfit/fitters/include
${VMC_INCLUDE_DIRS}
)

include_directories(${INCLUDE_DIRECTORIES} ${FAIRROOT_INCLUDE_DIR})
include_directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
${CMAKE_SOURCE_DIR}/shipdata

)

link_directories( ${LINK_DIRECTORIES})

set(SRCS
ShipHit.cxx
ShipStack.cxx
ShipMCTrack.cxx
ShipParticle.cxx
TrackInfo.cxx
)

Set(HEADERS )
Set(LINKDEF MCStackLinkDef.h)
Set(LIBRARY_NAME ShipData)
set(DEPENDENCIES Base EG Physics Core genfit FairLogger::FairLogger)

GENERATE_LIBRARY()
