# Create a library called "libShipGen" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/shipgen
${CMAKE_SOURCE_DIR}/generators
${CMAKE_SOURCE_DIR}/shipdata
${CMAKE_SOURCE_DIR}/veto
${genfit2_INCDIR}
${TPYTHIA6_INCLUDE_DIR}
)

set(SYSTEM_INCLUDE_DIRECTORIES
${SYSTEM_INCLUDE_DIRECTORIES}
${PYTHIA8_INCLUDE_DIR}
${PYTHIA8_INCLUDE_DIR}/..
${HEPMC_INCLUDE_DIR}
${EVTGEN_INCLUDE_DIR}
${EVTGEN_INCLUDE_DIR}/..
)

include_directories(${INCLUDE_DIRECTORIES} ${VMC_INCLUDE_DIRS} ${FAIRROOT_INCLUDE_DIR})
include_directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
${PYTHIA8_LIBRARY_DIR}
${EVTGEN_LIBRARY_DIR}
${genfit2_LIBDIR}
${SIMPATH}/lib
)

link_directories(${LINK_DIRECTORIES})

set(SRCS
Pythia8Generator.cxx
Pythia6Generator.cxx
tPythia6Generator.cxx
HNLPythia8Generator.cxx
DPPythia8Generator.cxx
GenieGenerator.cxx
NtupleGenerator.cxx
MuonBackGenerator.cxx
CosmicsGenerator.cxx
MuDISGenerator.cxx
FixedTargetGenerator.cxx
EvtCalcGenerator.cxx
)

set(LINKDEF GenLinkDef.h)
set(LIBRARY_NAME ShipGen)

set(DEPENDENCIES Base ShipData EvtGen EvtGenExternal FairLogger::FairLogger)

GENERATE_LIBRARY()
