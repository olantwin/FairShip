// -------------------------------------------------------------------------
// -----             Pythia6Generator source file                      -----
// -----          Created 08/08/08  by S. Spataro                      -----
// -------------------------------------------------------------------------
#include "Pythia6Generator.h"

#include "FairPrimaryGenerator.h"

#include <iostream>
#include <cstdio>

using std::cout;
using std::endl;
using std::max;

// -----   Default constructor   ------------------------------------------
Pythia6Generator::Pythia6Generator() {}
// ------------------------------------------------------------------------



// -----   Standard constructor   -----------------------------------------
Pythia6Generator::Pythia6Generator(const char* fileName) {
  fFileName  = fileName;
  fVerbose = 0;
  cout << "-I Pythia6Generator: Opening input file " << fileName << endl;
  if ((fInputFile = fopen(fFileName,"r"))==NULL)
    //  fInputFile = new ifstream(fFileName);
    //  if ( ! fInputFile->is_open() )
    Fatal("Pythia6Generator","Cannot open input file.");

  // fPDG=TDatabasePDG::Instance();
}
// ------------------------------------------------------------------------



// -----   Destructor   ---------------------------------------------------
Pythia6Generator::~Pythia6Generator() {
  CloseInput();
}
// ------------------------------------------------------------------------



// -----   Public method ReadEvent   --------------------------------------
Bool_t Pythia6Generator::ReadEvent(FairPrimaryGenerator* primGen) {

  // Check for input file
  if (!fInputFile) {
    // if ( ! fInputFile->is_open() ) {
    cout << "-E Pythia6Generator: Input file not open!" << endl;
    return kFALSE;
  }

  // Define event variable to be read from file
   Int_t ntracks = 0, eventID = 0, ncols = 0;

  // Define track variables to be read from file
  Int_t nLev = 0, pdgID = 0, nM1 = -1, nM2 = -1, nDF = -1, nDL = -1;
  Float_t fPx = 0., fPy = 0., fPz = 0., fM = 0., fE = 0.;
  Float_t fVx = 0., fVy = 0., fVz = 0., fT = 0.;

  // Read event header line from input file

  Int_t max_nr = 0;

  Text_t buffer[200];
  ncols = fscanf(fInputFile,"%d\t%d", &eventID, &ntracks);

  if (ncols && ntracks>0) {

    if (fVerbose>0) cout << "Event number: " << eventID << "\tNtracks: " << ntracks << endl;

    for (Int_t ll=0; ll<ntracks; ll++)
      {
	ncols = fscanf(fInputFile,"%d %d %d %d %d %d %f %f %f %f %f %f %f %f %f", &nLev, &pdgID, &nM1, &nM2, &nDF, &nDL, &fPx, &fPy, &fPz, &fE, &fM, &fVx, &fVy, &fVz, &fT);
	if (fVerbose>0) cout << nLev << "\t" << pdgID << "\t" << nM1 << "\t" << nM2 << "\t" << nDF << "\t" << nDL <<
	  "\t" << fPx << "\t" << fPy << "\t" << fPz << "\t" << fE << "\t" << fM << "\t" << fVx << "\t" << fVy << "\t" << fVz << "\t" << fT <<  endl;
	if (nLev==1)
	  primGen->AddTrack(pdgID, fPx, fPy, fPz, fVx, fVy, fVz);
      }
  }
  else {
    cout << "-I Pythia6Generator: End of input file reached " << endl;
    CloseInput();
    return kFALSE;
  }


  // If end of input file is reached : close it and abort run
  if ( feof(fInputFile) ) {
    cout << "-I Pythia6Generator: End of input file reached " << endl;
    CloseInput();
    return kFALSE;
  }

  /*
    cout << "-I Pythia6Generator: Event " << eventID << ",  vertex = ("
    << vx << "," << vy << "," << vz << ") cm,  multiplicity "
    << ntracks << endl;
  */

  return kTRUE;
}
// ------------------------------------------------------------------------



// -----   Private method CloseInput   ------------------------------------
void Pythia6Generator::CloseInput() {
  if ( fInputFile ) {
    //if ( fInputFile->is_open() ) {
    {
      cout << "-I Pythia6Generator: Closing input file "
	   << fFileName << endl;
      //  fInputFile->close();

      fclose(fInputFile);
    }
    delete fInputFile;
    fInputFile = NULL;
  }
}
// ------------------------------------------------------------------------
