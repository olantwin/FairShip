// $Id: basiclibs.C,v 1.1.1.1 2005/06/23 07:14:09 dbertini Exp $
//
// Macro for loading basic libraries used with both Geant3 and Geant4

void basiclibs()
{
 /*
  gSystem->Load("libPluto");
  */
  gSystem->Load("libRIO");
  gSystem->Load("libGeom");
  gSystem->Load("libGeomPainter");
  gSystem->Load("libVMC");
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libPythia6");
  gSystem->Load("libPhysics");
  gSystem->Load("libNet");
  gSystem->Load("libTree");
  gSystem->Load("libMinuit");
  gSystem->Load("libMathMore");
  gSystem->Load("libpythia8");
  gSystem->Load("libgenfit.so");
  gSystem->Load("libLHAPDF.so");
}
