
// Configuration macro for Geane VirtualMC

void Config()
{

  gMC3= new  TGeant3TGeo("C++ Interface to Geant3");
  cout << "-I- G3Config: Geant3 with TGeo has been created for Geane."
       << endl;
  // create Fair Specific Stack
  FairStack *st = new FairStack(10);
  gMC3->SetStack( st ) ;

  // ******* GEANEconfiguration for simulated Runs  *******
  gMC3->SetDEBU(0, 0, 1);
  gMC3->SetSWIT(4, 10);

  gMC3->SetDCAY(0);
  gMC3->SetPAIR(0);
  gMC3->SetCOMP(0);
  gMC3->SetPHOT(0);
  gMC3->SetPFIS(0);
  gMC3->SetDRAY(0);
  gMC3->SetANNI(0);
  gMC3->SetBREM(1);
  gMC3->SetMUNU(0);
  gMC3->SetCKOV(0);
  gMC3->SetHADR(0);         //Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3)//4 fluka 5 gcalor
  gMC3->SetLOSS(4);
  gMC3->SetMULS(1); 	    //1=Moliere,3=Gaussian
  gMC3->SetRAYL(0);
  gMC3->SetSTRA(0);

  gMC3->SetAUTO(1);         //Select automatic STMIN etc... calc. (AUTO 1) or manual (AUTO 0)
  gMC3->SetABAN(0);         //Restore 3.16 behaviour for abandoned tracks
  gMC3->SetOPTI(0);         //Select optimisation level for GEANT geometry searches (0,1,2)
  gMC3->SetERAN(5.e-7);


  // -------->>>>> PAY ATTENTION!!!!!
  // For a correct use of GEANE, you MUST use the cuts as set below!!!
  // i.e. Since GEANE is tracking only the primary particle, DCUTE, DCUTM, BCUTE and BCUTM must be put
  // at very high values (10 TeV) in order to calculate properly the energy loss.
  // For a more complete explanation of the chosen values, refer to GEANT manual

  Float_t cut =  1.e-3;                               // 1 MeV cut by default
  Float_t cutd = 1.e4 ;                               // 10 TeV - Threshold for delta-rays
  Float_t cutb = cutd;                                // 10 TeV - Cut for bremsstrahlung
  Float_t tofmax = 1.e10;                             // seconds
  Float_t usrcuts[5] = {0.,0.,0.,0.,0.};              // usercuts
  Float_t gcalpha = 0.999;                            // Optimal value for alpha


  cout<<"Energy straggling area parameter from user set to: "<<gcalpha<<endl;
  if(gcalpha<0.9)
    {
      gcalpha=0.9;
      cout<<"User alpha parameter too low: forced to 0.9"<<endl;
    }

  // set cuts here
  //             GAM ELEC NHAD CHAD MUON EBREM MUHAB EDEL MUDEL MUPA TOFMAX
  gMC3->SetCUTS(cut,  		// CUTGAM = gammas
	 	cut,   	        // CUTELE = electrons
		cut,   	        // CUTNEU = neutral hadrons
		cut,   	        // CUTHAD = charged hadrons
		cut,   	        // CUTMUO = muons
		cutb,  		// BCUTE  = electron bremsstrahlung
		cutb,  		// BCUTM  = muon bremsstrahlung
		cutd,  		// DCUTE  = delta rays by electrons
		cutd,  		// DCUTM  = delta rays by muons
		cutb,   	// PPCUTM = pair production by muons
		tofmax, 	// TOFMAX = time of flight cut
		usrcuts);

  gMC3->SetECut(gcalpha);
}
