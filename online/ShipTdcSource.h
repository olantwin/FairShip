#ifndef ONLINE_SHIPTDCSOURCE_H
#define ONLINE_SHIPTDCSOURCE_H

#include "FairOnlineSource.h"
#include "TObjArray.h"
#include "TFile.h"

#include "FairUnpack.h"

class FairEventHeader;

class ShipTdcSource : public FairOnlineSource {
public:
   ShipTdcSource();
   explicit ShipTdcSource(TString filename);
   ShipTdcSource(const ShipTdcSource &source);
   virtual ~ShipTdcSource();

   virtual Bool_t Init();
   virtual Int_t ReadEvent(UInt_t = 0); // Read frame by frame
   virtual void Close();
   void FillEventHeader(FairEventHeader *feh);

protected:
   Bool_t Unpack(Int_t *data, Int_t size, uint16_t partitionId);
   Int_t UnpackEventFrame(Int_t *data, Int_t size);
   TFile *fIn;
   unsigned char buffer[UINT16_MAX];
   Double_t fEventTime = 0;

   TString fFilename;

   ClassDef(ShipTdcSource, 1)
};

#endif
