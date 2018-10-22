#ifndef MUFLUXSPECTROMETERHIT_H
#define MUFLUXSPECTROMETERHIT_H 1

#include "ShipHit.h"
#include "MufluxSpectrometerPoint.h"
#include "ShipOnlineDataFormat.h"
#include "TVector3.h"

class MufluxSpectrometerHit : public ShipHit {
public:
   /** Default constructor **/
   MufluxSpectrometerHit() = default;

   /** Constructor with arguments
    *@param detID    Detector ID
    *@param digi     digitized/measured TDC
    *@param flags    collection of flags
    **/
   MufluxSpectrometerHit(Int_t detID, Float_t ftdc, Float_t signal_width, uint16_t flag, uint16_t ch);
   MufluxSpectrometerHit(MufluxSpectrometerPoint *p, Double_t t0);
   void MufluxSpectrometerEndPoints(TVector3 &vbot, TVector3 &vtop);
   /** Destructor **/
   virtual ~MufluxSpectrometerHit();

   /** Output to screen **/
   virtual void Print() const;
   Float_t tdc() const { return fdigi; }
   void setInvalid() { flags |= DriftTubes::InValid; }
   bool isValid() const { return !((flags & DriftTubes::InValid) == DriftTubes::InValid); }
   int GetTDC() const { return int((channel & 0xF00) >> 8); }
   bool TDCGood() const
   {
      auto TDC = GetTDC();
      auto AllOK = (flags & DriftTubes::All_OK) == DriftTubes::All_OK;
      uint16_t TDCNotOK = 1 << (TDC + 1);
      return AllOK || !((flags & TDCNotOK) == TDCNotOK);
   }
   bool hasDelay() const { return !((flags & DriftTubes::NoDelay) == DriftTubes::NoDelay); }
   Float_t GetWidth() const { return width; }

private:
   /** Copy constructor **/
   MufluxSpectrometerHit(const MufluxSpectrometerHit &point);
   MufluxSpectrometerHit operator=(const MufluxSpectrometerHit &point);

   Float_t width;
   uint16_t flags;
   uint16_t channel;

   ClassDef(MufluxSpectrometerHit, 7);
};

#endif
