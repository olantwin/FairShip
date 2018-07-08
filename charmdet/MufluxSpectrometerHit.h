#ifndef MUFLUXSPECTROMETERHIT_H
#define MUFLUXSPECTROMETERHIT_H 1

#include "ShipHit.h"
#include "MufluxSpectrometerPoint.h"
#include "TObject.h"
#include "TVector3.h"

class MufluxSpectrometerHit : public ShipHit {
public:
   /** Default constructor **/
   MufluxSpectrometerHit();

   /** Constructor with arguments
    *@param detID    Detector ID
    *@param digi      digitized/measured TDC
    *@param flag      True/False, false if there is another hit with smaller tdc
    **/
   MufluxSpectrometerHit(Int_t detID, Float_t tdc);
   MufluxSpectrometerHit(MufluxSpectrometerPoint *p, Double_t t0);
   void MufluxSpectrometerEndPoints(TVector3 &vbot, TVector3 &vtop);
   /** Destructor **/
   virtual ~MufluxSpectrometerHit();

   /** Output to screen **/
   virtual void Print() const;
   Float_t tdc() const { return fdigi; }
   void setInvalid() { flags &= ~DriftTubeFlag::Valid; }
   bool isValid() const { return flag & DriftTubeFlag::Valid; }

private:
   /** Copy constructor **/
   MufluxSpectrometerHit(const MufluxSpectrometerHit &point);
   MufluxSpectrometerHit operator=(const MufluxSpectrometerHit &point);

   uint16_t flags; ///< flag

   ClassDef(MufluxSpectrometerHit, 3);
};

enum DriftTubeFlag : uint16_t {
   None = 0,
   Valid = 1,
   UnmatchedHeaderTrailer = 1 << 15,
}

#endif
