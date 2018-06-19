#include "ShipTdcSource.h"
#include "FairLogger.h"
#include "ShipOnlineDataFormat.h"
#include "DriftTubeUnpack.h"

ShipTdcSource::ShipTdcSource() : FairOnlineSource(), fFilename("tdcdata.bin") {}

ShipTdcSource::ShipTdcSource(std::string filename) : FairOnlineSource(), fFilename(filename) {}

ShipTdcSource::ShipTdcSource(const ShipTdcSource &source) : FairOnlineSource(source) {}

ShipTdcSource::~ShipTdcSource() {}

Bool_t ShipTdcSource::Init()
{
   fIn = new std::ifstream(fFilename, std::ios::binary);
   return true;
}

void ShipTdcSource::Close()
{
   fIn->close();
}

Int_t ShipTdcSource::ReadEvent(UInt_t)
{
   auto df = new (buffer) DataFrame();
   if (fIn->read(reinterpret_cast<char *>(df), sizeof(DataFrame))) {
      auto nhits = df->getHitCount();
      if (fIn->read(reinterpret_cast<char *>(df->hits), nhits * sizeof(RawDataHit))) {

         if (Unpack(reinterpret_cast<Int_t *>(&buffer), sizeof(DataFrame) + nhits * sizeof(RawDataHit),
                    df->header.partitionId)) {
            return 0;
         }
         LOG(WARNING) << "ShipTdcSource: Failed to Unpack." << FairLogger::endl;
         return 3;
      }
      LOG(WARNING) << "ShipTdcSource: Failed to read hits." << FairLogger::endl;
      return 2;
   }
   return 1;
   LOG(WARNING) << "ShipTdcSource: Failed to read header." << FairLogger::endl;
}

Bool_t ShipTdcSource::Unpack(Int_t *data, Int_t size, uint16_t partitionId)
{

   for (TObject *item : *fUnpackers) {
      auto unpacker = static_cast<DriftTubeUnpack *>(item); // TODO need ShipUnpack with partitionId?
      if (unpacker->GetPartition() == partitionId) {
         return unpacker->DoUnpack(data, size);
      }
   }

   LOG(WARNING) << "ShipTdcSource: Failed to find suitable unpacker." << FairLogger::endl;
   return kFALSE;
}

ClassImp(ShipTdcSource)
