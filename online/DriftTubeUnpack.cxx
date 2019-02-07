#include <cassert>
#include <unordered_map>
#include <iostream>
#include <bitset>

// ROOT headers
#include "TClonesArray.h"
#include "ROOT/TSeq.hxx"
#include "ROOT/RVec.hxx"

// Fair headers
#include "FairRootManager.h"
#include "FairRunOnline.h"
#include "FairLogger.h"

// SHiP headers
#include "DriftTubeUnpack.h"
#include "MufluxSpectrometerHit.h"
#include "ScintillatorHit.h"
#include "ShipOnlineDataFormat.h"

using DriftTubes::ChannelId;

// DriftTubeUnpack: Constructor
DriftTubeUnpack::DriftTubeUnpack()
   : fRawTubes(new TClonesArray("MufluxSpectrometerHit")), fRawLateTubes(new TClonesArray("MufluxSpectrometerHit")),
     fRawScintillator(new TClonesArray("ScintillatorHit")), fRawBeamCounter(new TClonesArray("ScintillatorHit")),
     fRawMasterTrigger(new TClonesArray("ScintillatorHit")), fRawTriggers(new TClonesArray("ScintillatorHit")),
     fPartitionId(0x0C00)
{
}

// Virtual DriftTubeUnpack: Public method
DriftTubeUnpack::~DriftTubeUnpack() = default;

// Init: Public method
Bool_t DriftTubeUnpack::Init()
{
   Register();
   return kTRUE;
}

// Register: Protected method
void DriftTubeUnpack::Register()
{
   LOG(INFO) << "DriftTubeUnpack : Registering..." << FairLogger::endl;
   auto fMan = FairRootManager::Instance();
   if (!fMan) {
      return;
   }
   fMan->Register("Digi_MufluxSpectrometerHits", "DriftTubes", fRawTubes.get(), kTRUE);
   fMan->Register("Digi_LateMufluxSpectrometerHits", "DriftTubes", fRawLateTubes.get(), kTRUE);
   fMan->Register("Digi_Scintillators", "DriftTubes", fRawScintillator.get(), kTRUE);
   fMan->Register("Digi_BeamCounters", "DriftTubes", fRawBeamCounter.get(), kTRUE);
   fMan->Register("Digi_MasterTrigger", "DriftTubes", fRawMasterTrigger.get(), kTRUE);
   fMan->Register("Digi_Triggers", "DriftTubes", fRawTriggers.get(), kTRUE);
}

// DoUnpack: Public method
Bool_t DriftTubeUnpack::DoUnpack(Int_t *data, Int_t size)
{
   LOG(DEBUG) << "DriftTubeUnpack : Unpacking frame... size/bytes = " << size << FairLogger::endl;

   auto df = reinterpret_cast<DataFrame *>(data);
   assert(df->header.size == size);
   switch (df->header.frameTime) {
   case SoS:
      LOG(DEBUG) << "DriftTubeUnpacker: SoS frame." << FairLogger::endl;
      for (auto i : ROOT::MakeSeq(size)) {
         if (i % 4 == 0) {
            std::cout << ' ';
         } else if (i % 16 == 0) {
            std::cout << '\n';
         }
         std::cout << std::hex << +data[i] << std::dec;
      }
      std::cout << std::endl;
      return kTRUE;
   case EoS: LOG(DEBUG) << "DriftTubeUnpacker: EoS frame." << FairLogger::endl; return kTRUE;
   default: break;
   }
   LOG(DEBUG) << "Sequential trigger number " << df->header.timeExtent << FairLogger::endl;
   auto nhits = df->getHitCount();
   int nhitsScintillator = 0;
   int nhitsBeamCounter = 0;
   int nhitsMasterTrigger = 0;
   int nhitsTriggers = 0;
   auto flags = df->header.flags;
   int skipped = 0;
   int trigger = 0;
   int expected_triggers = 5;
   if ((flags & DriftTubes::All_OK) == DriftTubes::All_OK) {
      LOG(DEBUG) << "All TDCs are OK" << FairLogger::endl;
   } else {
      LOG(DEBUG) << "Not all TDCs are OK:" << std::bitset<16>(flags) << FairLogger::endl;
      for (auto i : ROOT::MakeSeq(5)) {
         if ((flags & 1 << (i + 1)) == 1 << (i + 1)) {
            expected_triggers--;
            LOG(WARNING) << "TDC " << i << " NOT OK" << FairLogger::endl;
         } else {
            LOG(DEBUG) << "TDC " << i << " OK" << FairLogger::endl;
         }
      }
   }
   std::vector<RawDataHit> hits(df->hits, df->hits + nhits);
   std::unordered_map<int, uint16_t> triggerTime;
   uint16_t master_trigger_time = 0;
   std::vector<std::pair<int, uint16_t>> trigger_times;
   for (auto &&hit : hits) {
      auto channel = reinterpret_cast<ChannelId *>(&(hit.channelId));
      auto TDC = channel->TDC;
      auto detectorId = channel->GetDetectorId();
      auto hit_time = hit.hitTime;
      if (!detectorId) {
         if (channel->edge == 0) {
            trigger++;
            triggerTime[TDC] =
               (triggerTime.find(TDC) != triggerTime.end()) ? std::min(hit_time, triggerTime[TDC]) : hit_time;
            trigger_times.emplace_back(TDC, hit_time);
         }
         new ((*fRawTriggers)[nhitsTriggers])
            ScintillatorHit(detectorId, 0.098 * Float_t(hit_time), flags, hit.channelId);
         nhitsTriggers++;
      } else if (detectorId == 1) {
         if (channel->edge == 0) {
            // Use the earliest if there are several
            if (nhitsMasterTrigger == 0 || hit_time < master_trigger_time) {
               master_trigger_time = hit_time;
            }
         }
         new ((*fRawMasterTrigger)[nhitsMasterTrigger])
            ScintillatorHit(detectorId, 0.098 * Float_t(hit_time), flags, hit.channelId);
         nhitsMasterTrigger++;
      } else if (detectorId == -1) {
         // beam counter
         new ((*fRawBeamCounter)[nhitsBeamCounter])
            ScintillatorHit(detectorId, 0.098 * Float_t(hit_time), flags, hit.channelId);
         nhitsBeamCounter++;
      } else if (detectorId == 6 || detectorId == 7) {
         // beam counter
         new ((*fRawScintillator)[nhitsScintillator])
            ScintillatorHit(detectorId, 0.098 * Float_t(hit_time), flags, hit.channelId);
         nhitsScintillator++;
      }
   }

   return kTRUE;
}

// Reset: Public method
void DriftTubeUnpack::Reset()
{
   LOG(DEBUG) << "DriftTubeUnpack : Clearing Data Structure" << FairLogger::endl;
   fRawTubes->Clear();
   fRawLateTubes->Clear();
   fRawScintillator->Clear();
   fRawBeamCounter->Clear();
   fRawMasterTrigger->Clear();
   fRawTriggers->Clear();
}

ClassImp(DriftTubeUnpack)
