/***************************************************************************
 *
 * $Id$
 *
 * Author: Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *
 * $Log$
 *
 ***************************************************************************/

#include "TGeoMatrix.h"
#include "TGeoManager.h"

#include "St_base/StMessMgr.h"
#include "StarMagField/StarMagField.h"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StarClassLibrary/StParticleDefinition.hh"
#include "StarClassLibrary/StThreeVectorD.hh"
#include "StarClassLibrary/StPhysicalHelixD.hh"

#include "StEvent/StEvent.h"
#include "StMcEvent/StMcTpcHitCollection.hh"
#include "StMcEvent/StMcTrack.hh"
#include "StMcEvent/StMcPxlHit.hh"
#include "StMcEvent/StMcPxlHitCollection.hh"
#include "StMcEvent/StMcIstHit.hh"
#include "StMcEvent/StMcIstHitCollection.hh"
#include "StMcEvent/StMcEvent.hh"
#include "StPxlDbMaker/StPxlDb.h"
#include "StIstDbMaker/StIstDb.h"
// #include "StPxlUtil/StPxlConstants.h"
#include "StIstUtil/StIstConsts.h"

#include "../StHistograms/StHistograms.h"

#include "StHftMcHitMoverMaker.h"

// temp
float const kPxlActiveLengthX = 1.921;
float const kPxlActiveLengthY = 1.9872;

ClassImp(StHftMcHitMover);

StHftMcHitMover::StHftMcHitMover(const Char_t *name) :
   StMaker(name), mPxlDb(NULL), mIstDb(NULL),
   mPionsHists(NULL), mKaonsHists(NULL), mProtonsHists(NULL)
{}

Int_t StHftMcHitMover::Init()
{
   LOG_INFO << "StHftMcHitMover::Init() " << endm;

   if(!mOutFileName.Length()) mOutFileName = "hftMcHitMover";
   mOutFileName = mOutFileName.ReplaceAll(".root","");
   mOutFileName = mOutFileName.ReplaceAll(".event","");
   mOutFileName = mOutFileName.ReplaceAll(".daq","");

   mPionsHists = new StHistograms(Form("%s.Pions",mOutFileName.Data()));
   mKaonsHists = new StHistograms(Form("%s.Kaons",mOutFileName.Data()));
   mProtonsHists = new StHistograms(Form("%s.Protons",mOutFileName.Data()));

   TObjectSet* pxlDbDataSet = static_cast<TObjectSet*>(GetDataSet("pxl_db"));

   if (pxlDbDataSet)
   {
      mPxlDb = static_cast<StPxlDb*>(pxlDbDataSet->GetObject());
   }

   if (!pxlDbDataSet || !mPxlDb)
   {
      LOG_FATAL << "StHftMcHitMover::Init() : no pxlDb - Cannot proceed" << endm;
      return kStFatal;
   }

   TObjectSet* istDbDataSet = static_cast<TObjectSet*>(GetDataSet("ist_db"));

   if (istDbDataSet)
   {
      mIstDb = static_cast<StIstDb*>(istDbDataSet->GetObject());
   }

   if (!istDbDataSet || !mIstDb)
   {
      LOG_FATAL << "StHftMcHitMover::Init() : no istDb - Cannot proceed" << endm;
      return kStErr;
   }

      return kStOK;
}

Int_t StHftMcHitMover::Make()
{
   float center[3] = {0, 0, 0};
   float B[3] = {0, 0, 0};
   StarMagField::Instance()->BField(center, B);
   mBField   = B[2];

   StMcEvent* const mcEvent = (StMcEvent *) GetInputDS("StMcEvent");

   if (!mcEvent)
   {
      LOG_WARN << "StHftMcHitMover::Make() : No StMcEvent" << endm;
      return kStWarn;
   }

   StMcPxlHitCollection* const mcPxlProjCollection = new StMcPxlHitCollection();
   StMcIstHitCollection* const mcIstProjCollection = new StMcIstHitCollection();

   if(!gGeoManager) GetDataBase("VmcGeometry");
   if(!gGeoManager)
   {
     LOG_FATAL << "StHftMcHitMover - FATAL ERROR - gGeoManager is not available" <<endl;
     return kStFatal;
   }

   mPionsHists->addEvent(mcEvent);
   mKaonsHists->addEvent(mcEvent);
   mProtonsHists->addEvent(mcEvent);

   StSPtrVecMcTrack const& tracks = mcEvent->tracks();

   LOG_INFO << "Number of mcTracks = " << tracks.size() << endm;
   LOG_INFO << "Number of TPC Hits in event = " << mcEvent->tpcHitCollection()->numberOfHits() << endm;
   LOG_INFO << "Number of PXL Hits = " << mcEvent->pxlHitCollection()->numberOfHits() << endm;
   LOG_INFO << "Number of IST Hits = " << mcEvent->istHitCollection()->numberOfHits() << endm;

   for (size_t iTrk = 0; iTrk < tracks.size(); ++iTrk)
   {
      StMcTrack* const trk = tracks[iTrk];
      if (!trk) continue;

      // ----------------- Beginning of development cuts
      if (trk->pt() < 0.15 || fabs(trk->momentum().pseudoRapidity()) > 1.0) continue;

      StHistograms* hists = NULL;
      if (trk->geantId() == 8 || trk->geantId() == 9) hists = mPionsHists;
      else if (trk->geantId() == 11 || trk->geantId() == 12) hists = mKaonsHists;
      else if (trk->geantId() == 14 || trk->geantId() == 15) hists = mProtonsHists;
      else continue;
      // ----------------- End of development cuts

      StPtrVecMcPxlHit& pxlHits = trk->pxlHits();
      StPtrVecMcPxlHit newPxlHits;

      for(StMcPxlHitIterator iHit = pxlHits.begin(); iHit != pxlHits.end(); ++iHit)
      {
        StMcPxlHit const* const mcPxlHit = *iHit;

        if(mcPxlHit->localMomentum().mag() < 0.100) continue;

        gGeoManager->RestoreMasterVolume();
        gGeoManager->CdTop();
        gGeoManager->cd(Form("/HALL_1/CAVE_1/TpcRefSys_1/IDSM_1/PXMO_1/PXLA_%i/LADR_%i/PXSI_%i/PLAC_1", mcPxlHit->sector(), mcPxlHit->ladder(), mcPxlHit->sensor()));

        TGeoHMatrix const* volumeM = (TGeoHMatrix*)mPxlDb->geoHMatrixSensorOnGlobal(mcPxlHit->sector(), mcPxlHit->ladder(), mcPxlHit->sensor());
        if (!volumeM) continue;

        mcPxlHit->Print();
        double localProjection[3] = {999.,999.,999.};
        double localMomentum[3]   = {999.,999.,999.};
        double gPosition[3]   = {999.,999.,999.};
        double newGPosition[3]   = {999.,999.,999.};
        projectToVolume(trk,mcPxlHit,localProjection,localMomentum,volumeM,gPosition,newGPosition);

        int correctSensor[3] = {mcPxlHit->sector(),mcPxlHit->ladder(),mcPxlHit->sensor()};
        ProjectionStatus projectionStatus = isOnPxlSensor(localProjection,mcPxlHit,correctSensor);

        if(projectionStatus == DifferentSensor)
        {
          volumeM = (TGeoHMatrix*)mPxlDb->geoHMatrixSensorOnGlobal(correctSensor[0], correctSensor[1], correctSensor[2]);
          if (!volumeM) continue;

          projectToVolume(trk,mcPxlHit,localProjection,localMomentum,volumeM,gPosition,newGPosition);
          projectionStatus = isOnPxlSensor(localProjection,mcPxlHit,correctSensor);
        }

        StMcPxlHit* mcProj = NULL;
        if(projectionStatus == GoodProjection)
        {
          Long_t volumeId = correctSensor[2] * 100 + correctSensor[1] * 10000 + correctSensor[0] * 1000000;

          mcProj = new StMcPxlHit(localProjection, localMomentum, mcPxlHit->dE(),
              mcPxlHit->dS(), mcPxlHit->tof(),
              mcPxlHit->key(), volumeId, trk);

          newPxlHits.push_back(mcProj);
          mcProj->Print();
          mcPxlProjCollection->addHit(mcProj);
          hists->addHits(gPosition,newGPosition);
        }

        StHistograms::Layer layer = (int)mcPxlHit->ladder() == 1? StHistograms::kPxl1 : StHistograms::kPxl2;
        hists->addHits(layer,mcPxlHit,mcProj);
      }

      pxlHits = newPxlHits;

      StPtrVecMcIstHit& istHits = trk->istHits();
      StPtrVecMcIstHit newIstHits;

      for(StMcIstHitIterator iHit = istHits.begin(); iHit != istHits.end(); ++iHit)
      {
        StMcIstHit const* const mcIstHit = *iHit;

        if(mcIstHit->localMomentum().mag() < 0.100) continue;

        gGeoManager->RestoreMasterVolume();
        gGeoManager->CdTop();
        gGeoManager->cd(Form("/HALL_1/CAVE_1/TpcRefSys_1/IDSM_1/IBMO_1/IBAM_%i/IBLM_%i/IBSS_1", static_cast<int>(mcIstHit->ladder()), static_cast<int>(mcIstHit->wafer())));

        TGeoHMatrix const* volumeM = (TGeoHMatrix*)mIstDb->getHMatrixSensorOnGlobal(mcIstHit->ladder(), mcIstHit->wafer());
        if (!volumeM) continue;

        mcIstHit->Print();
        double localProjection[3] = {999.,999.,999.};
        double localMomentum[3]   = {999.,999.,999.};
        double gPosition[3]   = {999.,999.,999.};
        double newGPosition[3]   = {999.,999.,999.};
        projectToVolume(trk,mcIstHit,localProjection,localMomentum,volumeM,gPosition,newGPosition);

        int correctSensor[2] = {mcIstHit->ladder(),mcIstHit->wafer()};
        ProjectionStatus projectionStatus = isOnIstSensor(localProjection,mcIstHit,correctSensor);

        if(projectionStatus == DifferentSensor)
        {
          volumeM = (TGeoHMatrix*)mIstDb->getHMatrixSensorOnGlobal(correctSensor[0],correctSensor[1]);
          if (!volumeM) continue;

          projectToVolume(trk,mcIstHit,localProjection,localMomentum,volumeM,gPosition,newGPosition);
          projectionStatus = isOnIstSensor(localProjection,mcIstHit,correctSensor);
        }

        StMcIstHit* mcProj = NULL;
        if(projectionStatus == GoodProjection)
        {
          Long_t volumeId = correctSensor[1]*10000 + (1+correctSensor[0])*1000000;

          mcProj = new StMcIstHit(localProjection, localMomentum, mcIstHit->dE(),
              mcIstHit->dS(), mcIstHit->tof(),
              mcIstHit->key(), volumeId, trk);

          newIstHits.push_back(mcProj);
          mcProj->Print();
          mcIstProjCollection->addHit(mcProj);
          hists->addHits(gPosition,newGPosition);
        }

        hists->addHits(StHistograms::kIst,mcIstHit,mcProj);
      }

      istHits = newIstHits;
   }

   mcEvent->setPxlHitCollection(mcPxlProjCollection); // StMcEvent::setPxlHitCollection will delete the old collection
   mcEvent->setIstHitCollection(mcIstProjCollection); // StMcEvent::setIstHitCollection will delete the old collection

   LOG_INFO << "Number of new PXL Hits = " << mcEvent->pxlHitCollection()->numberOfHits() << endm;
   LOG_INFO << "Number of new IST Hits = " << mcEvent->istHitCollection()->numberOfHits() << endm;

   return kStOK;
}

void StHftMcHitMover::projectToVolume(StMcTrack const* const trk, StMcHit const* const mcHit,
    double* const localProjection, double* const localMomentum, TGeoHMatrix const* const volumeM,double* const gPosition,double* const newGPosition) const
{
  double const lPosition[3] = {mcHit->position().x(),mcHit->position().y(),mcHit->position().z()};
  double const lMomentum[3] = {mcHit->localMomentum().x(),mcHit->localMomentum().y(),mcHit->localMomentum().z()};
  // double gPosition[3] = {999.,999.,999.};
  double gMomentum[3] = {999.,999.,999.};

  // volumeM->LocalToMaster(lPosition,gPosition);
  // volumeM->LocalToMaster(lMomentum,gMomentum);

  // get global coordinates of hit position and local momentum using ideal geometry
  gGeoManager->GetCurrentMatrix()->LocalToMaster(lPosition,gPosition);
  gGeoManager->GetCurrentMatrix()->LocalToMaster(lMomentum,gMomentum);

  // construct a helix
  StPhysicalHelixD helix(gMomentum, gPosition, mBField * kilogauss, trk->particleDefinition()->charge());

  // project helix to real geometry sensor
  double const* rotation = volumeM->GetRotationMatrix();
  double const* translation = volumeM->GetTranslation();
  StThreeVectorD const sensorNormal(rotation[1], rotation[4], rotation[7]);
  StThreeVectorD const sensorCenter(translation);

  double const s = helix.pathLength(sensorCenter, sensorNormal);

  StThreeVectorD gMcProj = helix.at(s);
  newGPosition[0] = gMcProj.x();
  newGPosition[1] = gMcProj.y();
  newGPosition[2] = gMcProj.z();

  volumeM->MasterToLocal(helix.at(s).xyz(), localProjection);
  volumeM->MasterToLocal(helix.momentumAt(s,mBField * kilogauss).xyz(), localMomentum);
}

StHftMcHitMover::ProjectionStatus StHftMcHitMover::isOnPxlSensor(double const* const localPosition, StMcPxlHit const* const mcPxlHit,int* correctSensor) const
{
   LOG_INFO << " projection x/y/z = " << localPosition[0] << " / " << localPosition[1] << " / " << localPosition[2] << "\n" << endm;
   int shiftLadder = 0;
   if (localPosition[0] > kPxlActiveLengthX/2.) ++shiftLadder;
   else if (localPosition[0] < -kPxlActiveLengthX/2.) --shiftLadder;

   int shiftSensor = 0;
   if (localPosition[2] > kPxlActiveLengthY/2.) ++shiftSensor;
   else if (localPosition[2] < -kPxlActiveLengthY/2.) --shiftSensor;

   LOG_INFO << "Sector/Ladder/Sensor = " << (int)mcPxlHit->sector() <<" / "<<(int)mcPxlHit->ladder() << " / " << (int)mcPxlHit->sensor() <<endm;
   LOG_INFO << "shiftLadder/shiftSensor = " << shiftLadder << " / "<<shiftSensor << endm;

   ProjectionStatus projectionStatus = GoodProjection;

   int shiftSector = 0;
   if(shiftLadder)
   {
     projectionStatus = DifferentSensor;

     if(mcPxlHit->ladder() == 1)
     {
       shiftSector = shiftLadder; // PXL1 rPhi is clockwise, same as sector numbering
       correctSensor[1] = 1;
     }
     else
     {
       switch (mcPxlHit->ladder() - shiftLadder) // PXL2 rPhi is ccw, opposite to sector numbering
       {
         case 5:
           ++shiftSector;
           correctSensor[1] = 2;
           break;
         case 1:
           --shiftSector;
           correctSensor[1] = 4;
           break;
         default:
           correctSensor[1] = mcPxlHit->ladder() - shiftLadder;
           break;
       }
     }
   }

   if(shiftSector)
   {
     switch (mcPxlHit->sector() + shiftSector)
     {
       case 0:
         correctSensor[0] = 10;
         break;
       case 11:
         correctSensor[0] = 1;
         break;
       default:
         correctSensor[0] = mcPxlHit->sector() + shiftSector;
     }
   }

   if(shiftSensor)
   {
     projectionStatus = DifferentSensor;

     switch (mcPxlHit->sensor() + shiftSensor)
     {
       case 0:
       case 11:
         projectionStatus = OutOfAcceptance;
         break;
       default:
         correctSensor[2] = mcPxlHit->sensor() + shiftSensor;
     }
   }

   LOG_INFO << "New Sector/Ladder/Sensor = " << correctSensor[0] <<" / "<< correctSensor[1] << " / " << correctSensor[2] <<endm;
   LOG_INFO << "ProjectionStatus = " << projectionStatus << "\n" << endm;
   return projectionStatus;
}

StHftMcHitMover::ProjectionStatus StHftMcHitMover::isOnIstSensor(double const* localPosition,StMcIstHit const* const mcIstHit,int* correctSensor) const
{
   LOG_INFO << " projection x/y/z = " << localPosition[0] << " / " << localPosition[1] << " / " << localPosition[2] << "\n" << endm;
   int shiftLadder = 0;
   if (localPosition[0] > kIstSensorActiveSizeRPhi/2.) ++shiftLadder;
   else if (localPosition[0] < -kIstSensorActiveSizeRPhi/2.) --shiftLadder;

   int shiftSensor = 0;
   if (localPosition[2] > kIstSensorActiveSizeZ/2.) ++shiftSensor;
   else if (localPosition[2] < -kIstSensorActiveSizeZ/2.) --shiftSensor;

   LOG_INFO << "Ladder/Wafer " <<(int)mcIstHit->ladder() << " / " << (int)mcIstHit->wafer() <<endm;
   LOG_INFO << "shiftLadder/shiftSensor = " << shiftLadder << " / "<<shiftSensor << endm;

   ProjectionStatus projectionStatus = GoodProjection;

   if(shiftLadder)
   {
     projectionStatus = DifferentSensor;

     switch (mcIstHit->ladder() - shiftLadder) // rPhi goes +ve is opposite to ladder numbering
     {
       case 0:
         correctSensor[0] = 24;
         break;
       case 25:
         correctSensor[0] = 1;
         break;
       default:
         correctSensor[0] = mcIstHit->ladder() - shiftLadder;
     }
   }

   if(shiftSensor)
   {
     projectionStatus = DifferentSensor;

     switch (mcIstHit->wafer() + shiftSensor)
     {
       case 0:
       case 7:
         projectionStatus = OutOfAcceptance;
         break;
       default:
         correctSensor[1] = mcIstHit->wafer() + shiftSensor;
     }
   }

   LOG_INFO << "New Ladder/Wafer = " << correctSensor[0] <<" / "<< correctSensor[1] <<endm;
   LOG_INFO << "ProjectionStatus = " << projectionStatus << "\n" << endm;
   return projectionStatus;

  // return fabs(localPosition[0]) < kIstSensorActiveSizeRPhi/2. && fabs(localPosition[2]) < kIstSensorActiveSizeZ/2.;
}

Int_t StHftMcHitMover::Finish()
{
   mPionsHists->closeFile();
   mKaonsHists->closeFile();
   mProtonsHists->closeFile();

   return kStOK;
}
