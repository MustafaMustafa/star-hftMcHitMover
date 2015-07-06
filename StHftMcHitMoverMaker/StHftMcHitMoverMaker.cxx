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
#include "StMcEvent/StMcTrack.hh"
#include "StMcEvent/StMcPxlHit.hh"
#include "StMcEvent/StMcPxlHitCollection.hh"
#include "StMcEvent/StMcIstHit.hh"
#include "StMcEvent/StMcIstHitCollection.hh"
#include "StMcEvent/StMcEvent.hh"
#include "StPxlDbMaker/StPxlDb.h"
#include "StIstDbMaker/StIstDb.h"

#include "../StHistograms/StHistograms.h"

#include "StHftMcHitMoverMaker.h"

ClassImp(StHftMcHitMover);

StHftMcHitMover::StHftMcHitMover(const Char_t *name) :
   StMaker(name), mPxlDb(NULL), mIstDb(NULL),
   mPionsHists(NULL), mKaonsHists(NULL), mProtonsHists(NULL)
{
   mPionsHists = new StHistograms("idealPions");
   mKaonsHists = new StHistograms("idealKaons");
   mProtonsHists = new StHistograms("idealProtons");
}

Int_t StHftMcHitMover::Init()
{
   LOG_INFO << "StHftMcHitMover::Init() " << endm;
   return kStOK;
}

Int_t StHftMcHitMover::InitRun(Int_t runNumber)
{
   TObjectSet* pxlDbDataSet = static_cast<TObjectSet*>(GetDataSet("pxl_db"));

   if (pxlDbDataSet)
   {
      mPxlDb = static_cast<StPxlDb*>(pxlDbDataSet->GetObject());
   }

   if (!pxlDbDataSet || !mPxlDb)
   {
      LOG_FATAL << "StHftMcHitMover::InitRun() : no pxlDb - Cannot proceed" << endm;
      return kStFatal;
   }

   TObjectSet* istDbDataSet = static_cast<TObjectSet*>(GetDataSet("ist_db"));

   if (istDbDataSet)
   {
      mIstDb = static_cast<StIstDb*>(istDbDataSet->GetObject());
   }

   if (!istDbDataSet || !mIstDb)
   {
      LOG_FATAL << "StHftMcHitMover::InitRun() : no istDb - Cannot proceed" << endm;
      return kStErr;
   }

   float center[3] = {0, 0, 0};
   float B[3] = {0, 0, 0};
   StarMagField::Instance()->BField(center, B);
   mBField   = B[2];

   return kStOK;
}

Int_t StHftMcHitMover::Make()
{
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

   for (size_t iTrk = 0; iTrk < tracks.size(); ++iTrk)
   {
      StMcTrack* const trk = tracks[iTrk];
      if (!trk) continue;

      if (trk->pt() < 0.15 || fabs(trk->momentum().pseudoRapidity()) > 1.0) continue;

      StHistograms* hists = NULL;
      if (trk->geantId() == 8 || trk->geantId() == 9) hists = mPionsHists;
      else if (trk->geantId() == 11 || trk->geantId() == 12) hists = mKaonsHists;
      else if (trk->geantId() == 14 || trk->geantId() == 15) hists = mProtonsHists;
      else continue;

      StPtrVecMcPxlHit& pxlHits = trk->pxlHits();
      StPtrVecMcPxlHit newPxlHits;

      for(StMcPxlHitIterator iHit = pxlHits.begin(); iHit != pxlHits.end(); ++iHit)
      {
        StMcPxlHit const* const mcPxlHit = *iHit;

        if(mcPxlHit->localMomentum().mag() < 0.100) continue;

        TGeoHMatrix const* const volumeM = (TGeoHMatrix*)mPxlDb->geoHMatrixSensorOnGlobal(mcPxlHit->sector(), mcPxlHit->ladder(), mcPxlHit->sensor());
        if (!volumeM) continue;

        gGeoManager->RestoreMasterVolume();
        gGeoManager->CdTop();
        gGeoManager->cd(Form("/HALL_1/CAVE_1/TpcRefSys_1/IDSM_1/PXMO_1/PXLA_%i/LADR_%i/PXSI_%i/PLAC_1", mcPxlHit->sector(), mcPxlHit->ladder(), mcPxlHit->sensor()));

        double localProjection[3] = {999.,999.,999.};
        double localMomentum[3]   = {999.,999.,999.};
        projectToVolume(trk,mcPxlHit,localProjection,localMomentum,volumeM);

        if(!isOnPxlSensor(localProjection)) continue;

        StMcPxlHit* const mcProj = new StMcPxlHit(localProjection, localMomentum, mcPxlHit->dE(),
                                            mcPxlHit->dS(), mcPxlHit->tof(),
                                            mcPxlHit->key(), mcPxlHit->volumeId(), trk);

        StHistograms::Layer layer = (int)mcPxlHit->ladder() == 1? StHistograms::kPxl1 : StHistograms::kPxl2;
        hists->addHits(layer,mcPxlHit,mcProj);

        newPxlHits.push_back(mcProj);
        mcPxlProjCollection->addHit(mcProj);
      }

      pxlHits = newPxlHits;

      StPtrVecMcIstHit& istHits = trk->istHits();
      StPtrVecMcIstHit newIstHits;

      for(StMcIstHitIterator iHit = istHits.begin(); iHit != istHits.end(); ++iHit)
      {
        StMcIstHit const* const mcIstHit = *iHit;

        if(mcIstHit->localMomentum().mag() < 0.100) continue;

        TGeoHMatrix const* const volumeM = (TGeoHMatrix*)mIstDb->getHMatrixSensorOnGlobal(mcIstHit->ladder(), mcIstHit->wafer());
        if (!volumeM) continue;

        gGeoManager->RestoreMasterVolume();
        gGeoManager->CdTop();
        gGeoManager->cd(Form("/HALL_1/CAVE_1/TpcRefSys_1/IDSM_1/IBMO_1/IBAM_%i/IBLM_%i/IBSS_1", static_cast<int>(mcIstHit->ladder()), static_cast<int>(mcIstHit->wafer())));

        double localProjection[3] = {999.,999.,999.};
        double localMomentum[3]   = {999.,999.,999.};
        projectToVolume(trk,mcIstHit,localProjection,localMomentum,volumeM);

        if(!isOnIstSensor(localProjection)) continue;

        StMcIstHit* const mcProj = new StMcIstHit(localProjection, localMomentum, mcIstHit->dE(),
                                            mcIstHit->dS(), mcIstHit->tof(),
                                            mcIstHit->key(), mcIstHit->volumeId(), trk);

        hists->addHits(StHistograms::kIst,mcIstHit,mcProj);

        newIstHits.push_back(mcProj);
        mcIstProjCollection->addHit(mcProj);
      }

      istHits = newIstHits;
   }

   mcEvent->setPxlHitCollection(mcPxlProjCollection); // StMcEvent::setPxlHitCollection will delete the old collection
   mcEvent->setIstHitCollection(mcIstProjCollection); // StMcEvent::setIstHitCollection will delete the old collection

   return kStOK;
}

void StHftMcHitMover::projectToVolume(StMcTrack const* const trk, StMcHit const* const mcHit,
    double* const localProjection, double* const localMomentum, TGeoHMatrix const* const volumeM) const
{
  double const lPosition[3] = {mcHit->position().x(),mcHit->position().y(),mcHit->position().z()};
  double const lMomentum[3] = {mcHit->localMomentum().x(),mcHit->localMomentum().y(),mcHit->localMomentum().z()};
  double gPosition[3] = {999.,999.,999.};
  double gMomentum[3] = {999.,999.,999.};

  // volumeM->LocalToMaster(lPosition,gPosition);
  // volumeM->LocalToMaster(lMomentum,gMomentum);
  gGeoManager->GetCurrentMatrix()->LocalToMaster(lPosition,gPosition);
  gGeoManager->GetCurrentMatrix()->LocalToMaster(lMomentum,gMomentum);

  StPhysicalHelixD helix(gMomentum, gPosition, mBField * kilogauss, trk->particleDefinition()->charge());

  double const* rotation = volumeM->GetRotationMatrix();
  double const* translation = volumeM->GetTranslation();
  StThreeVectorD const sensorNormal(rotation[1], rotation[4], rotation[7]);
  StThreeVectorD const sensorCenter(translation);

  double const s = helix.pathLength(sensorCenter, sensorNormal);

  volumeM->MasterToLocal(helix.at(s).xyz(), localProjection);
  volumeM->MasterToLocal(helix.momentumAt(s,mBField * kilogauss).xyz(), localMomentum);
}

Int_t StHftMcHitMover::Finish()
{
   mPionsHists->closeFile();
   mKaonsHists->closeFile();
   mProtonsHists->closeFile();

   return kStOK;
}
