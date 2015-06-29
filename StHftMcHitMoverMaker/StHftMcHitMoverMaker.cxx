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

   mPionsHists->addEvent(mcEvent);
   mKaonsHists->addEvent(mcEvent);
   mProtonsHists->addEvent(mcEvent);

   return kStOK;
}

Int_t StHftMcHitMover::Finish()
{
   mPionsHists->closeFile();
   mKaonsHists->closeFile();
   mProtonsHists->closeFile();

   return kStOK;
}
