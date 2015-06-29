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

#ifndef StHftMcHitMover_hh
#define StHftMcHitMover_hh

#include "StMaker.h"
#include "StPxlUtil/StPxlConstants.h"
#include "StIstUtil/StIstConsts.h"

class StPxlDb;
class StIstDb;
class TGeoHMatrix;
class StMcHit;
class StMcTrack;
class StHistograms;

class StHftMcHitMover : public StMaker
{
public:
   StHftMcHitMover(const Char_t *name = "hftMcHitMover");
   ~StHftMcHitMover() {}

   virtual Int_t  Init();
   virtual Int_t  InitRun(Int_t runNumber);
   virtual Int_t  Make();
   virtual Int_t  Finish();

private:
   void projectToVolume(StMcTrack const*,StMcHit const*, double* localProjection, double* localMomentum, TGeoHMatrix const*) const;
   bool isOnPxlSensor(double const* localPosition) const;
   bool isOnIstSensor(double const* localPosition) const;

   float         mBField;
   StPxlDb*      mPxlDb;
   StIstDb*      mIstDb;

   // temporary for testing in DEV mode
   StHistograms* mPionsHists;
   StHistograms* mKaonsHists;
   StHistograms* mProtonsHists;

   ClassDef(StHftMcHitMover, 0)
};

inline bool StHftMcHitMover::isOnPxlSensor(double const* const localPosition) const
{
   return fabs(localPosition[0]) < StPxlConsts::kPxlActiveLengthX/2. && fabs(localPosition[2]) < StPxlConsts::kPxlActiveLengthY/2.;
}

inline bool StHftMcHitMover::isOnIstSensor(double const* const localPosition) const
{
  return fabs(localPosition[0]) < kIstSensorActiveSizeRPhi/2. && fabs(localPosition[2]) < kIstSensorActiveSizeZ/2.;
}

#endif
