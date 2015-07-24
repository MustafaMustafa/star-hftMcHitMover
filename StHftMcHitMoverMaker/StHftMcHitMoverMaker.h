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

#include "TString.h"

#include "StMaker.h"
// #include "StPxlUtil/StPxlConstants.h"
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
   virtual Int_t  Make();
   virtual Int_t  Finish();

   virtual char const* GetCVS() const;
   void setOutFileName(TString in) { mOutFileName = in; }

private:
   void projectToVolume(StMcTrack const*,StMcHit const*, double* localProjection, double* localMomentum, TGeoHMatrix const*) const;
   bool isOnPxlSensor(double const* localPosition) const;
   bool isOnIstSensor(double const* localPosition) const;

   TString       mOutFileName;
   float         mBField;
   StPxlDb*      mPxlDb;
   StIstDb*      mIstDb;

   // temporary for testing in DEV mode
   StHistograms* mPionsHists;
   StHistograms* mKaonsHists;
   StHistograms* mProtonsHists;

   ClassDef(StHftMcHitMover, 0)
};

inline char const* StHftMcHitMover::GetCVS() const
{
   static const char cvs[]="Tag $Name:$ $Id:$ built "__DATE__" " __TIME__ ; 
   return cvs;
}

inline bool StHftMcHitMover::isOnPxlSensor(double const* const localPosition) const
{
   // return fabs(localPosition[0]) < StPxlConsts::kPxlActiveLengthX/2. && fabs(localPosition[2]) < StPxlConsts::kPxlActiveLengthY/2.;
   return fabs(localPosition[0]) < 1.921/2. && fabs(localPosition[2]) < 1.9872/2.;
}

inline bool StHftMcHitMover::isOnIstSensor(double const* const localPosition) const
{
  return fabs(localPosition[0]) < kIstSensorActiveSizeRPhi/2. && fabs(localPosition[2]) < kIstSensorActiveSizeZ/2.;
}

#endif
