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

class StPxlDb;
class StIstDb;
class TGeoHMatrix;
class StMcHit;
class StMcPxlHit;
class StMcIstHit;
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
   void setOutFileName(TString in);

private:
   void projectToVolume(StMcTrack const*,StMcHit const*, double* localProjection, double* localMomentum, TGeoHMatrix const*,double* gPosition,double* newGPosition) const;
   enum ProjectionStatus {GoodProjection,DifferentSensor,OutOfAcceptance};
   ProjectionStatus isOnPxlSensor(double const* localPosition,StMcPxlHit const*,int* correctSensor) const;
   ProjectionStatus isOnIstSensor(double const* localPosition,StMcIstHit const*,int* correctSensor) const;

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

inline void StHftMcHitMover::setOutFileName(TString in) { mOutFileName = in; }

inline char const* StHftMcHitMover::GetCVS() const
{
   static const char cvs[]="" ; 
   return cvs;
}
#endif
