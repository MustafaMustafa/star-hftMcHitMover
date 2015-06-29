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

class StPxlDb;
class StIstDb;
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
   float         mBField;
   StPxlDb*      mPxlDb;
   StIstDb*      mIstDb;

   // temporary for testing in DEV mode
   StHistograms* mPionsHists;
   StHistograms* mKaonsHists;
   StHistograms* mProtonsHists;

   ClassDef(StHftMcHitMover, 0)
};
#endif
