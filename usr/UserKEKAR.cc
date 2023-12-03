// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iostream>
#include <sstream>

#include <UnpackerManager.hh>

#include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoRawHit.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "DCGeomMan.hh"

#define FHitBranch 0 // make FiberHit branches (becomes heavy)

namespace
{
using namespace root;
const auto qnan = TMath::QuietNaN();
auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
auto& gRM = RMAnalyzer::GetInstance();
auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
class UserKEKAR : public VEvent
{
private:
  RawData*      rawData;

public:
  UserKEKAR();
  ~UserKEKAR();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserKEKAR::ClassName()
{
  static TString s_name("UserKEKAR");
  return s_name;
}

//_____________________________________________________________________________
UserKEKAR::UserKEKAR()
  : VEvent(),
    rawData(new RawData)
{
}

//_____________________________________________________________________________
UserKEKAR::~UserKEKAR()
{
  if(rawData) delete rawData;
}

//_____________________________________________________________________________
struct Event
{
 
  Int_t evnum;
  Int_t spill;

  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];
  
  Int_t t1nhits;
  Int_t t1hitpat[MaxHits];
  Double_t t1a[NumOfSegKEKART1];
  Double_t t1t[NumOfSegKEKART1][MaxDepth];

  Int_t t2nhits;
  Int_t t2hitpat[MaxHits];
  Double_t t2a[NumOfSegKEKART2];
  Double_t t2t[NumOfSegKEKART2][MaxDepth];

  Int_t t3nhits;
  Int_t t3hitpat[MaxHits];
  Double_t t3a[NumOfSegKEKART3];  
  Double_t t3t[NumOfSegKEKART3][MaxDepth];

  Int_t t4nhits;
  Int_t t4hitpat[MaxHits];
  Double_t t4a[NumOfSegKEKART4];
  Double_t t4t[NumOfSegKEKART4][MaxDepth];

  Int_t sacnhits;
  Int_t sachitpat[MaxHits];
  Double_t saca[NumOfSegKEKARE90SAC];
  Double_t sact[NumOfSegKEKARE90SAC][MaxDepth];
  
  Int_t sacsumnhits;
  Double_t sacsuma[NumOfSegKEKARE90SACSUM];
  Double_t sacsumt[NumOfSegKEKARE90SACSUM][MaxDepth];
  
  Int_t bacnhits;
  Int_t bachitpat[MaxHits]; 
  Double_t baca[NumOfSegKEKARE72BAC];
  Double_t bact[NumOfSegKEKARE72BAC][MaxDepth];

  Int_t bacsumnhits;
  Double_t bacsuma[NumOfSegKEKARE72BACSUM];
  Double_t bacsumt[NumOfSegKEKARE72BACSUM][MaxDepth];

  Int_t kvcnhits;
  Int_t kvchitpat[MaxHits];
  Double_t kvca[NumOfSegKEKARE72KVC];
  Double_t kvct[NumOfSegKEKARE72KVC][MaxDepth];
  
  Int_t kvcsumnhits;
  Int_t kvcsumhitpat[MaxHits];
  Double_t kvcsuma[NumOfSegKEKARE72KVCSUM];
  Double_t kvcsumt[NumOfSegKEKARE72KVCSUM][MaxDepth];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum      = 0;
  spill      = 0;
  t1nhits   = 0;
  t2nhits   = 0;
  t3nhits   = 0;
  t4nhits   = 0;
  sacnhits   = 0;
  sacsumnhits   = 0;
  bacnhits   = 0;
  bacsumnhits   = 0;
  kvcnhits   = 0;
  kvcsumnhits   = 0;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<MaxHits; ++it){
    t1hitpat[it]   = -1;
    t2hitpat[it]   = -1;
    t3hitpat[it]   = -1;
    t4hitpat[it]   = -1;
    bachitpat[it]   = -1;
    sachitpat[it]   = -1;
    kvchitpat[it]   = -1;
    kvcsumhitpat[it]   = -1;
  }

  for(Int_t it=0; it<NumOfSegKEKART1; ++it){
    t1a[it] = qnan;
    for(Int_t m = 0; m<MaxDepth; ++m){
      t1t[it][m] = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegKEKART2; ++it){
    t2a[it] = qnan;
    for(Int_t m = 0; m<MaxDepth; ++m){
      t2t[it][m] = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegKEKART3; ++it){
    t3a[it] = qnan;    
    for(Int_t m = 0; m<MaxDepth; ++m){
      t3t[it][m] = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegKEKART4; ++it){
    t4a[it] = qnan;
    for(Int_t m = 0; m<MaxDepth; ++m){
      t4t[it][m] = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegKEKARE90SAC; ++it){
    saca[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      sact[it][m] = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegKEKARE90SACSUM; ++it){
    sacsuma[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      sacsumt[it][m] = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegKEKARE72BAC; ++it){
    baca[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      bact[it][m] = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegKEKARE72BACSUM; ++it){
    bacsuma[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      bacsumt[it][m] = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegKEKARE72KVC; ++it){
    kvca[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      kvct[it][m] = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegKEKARE72KVCSUM; ++it){
    kvcsuma[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      kvcsumt[it][m] = qnan;
    }
  }

}

//_____________________________________________________________________________
namespace root
{
Event  event;
TH1*   h[MaxHist];
TTree* tree;
enum eDetHid {
  T1Hid     = 10000,
  T2Hid     = 20000,
  T3Hid     = 30000,
  T4Hid     = 40000,
  SACHid    = 50000,
  BACHid    = 60000,
  KVCHid    = 70000
};
}

//_____________________________________________________________________________
Bool_t
UserKEKAR::ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserKEKAR::ProcessingNormal()
{
  static const auto MinTdcT1 = gUser.GetParameter("TdcKEKART1", 0);
  static const auto MaxTdcT1 = gUser.GetParameter("TdcKEKART1", 1);
  static const auto MinTdcT2 = gUser.GetParameter("TdcKEKART2", 0);
  static const auto MaxTdcT2 = gUser.GetParameter("TdcKEKART2", 1);
  static const auto MinTdcT3 = gUser.GetParameter("TdcKEKART3", 0);
  static const auto MaxTdcT3 = gUser.GetParameter("TdcKEKART3", 1);
  static const auto MinTdcT4 = gUser.GetParameter("TdcKEKART4", 0);
  static const auto MaxTdcT4 = gUser.GetParameter("TdcKEKART4", 1);
  static const auto MinTdcSAC = gUser.GetParameter("TdcKEKARE90SAC", 0);
  static const auto MaxTdcSAC = gUser.GetParameter("TdcKEKARE90SAC", 1);
  static const auto MinTdcBAC = gUser.GetParameter("TdcKEKARE72BAC", 0);
  static const auto MaxTdcBAC = gUser.GetParameter("TdcKEKARE72BAC", 1);
  static const auto MinTdcKVC = gUser.GetParameter("TdcKEKARE72KVC", 0);
  static const auto MaxTdcKVC = gUser.GetParameter("TdcKEKARE72KVC", 1);

  rawData->DecodeHits();

  gRM.Decode();

  // event.evnum = gRM.EventNumber();
  event.evnum = gUnpacker.get_event_number();
  event.spill = gRM.SpillNumber();

  HF1(1, 0);

  //**************************************************************************
  //****************** RawData

  // Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  for(const auto& hit: rawData->GetTrigRawHC()){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc1();
    if(tdc > 0){
      event.trigpat[trigger_flag.count()] = seg;
      event.trigflag[seg] = tdc;
      trigger_flag.set(seg);
    }
  }

  if(trigger_flag[trigger::kSpillEnd]) return true;

  HF1(1, 1);

  ///// T1
  {
    Int_t t1_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetKEKART1RawHC();
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t A = hit->GetAdcUp();
      event.t1a[seg-1] = A;
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.t1t[seg-1][m] = T;
        if(MinTdcT1 < T && T < MaxTdcT1) is_hit = true;
      }
      // Hitpat
      if(is_hit){
        event.t1hitpat[t1_nhits++] = seg;
      }
    }
    event.t1nhits = t1_nhits;
  }
  ///// T2
  {
    Int_t t2_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetKEKART2RawHC();
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t A = hit->GetAdcUp();
      event.t2a[seg-1] = A;
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.t2t[seg-1][m] = T;
        if(MinTdcT2 < T && T < MaxTdcT2) is_hit = true;
      }
      // Hitpat
      if(is_hit){
        event.t2hitpat[t2_nhits++] = seg;
      }
    }
    event.t2nhits = t2_nhits;
  }
  ///// T3
  {
    Int_t t3_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetKEKART3RawHC();
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t A = hit->GetAdcUp();
      event.t3a[seg-1] = A;
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.t3t[seg-1][m] = T;
        if(MinTdcT3 < T && T < MaxTdcT3) is_hit = true;
      }
      // Hitpat
      if(is_hit){
        event.t3hitpat[t3_nhits++] = seg;
      }
    }
    event.t3nhits = t3_nhits;
  }
  ///// T4
  {
    Int_t t4_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetKEKART4RawHC();
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t A = hit->GetAdcUp();
      event.t4a[seg-1] = A;
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.t4t[seg-1][m] = T;
        if(MinTdcT4 < T && T < MaxTdcT4) is_hit = true;
      }
      // Hitpat
      if(is_hit){
        event.t4hitpat[t4_nhits++] = seg;
      }
    }
    event.t4nhits = t4_nhits;
  }
  ///// SAC
  {
    Int_t sac_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetKEKARE90SACRawHC();
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t A = hit->GetAdcUp();
      event.saca[seg-1] = A;
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.sact[seg-1][m] = T;
        if(MinTdcSAC < T && T < MaxTdcSAC) is_hit = true;
      }
      // Hitpat
      if(is_hit){
        event.sachitpat[sac_nhits++] = seg;
      }
    }
    event.sacnhits = sac_nhits;
  }
  ///// SACSUM
  {
    Int_t sacsum_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetKEKARE90SACSUMRawHC();
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t A = hit->GetAdcUp();
      event.sacsuma[seg-1] = A;
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.sacsumt[seg-1][m] = T;
        if(MinTdcSAC < T && T < MaxTdcSAC) is_hit = true;
      }
      if(is_hit) sacsum_nhits++;
    }
    event.sacsumnhits = sacsum_nhits;
  }
  ///// BAC
  {
    Int_t bac_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetKEKARE72BACRawHC();
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t A = hit->GetAdcUp();
      event.baca[seg-1] = A;
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.bact[seg-1][m] = T;
        if(MinTdcBAC < T && T < MaxTdcBAC) is_hit = true;
      }
      // Hitpat
      if(is_hit){
        event.bachitpat[bac_nhits++] = seg;
      }
    }
    event.bacnhits = bac_nhits;
  }
  ///// BACSUM
  {
    Int_t bacsum_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetKEKARE72BACSUMRawHC();
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t A = hit->GetAdcUp();
      event.bacsuma[seg-1] = A;
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.bacsumt[seg-1][m] = T;
        if(MinTdcBAC < T && T < MaxTdcBAC) is_hit = true;
      }
      if(is_hit) bacsum_nhits++;
    }
    event.bacsumnhits = bacsum_nhits;
  }
  ///// KVC
  {
    Int_t kvc_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetKEKARE72KVCRawHC();
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t A = hit->GetAdcUp();
      event.kvca[seg-1] = A;
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.kvct[seg-1][m] = T;
        if(MinTdcKVC < T && T < MaxTdcKVC) is_hit = true;
      }
      // Hitpat
      if(is_hit){
        event.kvchitpat[kvc_nhits++] = seg;
      }
    }
    event.kvcnhits = kvc_nhits;
  }
  ///// KVCSUM
  {
    Int_t kvcsum_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetKEKARE72KVCSUMRawHC();
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t A = hit->GetAdcUp();
      event.kvcsuma[seg-1] = A;
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.kvcsumt[seg-1][m] = T;
        if(MinTdcKVC < T && T < MaxTdcKVC) is_hit = true;
      }
      // HitPat
      if(is_hit){
        event.kvcsumhitpat[kvcsum_nhits++] = seg;
      }
    }
    event.kvcsumnhits = kvcsum_nhits;
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
UserKEKAR::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserKEKAR;
}

//_____________________________________________________________________________
namespace
{
const Int_t    NbinAdc = 4096;
const Double_t MinAdc  =    0.;
const Double_t MaxAdc  = 4096.;

const Int_t    NbinTdc = 4096;
const Double_t MinTdc  =    0.;
const Double_t MaxTdc  = 4096.;

const Int_t    NbinTdcHr = 1e6/10;
const Double_t MinTdcHr  =  0.;
const Double_t MaxTdcHr  = 1e6;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  HB1(1, "Status", 20, 0., 20.);
  HB1(10, "Trigger HitPat", NumOfSegTrig, 0., Double_t(NumOfSegTrig));
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1(10+i+1, Form("Trigger Trig %d", i+1), 0x1000, 0, 0x1000);
  }

  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("spill",     &event.spill,     "spill/I");
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  //T1
  tree->Branch("t1nhits",   &event.t1nhits,    "t1nhits/I");
  tree->Branch("t1hitpat",   event.t1hitpat,   Form("t1hitpat[%d]/I",NumOfSegKEKART1));
  tree->Branch("t1a",        event.t1a,        Form("t1a[%d]/D",     NumOfSegKEKART1));
  tree->Branch("t1t",        event.t1t,        Form("t1t[%d][%d]/D", NumOfSegKEKART1, MaxDepth));
  //T2
  tree->Branch("t2nhits",   &event.t2nhits,    "t2nhits/I");
  tree->Branch("t2hitpat",   event.t2hitpat,   Form("t2hitpat[%d]/I",NumOfSegKEKART2));
  tree->Branch("t2a",        event.t2a,        Form("t2a[%d]/D",     NumOfSegKEKART2));
  tree->Branch("t2t",        event.t2t,        Form("t2t[%d][%d]/D", NumOfSegKEKART2, MaxDepth));
  //T3
  tree->Branch("t3nhits",   &event.t3nhits,    "t3nhits/I");
  tree->Branch("t3hitpat",   event.t3hitpat,   Form("t3hitpat[%d]/I",NumOfSegKEKART3));
  tree->Branch("t3a",        event.t3a,        Form("t3a[%d]/D",     NumOfSegKEKART3));
  tree->Branch("t3t",        event.t3t,        Form("t3t[%d][%d]/D", NumOfSegKEKART3, MaxDepth));
  //T4
  tree->Branch("t4nhits",   &event.t4nhits,    "t4nhits/I");
  tree->Branch("t4hitpat",   event.t4hitpat,   Form("t4hitpat[%d]/I",NumOfSegKEKART4));
  tree->Branch("t4a",        event.t4a,        Form("t4a[%d]/D",     NumOfSegKEKART4));
  tree->Branch("t4t",        event.t4t,        Form("t4t[%d][%d]/D", NumOfSegKEKART4, MaxDepth));

  //SAC
  tree->Branch("sacnhits",   &event.sacnhits,   "sacnhits/I");
  tree->Branch("sachitpat",   event.sachitpat,  Form("sachitpat[%d]/I", NumOfSegKEKARE90SAC));
  tree->Branch("saca",        event.saca,       Form("saca[%d]/D", NumOfSegKEKARE90SAC));
  tree->Branch("sact",        event.sact,       Form("sact[%d][%d]/D", NumOfSegKEKARE90SAC, MaxDepth));
  //SACSUM
  tree->Branch("sacsumnhits",   &event.sacsumnhits,   "sacsumnhits/I");
  tree->Branch("sacsuma",        event.sacsuma,       Form("sacsuma[%d]/D", NumOfSegKEKARE90SACSUM));
  tree->Branch("sacsumt",        event.sacsumt,       Form("sacsumt[%d][%d]/D", NumOfSegKEKARE90SACSUM, MaxDepth));
  //BAC
  tree->Branch("bacnhits",   &event.bacnhits,   "bacnhits/I");
  tree->Branch("bachitpat",   event.bachitpat,  Form("bachitpat[%d]/I", NumOfSegKEKARE72BAC));
  tree->Branch("baca",        event.baca,       Form("baca[%d]/D", NumOfSegKEKARE72BAC));
  tree->Branch("bact",        event.bact,       Form("bact[%d][%d]/D", NumOfSegKEKARE72BAC, MaxDepth));
  //BAC
  tree->Branch("bacsumnhits",   &event.bacsumnhits,   "bacsumnhits/I");
  tree->Branch("bacsuma",        event.bacsuma,       Form("bacsuma[%d]/D", NumOfSegKEKARE72BACSUM));
  tree->Branch("bacsumt",        event.bacsumt,       Form("bacsumt[%d][%d]/D", NumOfSegKEKARE72BACSUM, MaxDepth));
  //KVC
  tree->Branch("kvcnhits",  &event.kvcnhits,    "kvcnhits/I");
  tree->Branch("kvchitpat",  event.kvchitpat,  Form("kvchitpat[%d]/I", NumOfSegKEKARE72KVC));
  tree->Branch("kvca",       event.kvca,       Form("kvca[%d]/D", NumOfSegKEKARE72KVC));
  tree->Branch("kvct",       event.kvct,       Form("kvct[%d][%d]/D", NumOfSegKEKARE72KVC, MaxDepth));
  //KVCSUM
  tree->Branch("kvcsumnhits",   &event.kvcsumnhits,    "kvcsumnhits/I");
  tree->Branch("kvcsumhitpat",   event.kvcsumhitpat,   Form("kvcsumhitpat[%d]/I", NumOfSegKEKARE72KVCSUM));
  tree->Branch("kvcsuma",        event.kvcsuma,       Form("kvcsuma[%d]/D", NumOfSegKEKARE72KVCSUM));
  tree->Branch("kvcsumt",        event.kvcsumt,       Form("kvcsumt[%d][%d]/D", NumOfSegKEKARE72KVCSUM, MaxDepth));


  // HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO") &&
     InitializeParameter<HodoParamMan>("HDPRM") &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
