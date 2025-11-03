// -*- C++ -*-

#include "FiberHit.hh"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <limits>

#include <std_ostream.hh>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "FuncName.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "PrintHelper.hh"
#include "RawData.hh"

namespace
{
const auto qnan = TMath::QuietNaN();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gHodo = HodoParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
}

//_____________________________________________________________________________
FiberHit::FiberHit(HodoRawHit *rhit)
  : HodoHit(rhit),
    m_position(qnan),
    m_adccor_hi(m_n_ch),
    m_adccor_low(m_n_ch),
    m_mip_hi(m_n_ch),
    m_mip_low(m_n_ch)
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
FiberHit::~FiberHit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
bool
FiberHit::Calculate()
{
  if(!HodoHit::Calculate())
    ;
    //return false;

  m_is_clustered.clear();

  Int_t id    = m_raw->DetectorId();
  Int_t plane = m_raw->PlaneId();
  Int_t seg   = m_raw->SegmentId();

  data_t trailing(m_n_ch);
  for(Int_t ch=0; ch<m_n_ch; ++ch){
    const auto& l_cont = m_time_leading[ch];
    m_ctime_leading[ch].clear();
    m_ctime_trailing[ch].clear();
    for(Int_t il=0, nl=l_cont.size(); il<nl; ++il){
      Double_t l = l_cont[il];
      Double_t l_next = (il+1) != nl ? l_cont[il+1] : DBL_MAX;
      Double_t buf = qnan;
      for(const auto& t: m_time_trailing[ch]){
        if(l<t && t<l_next){
          buf = t;
          break;
        }
      }
      trailing[ch].push_back(buf);
      Double_t ctime = qnan;
      Double_t tot = buf - l;
      gPHC.DoCorrection(id, plane, seg, ch, l, tot, ctime);
      m_ctime_leading[ch].push_back(ctime);
      m_ctime_trailing[ch].push_back(ctime + tot); // no use
      m_is_clustered.push_back(false);
    }
  }
  m_time_trailing = trailing;

  if (id != 31) {
    Int_t layer = gGeom.GetDetectorId(DetectorName()+"-"+PlaneName());
    m_position = gGeom.CalcWirePosition(layer, seg);
    m_dxdw     = gGeom.dXdW(layer);
  }
  //return true;
  
#if 1
  if(id==31){// CFT ADC
    Int_t ch=0;
    m_adccor_hi[ch].clear();
    m_adccor_low[ch].clear();
    m_mip_hi[ch].clear();
    m_mip_low[ch].clear();                  

    Double_t nhit_adc = m_raw->GetSizeAdcHigh();
    if(nhit_adc>0){
      Double_t hi  =  m_raw->GetAdcHigh();
      Double_t low =  -9999.;
      if (m_raw->GetSizeAdcLow()>0)
	low =  m_raw->GetAdcLow();

      Double_t pedeHi  = gHodo.GetP0(id, plane, seg, 0);
      Double_t pedeLow = gHodo.GetP0(id, plane, seg, 1);
      Double_t gainHi  = gHodo.GetP1(id, plane, seg, 0);// pedestal+mip(or peak)
      Double_t gainLow = gHodo.GetP1(id, plane, seg, 1);// pedestal+mip(or peak)
      Double_t Alow = gHodo.GetP0(id, plane, 0, 3);// same value for the same layer
      Double_t Blow = gHodo.GetP1(id, plane, 0, 3);// same value for the same layer

      if (hi>0) {
	Double_t adccor_hi  = hi  - pedeHi;
	if (m_pedcor_hg>-2000 && hi >0)
	  adccor_hi  = hi  + m_pedcor_hg;

	m_adccor_hi[ch].push_back(adccor_hi);	

	Double_t mip_hi = adccor_hi/gainHi ;
	m_mip_hi[ch].push_back(mip_hi);	
      }

      if (low>0) {
	Double_t adccor_low = low - pedeLow;

	if (m_pedcor_lg>-2000 && low >0)
	  adccor_low  = low  + m_pedcor_lg;

	m_adccor_low[ch].push_back(adccor_low);		

	Double_t mip_low = adccor_low/gainLow;	
	m_mip_low[ch].push_back(mip_low);

	Double_t de_low = 0; // MeV
	if(mip_low>0){
	  de_low = -(Alow/Blow) * log(1. - mip_low/Alow);// [MeV]
	  if(1-mip_low/Alow<0){ // when pe is too big
	    de_low = -(Alow/Blow) * log(1. - (Alow-0.001)/Alow);// [MeV] almost max
	  }
	}
	m_de_low[ch].push_back(de_low);
      }

#if 0
      if(m_dE_lg > 10){
        std::cout << "layer = " << plane << ", seg = " << seg << ", adcLow = "
                  << m_adc_lg << ", dE = " << m_dE_lg << ", mip_lg = "
                  << m_mip_lg << ", gainLow = " << gainLow << ", gainHi = "
                  << gainHi << std::endl;
      }
#endif

#if 0      
      for(auto& pair: m_pair_cont){
	Double_t time= pair.time_l;
	Double_t ctime = -100;
	if (m_adc_hg>20) {
	  gPHC.DoCorrection(id, plane, seg, m_ud, time, m_adc_hg, ctime);
	  pair.ctime_l = ctime;
	} else
	  pair.ctime_l = time;
      }
#endif      
    }
  }
#endif
  m_is_calculated = true;

  //HodoHit::Print();
  return true;
}

//_____________________________________________________________________________
Double_t
FiberHit::MeanTimeOverThreshold(Int_t j) const
{
  try {
    if(m_n_ch == 1){
      return TimeOverThreshold(HodoRawHit::kUp, j);
    }else{
      return TMath::Sqrt(
        TMath::Abs(TimeOverThreshold(HodoRawHit::kUp, j) *
                   TimeOverThreshold(HodoRawHit::kDown, j)));
    }
  }catch(const std::out_of_range&){
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
void
FiberHit::Print(Option_t* arg) const
{
  PrintHelper helper(3, std::ios::fixed);
  hddaq::cout << FUNC_NAME << " " << arg << std::endl
	      << "detector_name = " << m_raw->DetectorName() << std::endl
	      << "detector_id   = " << m_raw->DetectorId() << std::endl
	      << "plane_name    = " << m_raw->PlaneName()  << std::endl
	      << "plane_id      = " << m_raw->PlaneId()    << std::endl
	      << "segment_id    = " << m_raw->SegmentId()  << std::endl
              << "n_ch          = " << m_n_ch              << std::endl
              << "de            = " << DeltaE() << std::endl
              << "mt/cmt        = " << MeanTime()
              << " / " << CMeanTime() << std::endl
              << "mtot          = " << MeanTOT() << std::endl
              << "tdiff/ctdiff  = " << TimeDiff()
              << " / " << CTimeDiff() << std::endl;
  for(const auto data_map: std::map<TString, data_t>
        {{"de-hi  ", m_de_high},      {"de-lo  ", m_de_low},
         {"time-l ", m_time_leading}, {"time-t ", m_time_trailing},
         {"ctime-l", m_ctime_leading}, {"ctime-t", m_ctime_trailing}
        }){
    for(const auto& cont: data_map.second){
      hddaq::cout << " " << data_map.first << ":" << cont.size()
                  << " ";
      std::copy(cont.begin(), cont.end(),
                std::ostream_iterator<Double_t>(hddaq::cout," "));
    }
    hddaq::cout << std::endl;
  }
  hddaq::cout << " tot    :";
  for(Int_t ch=0; ch<m_n_ch; ++ch){
    for(Int_t j=0, n=GetEntries(ch); j<n; ++j){
      hddaq::cout << " " << TOT(ch, j);
    }
  }
  hddaq::cout << std::endl;
}
