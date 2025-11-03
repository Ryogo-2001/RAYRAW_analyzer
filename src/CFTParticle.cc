#include "CFTParticle.hh"
#include "DCLTrackHit.hh"
#include "HodoWaveformHit.hh"
#include "FiberHit.hh"
#include <TRandom.h>
#include "HodoAnalyzer.hh"
#include "RawData.hh"
//#include "BGOPosCorMan.hh"
#include "CATCHPidMan.hh"
#include "UserParamMan.hh"
//#include "CFTdECorrMan.hh"

namespace
{
  //const BGOPosCorMan& gBGOPosMan = BGOPosCorMan::GetInstance();
  const CATCHPidMan&  gCATCHPidMan = CATCHPidMan::GetInstance();
  //const CFTdECorrMan& gCFTdECorMan = CFTdECorrMan::GetInstance();

  const auto& gUser = UserParamMan::GetInstance();
  const auto& U = HodoRawHit::kUp;  
}

//CFTParticle::CFTParticle(DCLocalTrack *track, RawData *rawData)
CFTParticle::CFTParticle(const CFTLocalTrack *track, HodoAnalyzer *hodoAna)
  : m_Track(track),// RawData_(rawData), 
    m_hodoAna(hodoAna),
    m_CFTVtx(-999, -999, -999),
    m_bgo_seg(-1),
    m_piid_seg(-1),
    m_bgo_energy(0),
    m_cft_Total_dE(0),
    m_cft_Total_dE_max(0),
    m_norm_Total_dE_max(0),
    m_TotalE(0),
    m_TotalE_max(0),
    m_Mass(-999.),
    m_dist_xy(0.),
    m_bgo_z_surface(-999.),
    m_bgo_adc(-999),
    m_bgo_adc_cor(-999)
{  
  FindBGO();
  //Calculate();
}

CFTParticle::~CFTParticle()
{
}

const HodoWaveformHit * CFTParticle::GetBGOHit(Int_t i)
{
  if (i>=0 && i<m_BGOCont.size())
    return m_BGOCont[i];
  else
    return 0;
}

HodoHit * CFTParticle::GetPiVHit(Int_t i)
{
  if (i>=0 && i<m_PiVCont.size())
    return m_PiVCont[i];
  else
    return 0;
}


Bool_t CFTParticle::FindBGO()
{
  static const auto MinTimeBGO = gUser.GetParameter("TimeBGO", 0);   // ns
  static const auto MaxTimeBGO = gUser.GetParameter("TimeBGO", 1);   // ns
  static const auto MinPulseTimeBGO = gUser.GetParameter("PulseTimeBGO", 0);   // ns
  static const auto MaxPulseTimeBGO = gUser.GetParameter("PulseTimeBGO", 1);   // ns

  Int_t xyFlag = m_Track->GetCFTxyFlag();
  Int_t zFlag  = m_Track->GetCFTzFlag() ;
  Double_t Axy = m_Track->GetAxy(); Double_t Bxy = m_Track->GetBxy();
  Double_t Az  = m_Track->GetAz() ; Double_t Bz  = m_Track->GetBz();
  ThreeVector Pos0 = m_Track->GetPos0();
  ThreeVector Dir  = m_Track->GetDir();

  // track to BGO segment
  Double_t distBGO=1000.;
  Double_t min_distBGO=1000.;
  Int_t    min_seg = -1;

  // BGO
  Int_t nhBGO = m_hodoAna->GetNHits("BGO");
  Double_t max_de = 0;
  for (Int_t i=0; i<nhBGO; i++) {
    const auto& hit = m_hodoAna->GetHit<HodoWaveformHit>("BGO", i);
    //Hodo1Hit* hit = m_hodoAna->GetHitBGO(i);
    Int_t seg = hit->SegmentId();
    Double_t de = -999.;
    Bool_t hit_flag = false;

    Int_t Npulse = hit->GetNPulse(U);
    for(Int_t m = 0; m<Npulse; ++m){
      Double_t pulse_time   = hit->GetPulseTime(m);
      Double_t pulse_de     = hit->DeltaE(m);		
      if(MinPulseTimeBGO<pulse_time&&pulse_time<MaxPulseTimeBGO){
	de = pulse_de;
      }
    }

    Int_t NhitT = hit->GetEntries(U);
    if(de>0){
      for(Int_t m=0; m<NhitT; ++m){
	Double_t time = hit->GetTUp(m);
	if(MinTimeBGO<time&&time<MaxTimeBGO){
	  hit_flag = true;
	}
      }
    }

    if(!hit_flag)continue;
    Double_t x, y;
    BGOPos(seg, &x, &y);

    if     (xyFlag==0){distBGO = (Axy*x-y+Bxy)/sqrt(Axy*Axy+1.*1.);}
    else if(xyFlag==1){distBGO = (Axy*y-x+Bxy)/sqrt(Axy*Axy+1.*1.);}
    Double_t u=Dir.x(), v=Dir.y();
    Double_t x0=Pos0.x(), y0=Pos0.y();
    Double_t t = (u*(x-x0)+v*(y-y0))/(u*u+v*v);
    if (t>=0) {	
      Double_t z = BGOZPosAtSurface(seg);
      if (fabs(distBGO)<25 && z >= -80 && z <= 330 ) {
	if (fabs(distBGO) <= fabs(min_distBGO) ) {
	  min_distBGO = distBGO;
	  min_seg = seg;
	}
      }
    }
  }


  if (min_seg >= 0 ) {
    for (Int_t i=0; i<nhBGO; i++) {
      const auto& hit = m_hodoAna->GetHit<HodoWaveformHit>("BGO", i);
      Int_t seg = hit->SegmentId();

      if (seg == min_seg) {
	  m_bgo_seg = seg;
	  hit->SetTrackJoined();
	  break;
      }
    }
  }

  return true;
}



Bool_t CFTParticle::Calculate()
{
  /* here!
  FiberTotal_E_ = m_Track->TotalDEHiGain();
  FiberMax_E_   = m_Track->MaxDEHiGain();
  m_PathLength   = m_Track->GetTotalPathLength();
  NormalizedFiberTotal_E_ = m_Track->NormalizedTotalDEHiGain();
  NormalizedFiberMax_E_   = m_Track->NormalizedMaxDEHiGain();
  */
  static const auto MinTimeBGO = gUser.GetParameter("TimeBGO", 0);   // ns
  static const auto MaxTimeBGO = gUser.GetParameter("TimeBGO", 1);   // ns
  static const auto MinPulseTimeBGO = gUser.GetParameter("PulseTimeBGO", 0);   // ns
  static const auto MaxPulseTimeBGO = gUser.GetParameter("PulseTimeBGO", 1);   // ns


  m_cft_Total_dE = m_Track->GetTotalSumdE();
  m_cft_Total_dE_max = m_Track->GetTotalMaxdE();

  Int_t nhit_phi = m_Track->GetNHit();
  Int_t nhit_uv  = m_Track->GetNHitUV();

  Double_t Total_dEphi_max = m_Track->GetTotalMaxdEphi();
  Double_t Total_dEuv_max  = m_Track->GetTotalMaxdEuv ();
  Double_t theta =m_Track->GetThetaCFT();

  m_norm_Total_dE_max = Total_dEphi_max * sin(theta*TMath::DegToRad())/nhit_phi 
    + Total_dEuv_max/nhit_uv;
    

  Int_t xyFlag = m_Track->GetCFTxyFlag();
  Int_t zFlag  = m_Track->GetCFTzFlag() ;
  Double_t Axy = m_Track->GetAxy(); Double_t Bxy = m_Track->GetBxy();
  Double_t Az  = m_Track->GetAz() ; Double_t Bz  = m_Track->GetBz();
  ThreeVector Pos0 = m_Track->GetPos0();
  ThreeVector Dir  = m_Track->GetDir();

  Double_t dist_xy = std::abs(-Dir.y()*Pos0.x() + Dir.x()*Pos0.y())/sqrt(Dir.x()*Dir.x() + Dir.y()*Dir.y());
  m_dist_xy = dist_xy;


  // track to BGO segment
  Double_t distBGO=1000.;

  // BGO
  Int_t nhBGO = m_hodoAna->GetNHits("BGO");
  Double_t max_de = 0;
  for (Int_t i=0; i<nhBGO; i++) {
    const auto& hit = m_hodoAna->GetHit<HodoWaveformHit>("BGO", i);
    Int_t seg = hit->SegmentId();
    
    if (seg == m_bgo_seg) {
      Double_t de = -999.;
      Double_t pulse_height = -999.;      
      Int_t Npulse = hit->GetNPulse(U);
      for(Int_t m = 0; m<Npulse; ++m){
	Double_t pulse_time   = hit->GetPulseTime(m);
	if(MinPulseTimeBGO<pulse_time&&pulse_time<MaxPulseTimeBGO){
	  de = hit->DeltaE(m);		
	  pulse_height = hit->GetPulseHeight(m);	  
	}
      }

      Bool_t hit_flag = false;
      Int_t NhitT = hit->GetEntries(U);
      if(de>0){
	for(Int_t m=0; m<NhitT; ++m){
	  Double_t time = hit->GetTUp(m);
	  if(MinTimeBGO<time&&time<MaxTimeBGO){
	    hit_flag = true;
	  }
	}
      }

      if(!hit_flag)continue;

      Double_t z = BGOZPosAtSurface(seg);

      Double_t adc = pulse_height;
      Double_t cor_adc = adc;
      /*
      if (z>0) {
	cor_adc = gBGOPosMan.CorAdc(seg, z, adc);
      }

      Double_t de_cor = -999.;
      if (gCalib.GetEnergy(seg, cor_adc, de_cor))
	m_bgo_energy += de_cor;
      else
	m_bgo_energy += de;
      */
      m_bgo_energy += de;      

      AddBGOHit(hit);
      max_de = de;
      m_bgo_z_surface = z;
      m_bgo_adc = adc;
      m_bgo_adc_cor = cor_adc;
    }
  }

  for (Int_t i=0; i<nhBGO; i++) {
    const auto& hit = m_hodoAna->GetHit<HodoWaveformHit>("BGO", i);
    Int_t seg = hit->SegmentId();
    Int_t Npulse = hit->GetNPulse(U);
      
    if (seg != m_bgo_seg && Npulse>0 && !hit->IsTrackJoined()) {
      Double_t de = -999.;
      Double_t pulse_height = -999.;      
      Int_t Npulse = hit->GetNPulse(U);
      for(Int_t m = 0; m<Npulse; ++m){
	Double_t pulse_time   = hit->GetPulseTime(m);
	if(MinPulseTimeBGO<pulse_time&&pulse_time<MaxPulseTimeBGO){
	  de = hit->DeltaE(m);		
	  pulse_height = hit->GetPulseHeight(m);	  
	}
      }

      Bool_t hit_flag = false;
      Int_t NhitT = hit->GetEntries(U);
      if(de>0){
	for(Int_t m=0; m<NhitT; ++m){
	  Double_t time = hit->GetTUp(m);
	  if(MinTimeBGO<time&&time<MaxTimeBGO){
	    hit_flag = true;
	  }
	}
      }

      if(!hit_flag)continue;
      Double_t x, y;
      BGOPos(seg, &x, &y);

      if     (xyFlag==0){distBGO = (Axy*x-y+Bxy)/sqrt(Axy*Axy+1.*1.);}
      else if(xyFlag==1){distBGO = (Axy*y-x+Bxy)/sqrt(Axy*Axy+1.*1.);}
      Double_t u=Dir.x(), v=Dir.y();
      Double_t x0=Pos0.x(), y0=Pos0.y();
      Double_t t = (u*(x-x0)+v*(y-y0))/(u*u+v*v);
      if (t>=0) {	
	Double_t z = BGOZPosAtSurface(seg);
	if (fabs(distBGO)<25 && z >= -80 && z <= 330) {
	  //Double_t adc = (Double_t)hit->GetRawHit()->GetAdc2();
	  Double_t adc = pulse_height;
	  Double_t cor_adc = adc;
	  //if (z>0) {
	  //cor_adc = gBGOPosMan.CorAdc(seg, z, adc);
	  //}

	  //Double_t de_cor = -999.;
	  //if (gCalib.GetEnergy(seg, cor_adc, de_cor))
	  //m_bgo_energy += de_cor;
	  //else
	  //m_bgo_energy += de;

	  m_bgo_energy += de;	  

	  AddBGOHit(hit);
	}
      }
    }
  }

#if 1  
  // track to PiID segment
  Double_t distPiID=1000.;
  Int_t nhPiID = m_hodoAna->GetNHits("PiID");
  for(Int_t i = 0; i<nhPiID; ++i){
    const auto& hit = m_hodoAna->GetHit<CFTFiberHit>("PiID", i);
    Int_t seg   = hit->SegmentId();
    Int_t NhitT = hit->GetEntries(U);
    Bool_t hit_flag = false;
    
    for(Int_t m = 0; m<NhitT; ++m){
      Double_t time = hit->GetTUp(m);
      if(time>-20&&time<50){hit_flag=true;}    
    }

    if(!hit_flag)continue;
    Double_t x, y;
    PiIDPos(seg, &x, &y);
    if     (xyFlag==0){distPiID = (Axy*x-y+Bxy)/sqrt(Axy*Axy+1.*1.);}
    else if(xyFlag==1){distPiID = (Axy*y-x+Bxy)/sqrt(Axy*Axy+1.*1.);}
    Double_t u=Dir.x(), v=Dir.y();
    Double_t x0=Pos0.x(), y0=Pos0.y();
    Double_t t = (u*(x-x0)+v*(y-y0))/(u*u+v*v);
    if (t>=0) {	
      if (fabs(distPiID)<40) {
	m_piid_seg = seg;
      }
    }
  }
#endif
  
#if 1    
  Int_t nc = m_BGOCont.size();
  for (Int_t i=0; i<nc; i++) {
    const HodoWaveformHit *hitp = m_BGOCont[i];
  }
  Int_t ncPiV = m_PiVCont.size();
  for (Int_t i=0; i<ncPiV; i++) {
    HodoHit *hitp = m_PiVCont[i];
    //PiV_E_ += hitp->DeltaE();
  }
#endif

  //TotalE_ = FiberTotal_E_ + BGO_E_;

#if 1    
  if (nc ==1 || nc>= 3) {
    Double_t delta;
    if (gCATCHPidMan.CheckProton(m_bgo_energy, m_norm_Total_dE_max, delta) && m_piid_seg < 0) {
      m_Mass = 0.9382720;
    } else if (gCATCHPidMan.CheckPi(m_bgo_energy, m_norm_Total_dE_max)) {
      m_Mass = 0.1395701;
    }
  } else if (nc == 2) {
    //Hodo2Hit *hitp1 = BGOCont_[0];
    const HodoWaveformHit *hitp1 = m_BGOCont[0];
    Double_t BGO_E1 = 0.;
    
    Int_t Npulse1 = hitp1->GetNPulse(U);
    for(Int_t m = 0; m<Npulse1; ++m){
      Double_t pulse_time   = hitp1->GetPulseTime(m);
      if(-50<pulse_time&&pulse_time<50){
	BGO_E1 = hitp1->DeltaE(m);		
      }
    }

    //Hodo2Hit *hitp2 = BGOCont_[1];
    const HodoWaveformHit *hitp2 = m_BGOCont[1];
    Double_t BGO_E2 = hitp2->DeltaE();
    Int_t Npulse2 = hitp2->GetNPulse(U);
    for(Int_t m = 0; m<Npulse2; ++m){
      Double_t pulse_time   = hitp2->GetPulseTime(m);
      if(-50<pulse_time&&pulse_time<50){
	BGO_E2 = hitp2->DeltaE(m);		
      }
    }
    
    Bool_t flag0 = false,  flag1 = false, flag2 = false;
    Double_t delta0, delta1, delta2;

    flag0 = gCATCHPidMan.CheckProton(m_bgo_energy, m_norm_Total_dE_max, delta0);
    flag1 = gCATCHPidMan.CheckProton(BGO_E1, m_norm_Total_dE_max, delta1);
    flag2 = gCATCHPidMan.CheckProton(BGO_E2, m_norm_Total_dE_max, delta2);

    if (flag0  && m_piid_seg < 0) {
      m_Mass = 0.9382720;      
    } else if (flag1 && flag2  && m_piid_seg < 0) {
      m_Mass = 0.9382720;      
      if (fabs(delta1) < fabs(delta2))
	m_bgo_energy = BGO_E1;
      else
	m_bgo_energy = BGO_E2;
    } else if (flag1  && m_piid_seg < 0) {
      m_Mass = 0.9382720;      
      m_bgo_energy = BGO_E1;
    } else if (flag2  && m_piid_seg < 0) {
      m_Mass = 0.9382720;      
      m_bgo_energy = BGO_E2;
      //}  else if (checkPi(m_bgo_energy)){
    } else if (gCATCHPidMan.CheckPi(m_bgo_energy, m_norm_Total_dE_max)) {
      m_Mass = 0.1395701;
    }
  } else if (nc == 0) {
    Double_t delta;
    if (gCATCHPidMan.CheckProton(m_bgo_energy, m_norm_Total_dE_max, delta) && m_piid_seg < 0)
      m_Mass = 0.9382720;
    else if (gCATCHPidMan.CheckPi(m_bgo_energy, m_norm_Total_dE_max)) 
      m_Mass = 0.1395701;

  }
#endif  

  Double_t dE_uv  = m_Track->GetTotalSumdEuv();
  Double_t dE_phi = m_Track->GetTotalSumdEphi();

  //Double_t cft_total_dE_cor = gCFTdECorMan.CorCFTdE(theta, m_cft_Total_dE + m_bgo_energy, dE_uv, dE_phi, m_cft_Total_dE );

  
  m_TotalE = m_cft_Total_dE + m_bgo_energy;
  //m_TotalE = cft_total_dE_cor + m_bgo_energy;

  if (m_TotalE<=0)
    m_TotalE = m_cft_Total_dE + m_bgo_energy;

  m_TotalE_max = m_cft_Total_dE_max + m_bgo_energy;

  return true;

}


void CFTParticle::BGOPos(Int_t seg, Double_t *x, Double_t *y) const
{
  Int_t UnitNum = seg/(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);
  Int_t SegInUnit = seg%(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);

  Double_t theta = 22.5+(Double_t)UnitNum*45.;
  Double_t x0 = RadiusOfBGOSurface+BGO_Y/2;
  Double_t y0 = (Double_t)(SegInUnit-1)*BGO_X;
  Double_t xc = 0.;
  Double_t yc = 0.;

  if(seg%3==2){ // single BGO
    Double_t n=(seg+1)/3;
    Double_t angle = +22.5+45.*(n-1); // axis change
    xc = (120.+25./2.)*cos(angle*TMath::DegToRad());
    yc = (120.+25./2.)*sin(angle*TMath::DegToRad());    
#if 1 // new
  }else if(seg==0 || seg==1){
    xc = 100.0 + 25./2.;
    if(seg==0){yc = -30.0/2.;}
    else if(seg==1){yc = 30.0/2.;}
  }else if(seg==6 || seg==7){
    yc = 100.0 + 25./2.;
    if(seg==6){xc = 30.0/2.;}
    else if(seg==7){xc = -30.0/2.;}
  }else if(seg==12 || seg==13){
    xc = -100.0 - 25./2.;
    if(seg==12){yc = 30.0/2.;}
    else if(seg==13){yc = -30.0/2.;}  
  }else if(seg==18 || seg==19){
    yc = -100.0 -25./2.;
    if     (seg==18){xc = -30.0/2.;}
    else if(seg==19){xc = 30.0/2.;}
  }else if(seg==3 || seg==4){
    Double_t angle = 45.;
    x0 = (100. + 25./2.)*cos(angle*TMath::DegToRad());
    y0 = (100. + 25./2.)*sin(angle*TMath::DegToRad());
    if(seg==4){
      xc = x0 - 30./2.*cos(-angle*TMath::DegToRad());
      yc = y0 - 30./2.*sin(-angle*TMath::DegToRad());	
    }else if(seg==3){
      xc = x0 + 30./2.*cos(-angle*TMath::DegToRad());
      yc = y0 + 30./2.*sin(-angle*TMath::DegToRad());	
    }
  }else if(seg==9 || seg==10){
    //Double_t angle = -45.;
    Double_t angle = 135.;
    x0 = (100 + 25./2.)*cos(angle*TMath::DegToRad());
    y0 = (100 + 25./2.)*sin(angle*TMath::DegToRad());      
    if(seg==10){
      xc = x0 + 30./2.*cos(-angle*TMath::DegToRad());
      yc = y0 + 30./2.*sin(-angle*TMath::DegToRad());	
    }else if(seg==9){
      xc = x0 - 30./2.*cos(-angle*TMath::DegToRad());
      yc = y0 - 30./2.*sin(-angle*TMath::DegToRad());	
    }
  }else if(seg==15 || seg==16){
    Double_t angle = -135.;
    x0 = (100. + 25./2.)*cos(angle*TMath::DegToRad());
    y0 = (100. + 25./2.)*sin(angle*TMath::DegToRad());      
    if(seg==16){
      xc = x0 - 30./2.*cos(-angle*TMath::DegToRad());
      yc = y0 - 30./2.*sin(-angle*TMath::DegToRad());	
    }else if(seg==15){
      xc = x0 + 30./2.*cos(-angle*TMath::DegToRad());
      yc = y0 + 30./2.*sin(-angle*TMath::DegToRad());	
    }
  }else if(seg==21 || seg==22){
    //Double_t angle = 135.;
    Double_t angle = -45.;
    x0 = (100. + 25./2.)*cos(angle*TMath::DegToRad());
    y0 = (100. + 25./2.)*sin(angle*TMath::DegToRad());      
    if(seg==22){
      xc = x0 + 30./2.*cos(-angle*TMath::DegToRad());
      yc = y0 + 30./2.*sin(-angle*TMath::DegToRad());	
    }else if(seg==21){
      xc = x0 - 30./2.*cos(-angle*TMath::DegToRad());
      yc = y0 - 30./2.*sin(-angle*TMath::DegToRad());	
    }
  }
#endif
  *x = xc;
  *y = yc;
}

Double_t CFTParticle::BGOZPosAtSurface(Int_t seg) const
{
  Int_t UnitNum = seg/(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);
  Int_t SegInUnit = seg%(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);

  Double_t theta = 22.5+(Double_t)UnitNum*45.;
  Double_t x1 = 0.;
  Double_t y1 = 0.;
  Double_t angle = 0.;

  if(seg%3==2){ // single BGO
    Double_t n=(seg+1)/3;
    angle = +22.5+45.*(n-1); // axis change
    x1 = RadiusOfBGOSurface2*cos(angle*TMath::DegToRad());
    y1 = RadiusOfBGOSurface2*sin(angle*TMath::DegToRad());    
  } else {
    angle = 45.*UnitNum ; // axis change
    x1 = RadiusOfBGOSurface*cos(angle*TMath::DegToRad());
    y1 = RadiusOfBGOSurface*sin(angle*TMath::DegToRad());    
  }

  ThreeVector Pos0 = m_Track->GetPos0();
  ThreeVector Dir  = m_Track->GetDir();

  Double_t x0 = Pos0.x();
  Double_t y0 = Pos0.y();
  Double_t z0 = Pos0.z();

  Double_t u = Dir.x();
  Double_t v = Dir.y();
  Double_t w = Dir.z();

  Double_t t = (-cos(angle*TMath::DegToRad())*(x0-x1)-sin(angle*TMath::DegToRad())*(y0-y1))/(u*cos(angle*TMath::DegToRad())+v*sin(angle*TMath::DegToRad()));

  Double_t z = z0 + w * t;

  return z;
}

void CFTParticle::PiIDPos(Int_t seg, Double_t *x, Double_t *y) const
{
  Int_t UnitNum = seg/(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);
  Int_t SegInUnit = seg%(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);
  Double_t theta = 22.5+(Double_t)UnitNum*45.;
  Double_t x0 = RadiusOfBGOSurface+BGO_Y/2;
  Double_t y0 = (Double_t)(SegInUnit-1)*BGO_X;
  Double_t xc = 0.;
  Double_t yc = 0.;
  Double_t w  = 30.,  t = 15.;//width, thickness
  Double_t ww = 40., tt = 15.;//width, thickness for 45 deg.

  if(seg%4==3){ // single segment
    Double_t n=(seg+1)/4;
    Double_t angle = +22.5+45.*(n-1); // axis change
    xc = (164.1+tt/2.)*cos(angle*TMath::DegToRad());
    yc = (164.1+tt/2.)*sin(angle*TMath::DegToRad());    
  }else if(seg==0 || seg==1 || seg==2){      
    xc = 159.0 + t/2.;
    if     (seg==0){yc = -1.*w;}
    else if(seg==1){yc =  0.*w;}
    else if(seg==2){yc =  1.*w;}      

  }else if(seg==8 || seg==9 || seg==10){
    yc = 159.0 + t/2.;
    if     (seg==8) {xc = 1.*w;}
    else if(seg==9) {xc = 0.*w;}
    else if(seg==10){xc =-1.*w;}      

  }else if(seg==16 || seg==17 || seg==18){
    xc = -159.0 - t/2.;
    if     (seg==16){yc = 1.*w;}
    else if(seg==17){yc = 0.*w;}
    else if(seg==18){yc =-1.*w;}      

  }else if(seg==24 || seg==25 || seg==26){
    yc = -159.0 - t/2.;
    if     (seg==24) {xc =-1.*w;}
    else if(seg==25) {xc = 0.*w;}
    else if(seg==26) {xc = 1.*w;}      

  }else if(seg==4 || seg==5 || seg==6){ // Line
    Double_t angle = 45.;
    x0 = (159.+t/2.)*cos(angle*TMath::DegToRad());
    y0 = (159.+t/2.)*sin(angle*TMath::DegToRad());
    if(seg==6){
      xc = x0 - w*cos(45.*TMath::DegToRad());
      yc = y0 + w*cos(45.*TMath::DegToRad());
    }if(seg==5){
      xc = x0 - 0*cos(45.*TMath::DegToRad());
      yc = y0 + 0*cos(45.*TMath::DegToRad());
    }else if(seg==4){
      xc = x0 + w*cos(45.*TMath::DegToRad()); 
      yc = y0 - w*cos(45.*TMath::DegToRad()); 
    }   

  }else if(seg==12 || seg==13 || seg==14){ // Line
      Double_t angle = 135.;
      x0 = (159.+t/2.)*cos(angle*TMath::DegToRad());
      y0 = (159.+t/2.)*sin(angle*TMath::DegToRad());
      if(seg==12){
	xc = x0 + w*cos(45.*TMath::DegToRad());
	yc = y0 + w*cos(45.*TMath::DegToRad());
      }else if(seg==13){
	xc = x0 + 0*cos(45.*TMath::DegToRad());
	yc = y0 + 0*cos(45.*TMath::DegToRad());
      }else if(seg==14){
	xc = x0 - w*cos(45.*TMath::DegToRad());
	yc = y0 - w*cos(45.*TMath::DegToRad());
      }

  }else if(seg==20 || seg==21 || seg==22){ // Line
    Double_t angle = -135.;
    x0 = (159. +t/2.)*cos(angle*TMath::DegToRad());
    y0 = (159. +t/2.)*sin(angle*TMath::DegToRad());
    if(seg==22){
      xc = x0 + 1.*w*cos(45.*TMath::DegToRad()); 
      yc = y0 - 1.*w*cos(45.*TMath::DegToRad()); 
    }else if(seg==21){
      xc = x0 + 0.*w*cos(45.*TMath::DegToRad()); 
      yc = y0 - 0.*w*cos(45.*TMath::DegToRad()); 
    }else if(seg==20){
      xc = x0 - 1.*w*cos(45.*TMath::DegToRad());
      yc = y0 + 1.*w*cos(45.*TMath::DegToRad());
    }

  }else if(seg==28 || seg==29 || seg==30){ // Line
    Double_t angle = -45.;
    x0 = (159.+t/2.)*cos(angle*TMath::DegToRad());
    y0 = (159.+t/2.)*sin(angle*TMath::DegToRad());
    if(seg==28){
      xc = x0 - 1.*w*cos(45.*TMath::DegToRad());
      yc = y0 - 1.*w*cos(45.*TMath::DegToRad());
    }else if(seg==29){
      xc = x0 - 0.*w*cos(45.*TMath::DegToRad());
      yc = y0 - 0.*w*cos(45.*TMath::DegToRad());
    }else if(seg==30){
      xc = x0 + 1.*w*cos(45.*TMath::DegToRad());
      yc = y0 + 1.*w*cos(45.*TMath::DegToRad());
    }
  }

  *x = xc;
  *y = yc;
}
