// -*- C++ -*-

#include "CFTFiberCluster.hh"
#include "CFTPosParamMan.hh"
#include "DetectorID.hh"

#include <cmath>
#include <string>
#include <limits>

#include <TMath.h>

#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "FiberHit.hh"
#include "FuncName.hh"

namespace
{
const Bool_t reject_nan = false;
const CFTPosParamMan& gCFTPos = CFTPosParamMan::GetInstance();
}

//_____________________________________________________________________________
CFTFiberCluster::CFTFiberCluster(const HodoHC& cont,
                           const index_t& index)
  : HodoCluster(cont, index),
    m_total_adc_hi(),
    m_max_adc_hi(),  
    m_total_adc_low(),
    m_max_adc_low(),  
    m_total_mip_hi(),
    m_max_mip_hi(),  
    m_total_mip_low(),
    m_max_mip_low(),  
    m_total_de_hi(),
    m_max_de_hi(),  
    m_total_de_low(),
    m_max_de_low(),  
    m_max_ctime(),
    m_max_segment(),    
    m_mean_x(),
    m_mean_y(),  
    m_mean_phi(),
    m_mean_r(),  
    m_mean_z0(),
    m_slope(),
    m_z(),
    m_mean_x_cor(),
    m_mean_y_cor(),
    m_mean_phi_cor(),
    m_x_cal(),
    m_y_cal(),
    m_z_cal(),
    m_phi_cal(),
    m_pathlength()
{
  Calculate();
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
CFTFiberCluster::~CFTFiberCluster()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
CFTFiberCluster::Calculate()
{
  if(!m_is_good)
    return;

  m_total_adc_hi = 0.;
  m_max_adc_hi = 0.;  
  m_total_adc_low = 0.;
  m_max_adc_low = 0.;

  m_max_segment = -1;
  
  m_total_mip_hi = 0.;
  m_max_mip_hi = 0.;  
  m_total_mip_low = 0.;
  m_max_mip_low = 0.;  
  m_total_de_hi = 0.;
  m_max_de_hi = 0.;  
  m_total_de_low = 0.;
  m_max_de_low = 0.;
  
  m_de = 0.;
  
  m_max_ctime = 0.;
  m_mean_x = TMath::QuietNaN();
  m_mean_y = TMath::QuietNaN();
  m_mean_phi = TMath::QuietNaN();
  m_mean_r = TMath::QuietNaN();
  m_mean_z0 = TMath::QuietNaN();
  m_slope = TMath::QuietNaN();
  m_z = 0.;
  m_mean_x_cor = 0.;
  m_mean_y_cor = 0.;
  m_mean_phi_cor = 0.;
  m_x_cal = 0.;
  m_y_cal = 0.;
  m_z_cal = 0.;
  m_phi_cal = 0.;
  m_pathlength = 0.;
  
  m_tot           = 0.;
  m_mean_position = 0.;

  for(Int_t i=0; i<m_cluster_size; ++i){
    const auto& hit = dynamic_cast<CFTFiberHit*>(m_hit_container[i]);
    const auto& index = m_index[i];
    m_tot = TMath::Max(m_tot, hit->MeanTOT(index));
    // m_tot += hit->TOT(index);
    m_mean_position += hit->Position();

    Double_t adccorHi  = hit->GetAdcCorHigh();
    if (adccorHi>0) {
      m_total_adc_hi += adccorHi;
      if (adccorHi > m_max_adc_hi)
	m_max_adc_hi = adccorHi;
    }
    Double_t adccorLow  = hit->GetAdcCorLow();
    Double_t mipLow     = hit->GetMipLow();    
    if (adccorLow>0) {
      m_total_adc_low += adccorLow;
      m_total_mip_low   += mipLow;
      if (adccorLow > m_max_adc_low) {
	m_max_adc_low = adccorLow;
	m_max_segment = hit->SegmentId();
	m_max_mip_low = mipLow;
      }
    }

    Double_t dE_Hi  = hit->GetMipHigh();
    if (dE_Hi>0) {
      m_total_de_hi += dE_Hi;
      if (dE_Hi > m_max_de_hi) {
	m_max_de_hi = dE_Hi;
      }
    }

    Double_t dE_Low  = hit->DeltaELowGain();
    if (dE_Low>0) {
      m_total_de_low += dE_Low;
      m_de += dE_Low;      
      if (dE_Low > m_max_de_low) {
	m_max_de_low = dE_Low;
	Double_t max_ctime = hit->CMeanTime(index);
	m_max_ctime = max_ctime;
      }
    }
  }
  
  m_mean_position /= Double_t(m_cluster_size);
  m_mean_time = m_max_ctime;
  m_ctime = m_max_ctime;  

  Int_t n_true_hit=0;
  for(Int_t i=0; i<m_cluster_size; ++i){
    const auto& hit = dynamic_cast<CFTFiberHit*>(m_hit_container[i]);
    Double_t dE_Low  = hit->DeltaELowGain();

    if (dE_Low > m_max_de_low*0.3) {
      n_true_hit++;
      if (n_true_hit == 1) {
	m_mean_x   = hit->GetX();
	m_mean_y   = hit->GetY();
	m_mean_phi = hit->GetPhi();
	m_mean_r   = hit->GetR();
	m_mean_z0  = hit->GetZ0();
	m_slope    = hit->GetSlope();
      } else {
	m_mean_x   += hit->GetX();
	m_mean_y   += hit->GetY();
	m_mean_phi += hit->GetPhi();
	m_mean_r   += hit->GetR();
	m_mean_z0  += hit->GetZ0();
      }
    }
  }
  if (n_true_hit>0) {
    m_mean_x   /= double(n_true_hit);
    m_mean_y   /= double(n_true_hit);
    m_mean_phi /= double(n_true_hit);
    m_mean_r   /= double(n_true_hit);
    m_mean_z0  /= double(n_true_hit);
  }
  
}


void CFTFiberCluster::SetCorPhi()
{
  static const std::string funcname = "[CFTFiberCluster::SetCorPhi]";
  
  if(!gCFTPos.IsReady()){
    std::cerr << funcname << ": cannot get CFTPosParamManager" << std::endl; 
    return;
  }

  
  Int_t phi_layer = -1; // 0-3 only phi
  if (PlaneName()=="PHI1")
    phi_layer = 0;
  else if (PlaneName()=="PHI2")
    phi_layer = 1;
  else if (PlaneName()=="PHI3")
    phi_layer = 2;
  else if (PlaneName()=="PHI4")
    phi_layer = 3;
  
  Int_t angle = ((Int_t)m_mean_phi/10);

  //std::cout << "layer : " << phi_layer << ", angle : " << angle
  //<< ", z_cal : " << m_z_cal << ", " << m_mean_phi << " --> ";

  m_mean_phi_cor = m_mean_phi - gCFTPos.GetPhiShift(phi_layer, angle, m_z_cal);
  //std::cout << m_mean_phi_cor << std::endl;
  m_mean_x_cor = m_mean_r * TMath::Cos(m_mean_phi_cor*TMath::DegToRad());
  m_mean_y_cor = m_mean_r * TMath::Sin(m_mean_phi_cor*TMath::DegToRad());

}

Double_t CFTFiberCluster::GetCorZ(Double_t phi, Double_t mean_z, Double_t theta)
{
  static const std::string funcname = "[CFTFiberCluster::GetCorZ]";

  if(!gCFTPos.IsReady()){
    std::cerr << funcname << ": cannot get CFTPosParamManager" << std::endl; 
    return 0;
  }

  Int_t uv_layer = -1; // 0-3 only uv
  if (PlaneName()=="UV1")
    uv_layer = 0;
  else if (PlaneName()=="UV2")
    uv_layer = 1;
  else if (PlaneName()=="UV3")
    uv_layer = 2;
  else if (PlaneName()=="UV4")
    uv_layer = 3;

  Double_t meanSeg = m_segment;
  Int_t    seg1 = (Int_t)meanSeg;
  Int_t    seg2 = seg1 + 1;
  Double_t ratio = meanSeg - (Double_t)seg1;

  Int_t    plane = PlaneId(); // 0-7 including uv, phi
  if (seg2 >= NumOfSegCFT[plane])
    seg2 = 0;


  Double_t z1 = gCFTPos.GetZposU(uv_layer, seg1, phi, theta);  
  Double_t z2 = gCFTPos.GetZposU(uv_layer, seg2, phi, theta);  

  /* phi range determination */
  Double_t phi_range1 = 90.; 
  Double_t phi_range2 = 450.; 

  Double_t theta_range_diff = -999;

  if (uv_layer==0 || uv_layer==2) {
    phi_range1 -= 360.*(Double_t)seg1/(Double_t)NumOfSegCFT[plane];
    phi_range2 -= 360.*(Double_t)seg1/(Double_t)NumOfSegCFT[plane];

    Double_t diff1 = std::abs(phi_range1 - phi);
    Double_t diff2 = std::abs(phi_range2 - phi);
    Double_t diff = diff1;
    if (diff2 < diff1)
      diff = diff2;

    theta_range_diff = diff;

    if (phi <= phi_range1)
      z1 = gCFTPos.GetZposU(uv_layer, seg1, phi+360, theta);      
    else if (phi >= phi_range2)
      z1 = gCFTPos.GetZposU(uv_layer, seg1, phi-360, theta);      
    /*
    std::cout << "Layer " << uv_layer 
	      << ", Phi_range " << phi_range1 << " -- " << phi_range2 
	      << ", phi" << phi 
	      << ", z1 = " << z1
	      << std::endl;
    */
    /*
    if (z1<0) {
      z1 = gCFTPos.GetZposU(layer, seg1, phi+360);      
    } else if (z1 > 400) {
      z1 = gCFTPos.GetZposU(layer, seg1, phi-360);      
    }
    */
  } else if (uv_layer==1 || uv_layer==3) {
    phi_range1 = -270. + 360*(Double_t)seg1/(Double_t)NumOfSegCFT[plane];
    phi_range2 = 90. + 360.*(Double_t)seg1/(Double_t)NumOfSegCFT[plane];

    Double_t diff1 = std::abs(phi_range1 - phi);
    Double_t diff2 = std::abs(phi_range2 - phi);
    Double_t diff = diff1;
    if (diff2 < diff1)
      diff = diff2;

    theta_range_diff = diff;

    if (phi <= phi_range1)
      z1 = gCFTPos.GetZposU(uv_layer, seg1, phi+360, theta);      
    else if (phi >= phi_range2)
      z1 = gCFTPos.GetZposU(uv_layer, seg1, phi-360, theta);      
    /*
    std::cout << "Layer " << uv_layer 
	      << ", Phi_range " << phi_range1 << " -- " << phi_range2 
	      << ", phi" << phi 
	      << ", z1 = " << z1
	      << std::endl;
    */
    /*
    if (z1<0) {
      z1 = gCFTPos.GetZposU(layer, seg1, phi-360);      
    } else if (z1 > 400) {
      z1 = gCFTPos.GetZposU(layer, seg1, phi+360);      
    }
    */
  }

  if (mean_z>-150 && mean_z <450 && std::abs(z1-mean_z)>200 
      && std::abs(theta_range_diff)<5) {
    //std::cout << "z1 = " << z1 << ", mean_z = " << mean_z ;
    if (z1-mean_z < 0)
      z1 += 400;
    else
      z1 -= 400;
    //std::cout << " --> z1 = " << z1 << std::endl;
  }

  if (uv_layer==0 || uv_layer==2) {
    phi_range1 = 90. - 360.*(Double_t)seg2/(Double_t)NumOfSegCFT[plane];
    phi_range2 = 450. - 360.*(Double_t)seg2/(Double_t)NumOfSegCFT[plane];
    
    Double_t diff1 = std::abs(phi_range1 - phi);
    Double_t diff2 = std::abs(phi_range2 - phi);
    Double_t diff = diff1;
    if (diff2 < diff1)
      diff = diff2;

    theta_range_diff = diff;

    if (phi <= phi_range1)
      z2 = gCFTPos.GetZposU(uv_layer, seg2, phi+360, theta);      
    else if (phi >= phi_range2)
      z2 = gCFTPos.GetZposU(uv_layer, seg2, phi-360, theta);      
    /*
    std::cout << "Layer " << uv_layer 
	      << ", Phi_range " << phi_range1 << " -- " << phi_range2 
	      << ", phi" << phi 
	      << ", z2 = " << z2
	      << std::endl;
    */

    /*
    if (z2<0) {
      z2 = gCFTPos.GetZposU(layer, seg2, phi+360);      
    } else if (z2 > 400) {
      z2 = gCFTPos.GetZposU(layer, seg2, phi-360);      
    }
    */
  } else if (uv_layer==1 || uv_layer==3) {
    phi_range1 = -270. + 360*(Double_t)seg2/(Double_t)NumOfSegCFT[plane];
    phi_range2 = 90. + 360.*(Double_t)seg2/(Double_t)NumOfSegCFT[plane];

    Double_t diff1 = std::abs(phi_range1 - phi);
    Double_t diff2 = std::abs(phi_range2 - phi);
    Double_t diff = diff1;
    if (diff2 < diff1)
      diff = diff2;

    theta_range_diff = diff;

    if (phi <= phi_range1)
      z2 = gCFTPos.GetZposU(uv_layer, seg2, phi+360, theta);      
    else if (phi >= phi_range2)
      z2 = gCFTPos.GetZposU(uv_layer, seg2, phi-360, theta);      

    /*
    std::cout << "Layer " << uv_layer 
	      << ", Phi_range " << phi_range1 << " -- " << phi_range2 
	      << ", phi" << phi 
	      << ", z2 = " << z2
	      << std::endl;
    */

    /*
    if (z2<0) {
      z2 = gCFTPos.GetZposU(layer, seg2, phi-360);      
    } else if (z2 > 400) {
      z2 = gCFTPos.GetZposU(layer, seg2, phi+360);      
    }
    */
  }

  if (mean_z>-300 && mean_z <300 && std::abs(z2-mean_z)>200
      && std::abs(theta_range_diff)<5) {
    //std::cout << "z2 = " << z2 << ", mean_z = " << mean_z ;
    if (z2-mean_z < 0)
      z2 += 400;
    else
      z2 -= 400;
    //std::cout << " --> z2 = " << z2 << std::endl;
  }

  return z1 + (z2-z1)*ratio;

}
