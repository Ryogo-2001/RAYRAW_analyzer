// -*- C++ -*-

#include "HodoWaveformHit.hh"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <limits>
#include <utility>

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
#include "BGODiscriminator.hh"
#include "TemplateFitMan.hh"
#include "BGOCalibMan.hh"

namespace
{
const auto qnan = TMath::QuietNaN();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gHodo = HodoParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
const auto& gTempFit = TemplateFitMan::GetInstance();
const auto& gBGOCalib = BGOCalibMan::GetInstance();

const double graphStart = -3.0;
const double graphEnd   = 2.0;
const double y_err      = 30.0;
const double y_err_TAG  = 3.0;
// org
const double fitStart = -1.0;
const double fitEnd   = 1.0;
//const double fitStart = -0.2;
//const double fitEnd   = 0.3;
const int ParaMax = 64;
//const double TrigTime = 3.70;
const double TrigTimeReso = 1.00;
}

//_____________________________________________________________________________
HodoWaveformHit::HodoWaveformHit(HodoRawHit *rhit)
  : HodoHit(rhit),
    m_waveform(m_n_ch),
    m_pulse_height(m_n_ch),
    m_pulse_time(m_n_ch),
    m_position(qnan),
    m_adc_integral(qnan),
    m_n_discri_pulse(qnan),
    m_n_discri_diffpulse(qnan),
    m_func(0),
    m_JoinTrack(false)
{
  debug::ObjectCounter::increase(ClassName());

  if (DetectorName()=="BGO") {
    auto &cont = gTempFit.GetHitContainer(DetectorName());
    TemplateFitFunction *tempFunc = cont.at(SegmentId());
    m_func = new TF1(DetectorName()+"-"+std::to_string(SegmentId()),
		     //gTempFit.GetFitFunction(SegmentId()),
		     tempFunc,
		     fitStart, fitEnd, ParaMax );
  }else if (DetectorName()=="TAG-PL") {
    auto &cont = gTempFit.GetHitContainer(DetectorName());
    TemplateFitFunction *tempFunc = cont.at(SegmentId());
    m_func = new TF1(DetectorName()+"-"+std::to_string(SegmentId()),
		     //gTempFit.GetFitFunction(SegmentId()),
		     tempFunc,
		     fitStart, fitEnd, ParaMax );
  }

}

//_____________________________________________________________________________
HodoWaveformHit::~HodoWaveformHit()
{
  del::ClearContainer( m_TGraphC );

  if (m_func)
    delete m_func;

  m_func = 0;
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
bool
HodoWaveformHit::Calculate()
{
  if(!HodoHit::Calculate())
    ;
    //return false;

  m_is_clustered.clear();

  for(Int_t ch=0; ch<m_n_ch; ++ch){
    //reset m_de_high
    m_de_high.at(ch).clear();
  }


  Int_t id    = m_raw->DetectorId();
  Int_t plane = m_raw->PlaneId();
  Int_t seg   = m_raw->SegmentId();

  for(Int_t ch=0; ch<m_n_ch; ++ch){
    // adc

    Int_t ns=0;
    for(const auto& adc: m_raw->GetArrayAdcHigh(ch)){
      Double_t de   = TMath::QuietNaN();
      Double_t time = TMath::QuietNaN();
      if(adc > 0 && adc < 0xffff && ns > 0 &&
	 gHodo.GetTime2(id, plane, seg, ch, ns, time)
	 ){
	Double_t pedestal  = gHodo.GetP0(id, plane, seg, 0);
	de = (Double_t)adc - pedestal;

	std::pair<Double_t, Double_t> wf_pair(time, de);
        m_waveform.at(ch).push_back(wf_pair);

      }
      ns++;
    }

    Int_t Nped=0;
    Double_t event_pedestal=0;
    for(const auto& wf: m_waveform.at(ch)){
      Double_t time = wf.first;
      if (time >= -0.06 && time < -0.015) {
	event_pedestal = wf.second;
	Nped++;
      }
    }
    event_pedestal /= Nped;
    if (Nped>0)
      m_adc_integral = 0.;

    for(const auto& wf: m_waveform.at(ch)){
      Double_t time = wf.first;
      if (time >= -0.01 && time <=0.03) {
	Double_t de = wf.second;
	m_adc_integral += (de - event_pedestal) * (-1.0);
      }
    }

  }

  PulseSearch();

  m_is_calculated = true;

  //HodoHit::Print();
  return true;
}

//_____________________________________________________________________________
Bool_t HodoWaveformHit::PulseSearch( void )
{
  if (!MakeGraph())
    return false;

  // Original graph index = 0;
  Int_t index_original_graph = 0;

  MakeDifGraph(index_original_graph);
  Int_t index_diff_graph = 1;

  Double_t threshold = -100;
  Double_t width = 0.05;
  Double_t risetime = 0.1;
  if(DetectorName()=="TAG-PL") {
    TGraphErrors *gr = m_TGraphC[index_original_graph];
    Int_t n = gr->GetN() - 1;
    Double_t *refy = gr->GetY();
    Double_t baseline=0;
    Int_t Nbaseline=0;
    for (Int_t i=0; i<n; i++) {
      if (refy[i+1] - refy[i] >-5 && refy[i+1] - refy[i] < 0) {
	baseline += refy[i+1];
	Nbaseline++;
      }
    }
    baseline /= Nbaseline;
    //std::cout << "baseline = " << baseline << std::endl;
    threshold = -30 + baseline;
    width = 0.05;
    risetime = 0.012;
  }

  SearchParam sp1={"sp1", {index_original_graph, index_diff_graph},
    fitStart, fitEnd, fitStart, fitEnd,
    threshold, width, risetime};

  Bool_t flagPresearch = PreSearch(&sp1);
  if (!flagPresearch)
    return false;

  // if (DetectorName()!="BGO")
  //   return true;

  Int_t color = 4;
  FitParam fp1={"fp1", index_original_graph, color ,fitStart, fitEnd};

  Bool_t flagSetFitParam = SetFitParam(&fp1,sp1.foundx,sp1.foundy);
  if (!flagSetFitParam)
    return false;

  Fit1(&fp1);

  //std::cout << "BGO Fitting result ; segment " << SegmentId() << std::endl;
  //for (int i=0; i<fp1.ParaNum; i++) {
  //std::cout << "p" << i << " : " << fp1.FitParam[i] << std::endl;
  //}

  Double_t trigx = FittedTrigX(fp1,1.0);
  //trigx =1;/////////////////////
  if(std::abs(trigx)<TrigTimeReso){
    Double_t max_res = 0;
    fp1.Residual=  RisingResidual(index_original_graph, trigx, max_res);
    //fp1.Residual=1;/////////////

    //m_FitParamCont.push_back(fp1);

    if(fp1.Residual <5 || fabs(max_res)<50){
      //if(fp1.Residual <5000 || fabs(max_res)<500){
      Int_t waveNum = fp1.wavenum;
      for (Int_t nw=0; nw<waveNum; nw++) {
	Double_t time = fp1.FitParam[2*nw+2];
	Double_t height = fp1.FitParam[2*nw+3];
        m_pulse_time.at(HodoRawHit::kUp).push_back(time);
        m_pulse_height.at(HodoRawHit::kUp).push_back(height);
	Double_t energy=-999.;
	if (GetName() == "BGO") {
	  gBGOCalib.GetEnergy(SegmentId(), height, energy);
	  m_de_high.at(HodoRawHit::kUp).push_back(energy);
	}
      }
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t HodoWaveformHit::MakeGraph()
{
  Int_t nc = GetWaveformEntries(HodoRawHit::kUp);
  Int_t n_range = 0;
  for (Int_t i=0; i<nc; i++) {
    std::pair<Double_t, Double_t> fadc = GetWaveform(HodoRawHit::kUp, i);
    if ( fadc.first >= graphStart && fadc.first <= graphEnd )
      n_range++;
  }

  if (n_range<=0)
    return false;

  TGraphErrors *gr = new TGraphErrors(n_range);
  Int_t index = 0;
  for (Int_t i=0; i<nc; i++) {
    std::pair<Double_t, Double_t> fadc = GetWaveform(i);
    if ( fadc.first >= graphStart && fadc.first <= graphEnd ) {
      if ( index<n_range ) {
	gr->SetPoint(index, fadc.first, fadc.second);
	if (GetName() == "BGO") {
	  gr->SetPointError(index, 0, y_err);
	}else if(GetName() == "TAG-PL"){
	  gr->SetPointError(index, 0, y_err_TAG);
	} else {
	  gr->SetPointError(index, 0, 0);
	}
	index++;
      }
    }
  }

  m_TGraphC.push_back(gr);

  return true;

}

//_____________________________________________________________________________
Bool_t HodoWaveformHit::MakeDifGraph(Int_t index_org)
{
  TGraphErrors *gr = m_TGraphC[index_org];

  Int_t n = gr->GetN() - 1;

  Double_t *refx = gr->GetX();
  Double_t *refy = gr->GetY();

  Double_t x[n], y[n];

  for (Int_t i=0; i<n; i++) {
    x[i] = refx[i];
    y[i] = refy[i+1] - refy[i];
  }

  TGraphErrors *gr_diff = new TGraphErrors(n, x, y);

  m_TGraphC.push_back(gr_diff);

  return true;

}

//_____________________________________________________________________________
Bool_t HodoWaveformHit::PreSearch(struct SearchParam *sp)
{
  static const std::string func_name("["+ClassName()+"::"+__func__+"()]");

  Double_t begin = sp ->sbegin;
  Double_t end   = sp ->send;

  //BGODiscriminator BGODiscri();
  BGODiscriminator BGODiscri;
  BGODiscri.SetGraph(m_TGraphC[sp->tgen[0]]);

  std::vector<Double_t> OriR,OriF;
  Double_t threshold = -500.;
  if(DetectorName()=="TAG-PL") threshold = -20.;
  Bool_t   flagSelectRange = BGODiscri.SelectRange(threshold, begin, end);

  if (!flagSelectRange)
    return false;

  BGODiscri.GetRisingPoint(OriR);
  BGODiscri.GetFallingPoint(OriF);
  BGODiscri.AllClear();

  m_n_discri_pulse = OriR.size();

  //bool FF = FlagFrontWave(OriR,begin);
  Bool_t flagFrontWave = false;
  if (OriR.size()>0)
    if (OriR[0] - begin < 0.1)
      flagFrontWave = true;

  BGODiscri.SetGraph(m_TGraphC[sp->tgen[1]]);

  std::vector<Double_t> rise30,fall30,entry30;
  threshold = -30;
  if(DetectorName()=="TAG-PL") threshold = -15.;

  flagSelectRange = BGODiscri.SelectRange(threshold, begin, end);
  if (!flagSelectRange)
    return false;

  BGODiscri.GetRisingPoint(rise30);
  BGODiscri.GetFallingPoint(fall30);

  // if(DetectorName()=="TAG-PL") {
  //   std::cout << "segment (" << m_raw->SegmentId() << ") Vth=-30, ";
  //   for(int i=0;i<rise30.size();i++)
  //     std::cout<<rise30[i] << ", ";
  //   std::cout << std::endl;
  // }

  Double_t width_thr = 0.04;
  if(DetectorName()=="TAG-PL") width_thr = 0.001;

  Bool_t   flagWidthCut1 = WidthCut(rise30,fall30,width_thr,entry30);
  if (!flagWidthCut1) {
    std::cout << func_name << " : WidthCut1, Number of rise and fall points does not match" << std::endl;
    return false;
  }

  ///Bool_t FF2= FlagFrontWave(rise30,begin);
  Bool_t flagFrontWave2 = false;
  if (rise30.size()>0)
    if (rise30[0] - begin < 0.1)
      flagFrontWave2 = true;


  BGODiscri.Clear();

  Int_t flagAroundBegin = 0;

  if(flagFrontWave && flagFrontWave2)
    flagAroundBegin = 1;
  if(flagFrontWave && !flagFrontWave2)
    flagAroundBegin = 2;
  if(!flagFrontWave && flagFrontWave2){
    //    std::cout<<"Event:"<<CEI->GetEventNum()<<" ch:"<<CEI->GetChannel()<<"Search1 ahoooooo"<<std::endl;
  }
  //  std::cout<<"Event FF:"<<FF<<std::endl;

  std::vector<Double_t> rise100,fall100,entry100,entry30_100;
  threshold = -100;

  flagSelectRange = BGODiscri.SelectRange(threshold , begin, end);
  if (!flagSelectRange)
    return false;

  BGODiscri.GetRisingPoint(rise100);
  BGODiscri.GetFallingPoint(fall100);
  // if(DetectorName()=="TAG-PL") {
  //   std::cout << "segment (" << m_raw->SegmentId() << ") Vth=-100, ";
  //   for(int i=0;i<rise100.size();i++)
  //     std::cout<<rise100[i] << ", ";
  //   std::cout << std::endl;
  // }


  width_thr = 0.08;
  if(DetectorName()=="TAG-PL") width_thr = 0.001;
  WidthCut(rise100, fall100, width_thr, entry100);

  CompareRise(entry30,entry100,0.05,entry30_100);
  // if(DetectorName()=="TAG-PL") {
  //   std::cout << "segment (" << m_raw->SegmentId() << "), ";
  //   for(int i=0;i<entry30_100.size();i++)
  //     std::cout<<entry30_100[i] << ", ";
  //   std::cout << std::endl;
  // }
  m_n_discri_diffpulse = entry30_100.size();

  if (entry30_100.size()==0)
    return false;

  SetInitial(entry30_100,begin,end,sp->threshold,sp->risetime);

  Int_t index_original_graph = 0;
  for(unsigned int i=0;i<entry30_100.size();i++){
    if (-GXtoGY(index_original_graph, entry30_100[i]+sp->risetime)>1) {
      sp->foundx.push_back(entry30_100[i]+sp->risetime);
      sp->foundy.push_back(-GXtoGY(index_original_graph,
				   entry30_100[i]+sp->risetime)) ;
    }
    // if(DetectorName()=="TAG-PL") {
    //   std::cout << "( " << entry30_100[i]+sp->risetime << ", "
    // 		<< -GXtoGY(index_original_graph, entry30_100[i]+sp->risetime) << "), ";
    //   if (i==entry30_100.size()-1)
    // 	std::cout << std::endl;
    // }
  }

  if(flagAroundBegin==1){
    sp->foundx.insert(sp->foundx.begin(),begin);
    sp->foundy.insert(sp->foundy.begin(),
		      -GXtoGY(index_original_graph, begin+0.03));
  }
  if(flagAroundBegin==2){
    sp->foundx.insert(sp->foundx.begin(),begin-0.5);
    sp->foundy.insert(sp->foundy.begin(),
		      -GXtoGY(index_original_graph, begin+0.03));
  }

  return true;

}

//_____________________________________________________________________________
Bool_t HodoWaveformHit::WidthCut(std::vector<Double_t> rise,
				 std::vector<Double_t> fall,
				 Double_t width, std::vector<Double_t> &outrise)
{

  if(rise.size() != fall.size()){
    std::cout<<"CAlgorithm::WidthCut rise num != fall num"<<std::endl;
    return false;
  }

  for(unsigned int i=0;i<rise.size();i++){
    if( fall[i]-rise[i] > width){
      outrise.push_back(rise[i]);
    }
  }
  rise.clear();
  fall.clear();

  return true;
}

//_____________________________________________________________________________
void HodoWaveformHit::CompareRise(std::vector<Double_t> rise1,
				  std::vector<Double_t> rise2,
				  Double_t width, std::vector<Double_t> &outrise)
{
  for(Int_t i=0;i<rise2.size();i++){
    Int_t t=0;
    for(Int_t j=0;j<rise1.size();j++){
      if( rise2[i]-rise1[j] <width )
        t++;
    }
    if(t==0)
      outrise.push_back(rise2[i]);
  }

  for(unsigned int j=0;j<rise1.size();j++)
    outrise.push_back(rise1[j]);

  rise1.clear();
  rise2.clear();
}

//_____________________________________________________________________________
void HodoWaveformHit::SetInitial(std::vector<Double_t> &v,
				 Double_t begin, Double_t end,
				 Double_t thre, Double_t rise)
{
  Int_t size=v.size();
  Double_t SepaLimit = 0.08;
  if(DetectorName()=="TAG-PL") {
    SepaLimit = 0.004;
  }
  if(size>1)
    for(Int_t i=0;i<size;i++){
      if(v[i]<begin || v[i]>end)
        v[i]=-1;
    }

  if(size>1){
    std::sort(v.begin(),v.end());

    for(Int_t i=0;i<size-1;i++){
      if(v[i+1]-v[i]<SepaLimit){
	v[i+1] = (Double_t)((v[i+1]+v[i])/2);
	v[i]=-1;
      }
    }
  }

  // if(DetectorName()=="TAG-PL") {
  //   std::cout << "segment (" << m_raw->SegmentId() << "), ";
  //   for(int i=0;i<v.size();i++)
  //     std::cout<<v[i] << ", ";
  //   std::cout << std::endl;
  // }

  Int_t index_original_graph = 0;
  for(Int_t i=0;i<size;i++)
    if(v[i]!=-1){
      //      std::cout<<"GX "<<v[i]+0.15<<" GY "<<Original->GXtoGY(v[i]+0.15)<<std::endl;
      Bool_t flagOverThr = false;
      for (Double_t ratio =0; ratio <= 1.0; ratio += 0.1)
	if(GXtoGY(index_original_graph, v[i]+rise*ratio)<thre)
	  flagOverThr = true;

      if (!flagOverThr)
	v[i]=-1;
      /*
      if(GXtoGY(0,v[i]+rise)>thre)
	v[i]=-1;
      */
    }

  // if(DetectorName()=="TAG-PL") {
  //   std::cout << "segment (" << m_raw->SegmentId() << "), ";
  //   for(int i=0;i<v.size();i++)
  //     std::cout<<v[i] << ", ";
  //   std::cout << std::endl;
  // }

  std::sort(v.begin(),v.end(),std::greater<double>());

  std::vector<double>::iterator it = std::find(v.begin(),v.end(),-1);
  if(it!=v.end())
    v.erase(it,v.end());

  std::sort(v.begin(),v.end());
  size=v.size();

}

//_____________________________________________________________________________
Double_t HodoWaveformHit::GXtoGY(Int_t index_graph, Double_t gx)
{
  Int_t point=-1;
  Double_t k,l;
  Int_t sample_size = m_TGraphC[index_graph]->GetN();
  Double_t *GX = m_TGraphC[index_graph]->GetX();
  Double_t *GY = m_TGraphC[index_graph]->GetY();

  for(Int_t i=0;i<sample_size;i++){
    if(gx<=GX[i]){
      point = i ;
      break;
    }
  }

  if(point == -1)
    return 0;

  else{
    k = GX[point-1];
    if(point == sample_size-1)
      return GY[point];
    else
      l = GX[point];

    Double_t r =(gx-k)/(l-k);


    k = GY[point-1];
    l = GY[point];
    return k + (l-k)*r;
  }

  return 0;
}

//_____________________________________________________________________________
Bool_t HodoWaveformHit::SetFitParam(FitParam *fp, std::vector<Double_t> &inix,
				    std::vector<Double_t> &iniy)
{
  Int_t wm = inix.size();
  fp->wavenum = wm;
  fp->ParaNum = wm*2+2;

  if(fp -> ParaNum > ParaMax){
    std::cout<<"Too many Fitting Prameter"<<std::endl;
    return false;
  }

  fp->par[0]=wm;
  fp->par[1]=0;

  for(Int_t i=0;i<wm;i++){
    fp->par[2+2*i]=inix[i];
    fp->par[3+2*i]=iniy[i];
    //    std::cout<<"InitialY"<<InitialY[i]<<"fp y"<<fp->par[3+2*i]<<" iniy "<<iniy[i]<<std::endl;
  }
  return true;

}

//_____________________________________________________________________________
void HodoWaveformHit::Fit1(FitParam *fp)
{
  Int_t wavenum = fp->wavenum;
  Int_t ParaNum = fp->ParaNum;
  Double_t par[ParaNum];
  par[0]=wavenum;

  //std::cout << "wavenum = "  << wavenum << std::endl;

  for(Int_t i=0;i<ParaNum;i++){
    m_func -> ReleaseParameter(i);
    par[i] = fp->par[i];
    //if (seg==14 && m_flag_ch14 && (i>=2 && i%2 == 0))
    //par[i] += 0.08;
    if (i>=2 && i%2 == 1)
      m_func -> SetParLimits(i, 0, 100000);
    //std::cout << "par[" << i << "] = "  << par[i] << std::endl;
  }

  //m_func -> SetName(Name.str().c_str());
  m_func -> SetNpx(1000);
  m_func -> SetParameters(&par[0]);
  m_func -> FixParameter(0,par[0]);
  //m_func -> FixParameter(1,0);
  m_func -> SetLineColor(fp->color);
  for(Int_t nine= ParaNum;nine<ParaMax;nine++)
    m_func -> FixParameter(nine,0);

  m_TGraphC[fp->tgen] -> Fit(m_func,"qN","",fp->FitStart,fp->FitEnd);
  //std::cout<< "fit: " << fp->FitStart << "~" << fp->FitEnd << std::endl;

  for(Int_t i =0;i<ParaNum;i++){
    fp -> FitParam[i] = m_func -> GetParameter(i);
    //std::cout << i << " : " << m_func -> GetParameter(i)
    //<< std::endl;
  }

  //std::cout<< "fit: " << fp->FitStart << "~" << fp->FitEnd  << ", ch" << GetChannel() << ", tgen=" << fp->tgen << std::endl;

}

//_____________________________________________________________________________
Double_t HodoWaveformHit::FittedTrigX(FitParam fp, Double_t allowance)
{
  Int_t num = fp.wavenum;
  Double_t reso = TrigTimeReso *allowance;

  Int_t inrange=0;

  std::vector<Double_t> xx;
  for(Int_t i=0;i<num;i++){
    Double_t x =fp.FitParam[2+2*i];
    //std::cout<<"x "<<x<<"  ref "<<mean<<std::endl;
    if(x > - reso && x< reso ){
      inrange++;
      xx.push_back(x);
    }
  }
  if(inrange==0)
    return -9999.;

  else if(inrange==1)
    return xx[0];
  else{
    std::vector<Double_t> sub;
    for(Int_t i=0;i<inrange;i++)
      sub.push_back(std::abs(xx[i]));
    std::sort(sub.begin(),sub.end());
    for(Int_t i=0;i<inrange;i++){
      if(-sub[0] - xx[i] <0.0001)
        return xx[i];
      if(sub[0] - xx[i] <0.0001)
        return xx[i];
    }
    return -9999.;
  }

}

//_____________________________________________________________________________
Double_t HodoWaveformHit::RisingResidual(Int_t tge_No, Double_t trig, Double_t &res_max)
{

  if(trig <0)
    return -1;

  Double_t Residual = 0, max_res = 0;
  Double_t a,b;

  for(Int_t i =0; i<10;i++){
    Double_t x = trig - 0.1 + i*0.01;
    a = GXtoGY(tge_No, x);
    b = m_func->Eval(x);
    if(a !=0) {
      //Residual += sqrt((a-b)*(a-b))/ CEI->GetError() ;
      Residual += sqrt((a-b)*(a-b))/ 15 ; //kuso tekito

      if (std::abs(max_res) < std::abs(a-b))
        max_res = a-b;

    }
    //std::cout<<"Residual x:"<<x<<"  a:"<<a<<"  b:"<<b<<std::endl;
  }

  Residual /= -GXtoGY(tge_No, trig);
  Residual *= 100;

  res_max = max_res;

  return Residual;
}

//_____________________________________________________________________________
void
HodoWaveformHit::Print(Option_t* arg) const
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
    //              << "mtot          = " << MeanTOT() << std::endl
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
  hddaq::cout << std::endl;
}
