// -*- C++ -*-

#include "EventDisplay.hh"

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include <TApplication.h>
#include <TBRIK.h>
#include <TArc.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TEnv.h>
#include <TFile.h>
#include <TGeometry.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH2Poly.h>
#include <TF1.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TMarker3DBox.h>
#include <TMixture.h>
#include <TNode.h>
#include <TPad.h>
#include <TSystem.h>
#include <TPave.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TPolyLine.h>
#include <TPolyLine3D.h>
#include <TPolyMarker.h>
#include <TPolyMarker3D.h>
#include <TROOT.h>
#include <TRint.h>
#include <TRotMatrix.h>
#include <TStyle.h>
#include <TTRD1.h>
#include <TTRD2.h>
#include <TTUBS.h>
#include <TTUBS.h>
#include <TView.h>
#include <TLine.h>

#include <std_ostream.hh>
#include <UnpackerConfig.hh>
#include <UnpackerXMLReadDigit.hh>

#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "CFTLocalTrack.hh"
#include "CFTParticle.hh"
#include "Exception.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "DeleteUtility.hh"
#include "HodoParamMan.hh"
#include "DCTdcCalibMan.hh"
#include "AftHelper.hh"

#define BH2        0
#define BcOut      1
#define SdcIn      1
#define SdcOut     1
#define TOF        1
#define AC1        0
#define WC         0
#define Vertex     0
#define Hist       0
#define Hist_Timing 0
#define Hist_SdcOut 0
#define Hist_BcIn   0
#define DrawOneHypsTrack 1
#define CATCH        1
#define CATCH3d      1
#define CATCH_Timing 1
#define CATCH_ADC    1
#define BGO_WF       1
#define TAG_WF       1

namespace
{
const auto& gUnpackerConf = hddaq::unpacker::GConfig::get_instance();
const auto& gGeom = DCGeomMan::GetInstance();
  //auto& gAftHelper = AftHelper::GetInstance();
const Int_t& IdTarget  = gGeom.DetectorId("Target");
  //const Int_t& IdBH1     = gGeom.DetectorId("BH1");
  //const Int_t& IdBH2     = gGeom.DetectorId("BH2");
  const Int_t& IdSDC1X3  = gGeom.DetectorId("SDC1-X");
  const Int_t& IdSDC2V   = gGeom.DetectorId("SDC2-V");
  const Int_t& IdSDC2X   = gGeom.DetectorId("SDC2-X");
  const Int_t& IdSDC3X   = gGeom.DetectorId("SDC3-X");
  //const Int_t& IdRKINIT  = gGeom.DetectorId("RKINIT");
const Int_t& IdTOF     = gGeom.DetectorId("TOF-X");
  //const Int_t& IdAC1     = gGeom.DetectorId("AC1");
  //const Int_t& IdWC      = gGeom.DetectorId("WC");
const Double_t& zTarget = gGeom.LocalZ("Target");
  //const Double_t& zK18Target = gGeom.LocalZ("K18Target");
  //const Double_t& gzK18Target = gGeom.GlobalZ("K18Target");
// const Double_t& gxK18Target = gGeom.GetGlobalPosition("K18Target").x();
const Double_t& gxK18Target = -240.;
  //const Double_t& zBFT = gGeom.LocalZ("BFT");

//const Double_t BeamAxis = -150.; //E07
// const Double_t BeamAxis = -240.; //E40
// const Double_t BeamAxis = -50.; //E42
const Double_t BeamAxis = 0.; //E70

#if Vertex
const Double_t MinX = -50.;
const Double_t MaxX =  50.;
const Double_t MinY = -50.;
const Double_t MaxY =  50.;
const Double_t MinZ = -25.;
#endif
const Double_t MaxZ =  50.;

const HodoParamMan& gHodo = HodoParamMan::GetInstance();
const DCTdcCalibMan& gTdc = DCTdcCalibMan::GetInstance();

const double offsetCATCH = 155;
const double offsetBGO   = 60; // offset from CFT

}

//_____________________________________________________________________________
EventDisplay::EventDisplay()
  : m_is_ready(false),
    m_is_save_mode(),
    m_theApp(),
    m_geometry(),
    m_node(),
    m_canvas(),
    m_canvas_vertex(),
    m_canvas_hist(),
    m_canvas_hist2(),
    m_canvas_hist3(),
    m_canvas_hist4(),
    m_canvas_hist5(),
    m_canvas_hist6(),
    m_canvas_hist7(),
    m_canvas_hist8(),
    m_canvas_hist9(),
    m_canvas_hist10(),
    m_hist_vertex_x(),
    m_hist_vertex_y(),
    m_hist_p(),
    m_hist_m2(),
    m_hist_missmass(),
    m_hist_bh1(),
    m_hist_bft(),
    m_hist_bft_p(),
    m_hist_bcIn(),
    m_hist_bcOut(),
    m_BH1box_cont(),
    m_BH2box_cont(),
    m_hist_bh2(),
    m_hist_bcOut_sdcIn(),
    m_hist_sdcIn_predict(),
    m_hist_sdcIn_predict2(),
    m_TargetXZ_box2(),
    m_TargetYZ_box2(),
    m_hist_sch(),
    m_hist_tof(),
    m_hist_sdc1(),
    m_hist_sdc1p(),
    m_hist_sdc3_l(),
    m_hist_sdc3_t(),
    m_hist_sdc3p_l(),
    m_hist_sdc3p_t(),
    m_hist_sdc3y_l(),
    m_hist_sdc3y_t(),
    m_hist_sdc3yp_l(),
    m_hist_sdc3yp_t(),
    m_hist_sdc4_l(),
    m_hist_sdc4_t(),
    m_hist_sdc4p_l(),
    m_hist_sdc4p_t(),
    m_hist_sdc4y_l(),
    m_hist_sdc4y_t(),
    m_hist_sdc4yp_l(),
    m_hist_sdc4yp_t(),
    m_hist_bc3(),
    m_hist_bc3p(),
    m_hist_bc3u(),
    m_hist_bc3up(),
    m_hist_bc3v(),
    m_hist_bc3vp(),
    m_hist_bc4(),
    m_hist_bc4p(),
    m_hist_bc4u(),
    m_hist_bc4up(),
    m_hist_bc4v(),
    m_hist_bc4vp(),
    m_hist_bc3_time(),
    m_hist_bc3p_time(),
    m_hist_bc4_time(),
    m_hist_bc4p_time(),
    m_hist_cft1_l(),
    m_hist_cft1_t(),
    m_hist_cft1_hi(),
    m_hist_cft1_lo(),
    m_hist_cft2_l(),
    m_hist_cft2_t(),
    m_hist_cft2_hi(),
    m_hist_cft2_lo(),
    m_hist_cft3_l(),
    m_hist_cft3_t(),
    m_hist_cft3_hi(),
    m_hist_cft3_lo(),
    m_hist_cft4_l(),
    m_hist_cft4_t(),
    m_hist_cft4_hi(),
    m_hist_cft4_lo(),
    m_hist_cft5_l(),
    m_hist_cft5_t(),
    m_hist_cft5_hi(),
    m_hist_cft5_lo(),
    m_hist_cft6_l(),
    m_hist_cft6_t(),
    m_hist_cft6_hi(),
    m_hist_cft6_lo(),
    m_hist_cft7_l(),
    m_hist_cft7_t(),
    m_hist_cft7_hi(),
    m_hist_cft7_lo(),
    m_hist_cft8_l(),
    m_hist_cft8_t(),
    m_hist_cft8_hi(),
    m_hist_cft8_lo(),
    m_hist_bgo(),
    m_hist_piid_l(),
    m_hist_piid_t(),
    m_target_node(),
    m_BH2wall_node(),
    m_TOFwall_node(),
    m_TOFwall_0_8_node(),
    m_TOFwall_1_9_node(),
    m_TOFwall_10_node(),
    m_TOFwall_11_node(),
    m_TOFwall_12_22_node(),
    m_TOFwall_13_23_node(),
    m_TOFwall_24_34_node(),
    m_TOFwall_25_35_node(),
    m_TOFwall_36_node(),
    m_TOFwall_37_node(),
    m_TOFwall_38_46_node(),
    m_TOFwall_39_47_node(),
    m_AC1_node(),
    m_WCwall_node(),
    m_kurama_inner_node(),
    m_kurama_outer_node(),
    m_BcOutTrack(),
    m_SdcInTrack(),
    m_init_step_mark(),
    m_hs_step_mark(),
    m_Hyps_step_mark(),
    m_Hyps_step_mark_tolast(),
    m_TargetXZ_box(),
    m_TargetYZ_box(),
    m_VertexPointXZ(),
    m_VertexPointYZ(),
    m_HSMarkVertexXShs(),
    m_HypsMarkVertexXShs(),
    m_HypsMarkVertexX(),
    m_HypsMarkVertexY(),
    m_MissMomXZ_line(),
    m_MissMomYZ_line(),
    m_Tgt_Arc(),
    m_CFRP_Arc(),
    m_hbase_catch(),
    m_hbase_catch_zx(),
    m_hbase_catch_zy(),
    m_canvas_dE_E(),
    m_hist_dE_E(),
    m_canvas_catch3d(0),
    m_geometry_catch(0),    
    m_node_catch(0),    
    m_hbase_tagger()
{
}

//_____________________________________________________________________________
EventDisplay::~EventDisplay()
{
}

//_____________________________________________________________________________
// local function
void
ConstructionDone(const TString& name, std::ostream& ost=hddaq::cout)
{
  const Int_t n = 20;
  const Int_t s = name.Length();
  ost << " " << name << " ";
  for(Int_t i=0; i<n-s; ++i) ost << ".";
  ost << " done" << std::endl;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::Initialize()
{
  if(m_is_ready){
    hddaq::cerr << "#W " << FUNC_NAME
                << " already initialied" << std::endl;
    return false;
  }

  gStyle->SetOptStat(0);
  gStyle->SetStatH(0.040);
  gStyle->SetStatX(0.900);
  gStyle->SetStatY(0.900);
  Int_t myfont = 42;
  gStyle->SetTextFont(myfont);
  gStyle->SetLabelFont(myfont, "xyz");
  gStyle->SetLabelSize(0.025, "xyz");
  gStyle->SetTitleFont(myfont, "xyz");
  // gStyle->SetTitleFont(myfont, "p");
  gStyle->SetTitleSize(0.025, "xy");
  gStyle->SetTitleSize(0.016, "z");
  // gStyle->SetTitleSize(0.036, "p");
  gStyle->SetTitleOffset(1.5, "y");
  gStyle->SetTitleOffset(-1.36, "z");
  // gStyle->SetTitleOffset(-0.02, "p");

  gSystem->MakeDirectory("fig/evdisp");

  m_theApp = new TApplication("App", 0, 0);

  // gStyle->SetPalette(kCool);
  gStyle->SetNumberContours(255);

  m_geometry = new TGeometry("evdisp", "HYPS Event Display");

  ThreeVector worldSize(1000., 1000., 1000.); /*mm*/
  new TBRIK("world", "world", "void",
            worldSize.x(), worldSize.y(), worldSize.z());

  m_node = new TNode("node", "node", "world", 0., 0., 0.);
  m_geometry->GetNode("node")->SetVisibility(0);

#if BH2
  ConstructBH2();
#endif

  ConstructTarget();

  //ConstructS2S();
  ConstructKURAMA();

#if BcOut
  //ConstructBcOut();
#endif

#if SdcIn
  ConstructSdcIn();
#endif

#if SdcOut
  ConstructSdcOut();
#endif

#if TOF
  ConstructTOF();
#endif

#if AC1
  //ConstructAC1();
#endif

#if WC
  //ConstructWC();
#endif

  m_canvas = new TCanvas("canvas", "HYPS Event Display",
                         900, 350);
  m_canvas->Divide(2, 1);
  m_canvas->cd(1)->Divide(1, 2);
  m_canvas->cd(1)->cd(1)->SetPad(0.00, 0.92, 1.00, 1.00);
  m_canvas->cd(1)->cd(2)->SetPad(0.00, 0.00, 1.00, 0.92);
  // m_canvas->cd(1)->SetPad(0.00, 0.72, 1.00, 1.00);
  // m_canvas->cd(2)->SetPad(0.00, 0.00, 1.00, 0.72);
  m_canvas->cd(1)->cd(2);

  m_geometry->Draw();

  gPad->SetPhi(180.);
  gPad->SetTheta(-0.);
  gPad->GetView()->ZoomIn();

  m_canvas->cd(2);
  gPad->Update();

  m_canvas->Modified();
  m_canvas->Update();

  m_canvas->cd(2)->Divide(1, 2);
  m_canvas->cd(2)->cd(1)->SetPad(0.00, 0.37, 1.00, 1.00);
  m_canvas->cd(2)->cd(2)->SetPad(0.00, 0.00, 1.00, 0.37);
  ConstructTagger();

#if Vertex

  m_canvas_vertex = new TCanvas("canvas_vertex", "K1.8 Event Display (Vertex)",
                                1000, 800);
  m_canvas_vertex->Divide(1,2);
  m_canvas_vertex->cd(1);
  gPad->DrawFrame(MinZ, MinX, MaxZ, MaxX, "Vertex XZ projection");

  m_TargetXZ_box->Draw("L");

  gPad->Update();

  m_canvas_vertex->cd(2);
  gPad->DrawFrame(MinZ, MinY, MaxZ, MaxY, "Vertex YZ projection");
  m_TargetYZ_box->Draw("L");
  gPad->Update();
#endif

#if Hist
  m_canvas_hist = new TCanvas("canvas_hist", "EventDisplay Hist",
                              400, 800);
  m_canvas_hist->Divide(1,3);
  m_hist_p  = new TH1F("hist1", "Momentum", 100, 0., 3.);
  m_hist_m2 = new TH1F("hist2", "Mass Square", 200, -0.5, 1.5);
  m_hist_missmass = new TH1F("hist3", "Missing Mass", 100, 1, 1.4);
  m_canvas_hist->cd(1)->SetGrid();
  m_hist_p->Draw();
  m_canvas_hist->cd(2)->SetGrid();
  m_hist_m2->Draw();
  m_canvas_hist->cd(3)->SetGrid();
  m_hist_missmass->Draw();
  // gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1111110);
#endif

#if Hist_Timing
  m_canvas_hist2 = new TCanvas("canvas_hist2", "EventDisplay Detector Timing",
                               800, 1000);
  m_canvas_hist2->Divide(3,3);
  m_hist_bh2  = new TH2F("hist_bh2", "BH2", NumOfSegBH2, 0., NumOfSegBH2, 500, -500, 500);
  m_hist_bh2->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_tof = new TH2F("hist_tof", "TOF", NumOfSegTOF, 0, NumOfSegTOF, 500, -500, 500);
  m_hist_tof->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_sdc1 = new TH2F("hist_sdc1", "SDC1", MaxWireSDC1, 0, MaxWireSDC1, 500, -500, 500);
  m_hist_sdc1->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_sdc1p = new TH2F("hist_sdc1p", "SDC1 Xp", MaxWireSDC1, 0, MaxWireSDC1, 500, -500, 500);

  m_hist_sdc1p->SetFillColor(kBlack);
  m_canvas_hist2->cd(1)->SetGrid();
  m_hist_bh2->Draw("box");
  m_canvas_hist2->cd(2)->SetGrid();
  m_hist_sch->Draw("box");
  m_canvas_hist2->cd(8)->SetGrid();
  m_hist_tof->Draw("box");
  m_canvas_hist2->cd(3)->SetGrid();
  m_hist_sdc1->Draw("box");
  m_hist_sdc1p->Draw("samebox");

  m_canvas_hist3 = new TCanvas("canvas_hist3", "EventDisplay Detector Timing",
                               800, 1000);
  m_canvas_hist3->Divide(3,2);

  m_hist_bc3 = new TH2F("hist_bc3", "BC3 X", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500);
  m_hist_bc3p = new TH2F("hist_bc3p", "BC3 Xp", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500);
  m_hist_bc3p->SetFillColor(kBlack);
  m_hist_bc3->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc3p->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_bc3_time = new TH1F("hist_bc3_time", "BC3 X", 500, -500, 500);
  m_hist_bc3p_time = new TH1F("hist_bc3p_time", "BC3 X", 500, -500, 500);
  m_hist_bc3p_time->SetFillStyle(3001);
  m_hist_bc3p_time->SetFillColor(kBlack);

  m_hist_bc3u = new TH2F("hist_bc3u", "BC3 U", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500);
  m_hist_bc3up = new TH2F("hist_bc3up", "BC3 Up", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500);
  m_hist_bc3up->SetFillColor(kBlack);
  m_hist_bc3u->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc3up->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_bc3v = new TH2F("hist_bc3v", "BC3 V", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500);
  m_hist_bc3vp = new TH2F("hist_bc3vp", "BC3 Vp", MaxWireBC3, 0, MaxWireBC3, 500, -500, 500);
  m_hist_bc3vp->SetFillColor(kBlack);
  m_hist_bc3v->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc3vp->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_bc4 = new TH2F("hist_bc4", "BC4 X", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500);
  m_hist_bc4p = new TH2F("hist_bc4p", "BC4 Xp", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500);
  m_hist_bc4p->SetFillColor(kBlack);
  m_hist_bc4->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc4p->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_bc4_time = new TH1F("hist_bc4_time", "BC4 X", 500, -500, 500);
  m_hist_bc4p_time = new TH1F("hist_bc4p_time", "BC4 X", 500, -500, 500);
  m_hist_bc4p_time->SetFillStyle(3001);
  m_hist_bc4p_time->SetFillColor(kBlack);

  m_hist_bc4u = new TH2F("hist_bc4u", "BC4 U", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500);
  m_hist_bc4up = new TH2F("hist_bc4up", "BC4 Up", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500);
  m_hist_bc4up->SetFillColor(kBlack);
  m_hist_bc4u->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc4up->GetYaxis()->SetRangeUser(-100, 100);

  m_hist_bc4v = new TH2F("hist_bc4v", "BC4 V", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500);
  m_hist_bc4vp = new TH2F("hist_bc4vp", "BC4 Vp", MaxWireBC4, 0, MaxWireBC4, 500, -500, 500);
  m_hist_bc4vp->SetFillColor(kBlack);
  m_hist_bc4v->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bc4vp->GetYaxis()->SetRangeUser(-100, 100);

  m_canvas_hist3->cd(1)->SetGrid();
  m_hist_bc3->Draw("box");
  m_hist_bc3p->Draw("samebox");
  m_canvas_hist3->cd(2)->SetGrid();
  m_hist_bc3v->Draw("box");
  m_hist_bc3vp->Draw("samebox");
  m_canvas_hist3->cd(3)->SetGrid();
  m_hist_bc3u->Draw("box");
  m_hist_bc3up->Draw("samebox");
  m_canvas_hist3->cd(4)->SetGrid();
  m_hist_bc4->Draw("box");
  m_hist_bc4p->Draw("samebox");
  m_canvas_hist3->cd(5)->SetGrid();
  m_hist_bc4v->Draw("box");
  m_hist_bc4vp->Draw("samebox");
  m_canvas_hist3->cd(6)->SetGrid();
  m_hist_bc4u->Draw("box");
  m_hist_bc4up->Draw("samebox");
  /*
    m_canvas_hist3->cd(3)->SetGrid();
    m_hist_bc3_time->Draw("");
    m_hist_bc3p_time->Draw("same");
    m_canvas_hist3->cd(4)->SetGrid();
    m_hist_bc4_time->Draw();
    m_hist_bc4p_time->Draw("same");
  */
  // gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

#endif

#if Hist_SdcOut
  m_canvas_hist4 = new TCanvas("canvas_hist4", "EventDisplay Detector Timing (SdcOut)", 800, 800);
  m_canvas_hist4->Divide(2,4);

  m_hist_sdc3_l = new TH2F("hist_sdc3_l", "SDC3 (leading)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3_t = new TH2F("hist_sdc3_t", "SDC3 (trailing)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3_t->SetFillColor(kRed);

  m_hist_sdc3p_l = new TH2F("hist_sdc3p_l", "SDC3 Xp (leading)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3p_t = new TH2F("hist_sdc3p_t", "SDC3 Xp(trailing)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3p_l->SetFillColor(kBlack);
  m_hist_sdc3p_t->SetFillColor(kGreen);

  m_hist_sdc3y_l = new TH2F("hist_sdc3y_l", "SDC3 Y (leading)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3y_t = new TH2F("hist_sdc3y_t", "SDC3 Y (trailing)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3y_t->SetFillColor(kRed);

  m_hist_sdc3yp_l = new TH2F("hist_sdc3yp_l", "SDC3 Yp (leading)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3yp_t = new TH2F("hist_sdc3yp_t", "SDC3 Yp(trailing)", MaxWireSDC3, 0, MaxWireSDC3, 500, -500, 500);
  m_hist_sdc3yp_l->SetFillColor(kBlack);
  m_hist_sdc3yp_t->SetFillColor(kGreen);


  m_hist_sdc4_l = new TH2F("hist_sdc4_l", "SDC4 (leading)", MaxWireSDC4, 0, MaxWireSDC4, 500, -500, 500);
  m_hist_sdc4_t = new TH2F("hist_sdc4_t", "SDC4 (trailing)", MaxWireSDC4, 0, MaxWireSDC4, 500, -500, 500);
  m_hist_sdc4_t->SetFillColor(kRed);

  m_hist_sdc4p_l = new TH2F("hist_sdc4p_l", "SDC4 Xp(leading)", MaxWireSDC4, 0, MaxWireSDC4, 500, -500, 500);
  m_hist_sdc4p_t = new TH2F("hist_sdc4p_t", "SDC4 Xp(trailing)", MaxWireSDC4, 0, MaxWireSDC4, 500, -500, 500);
  m_hist_sdc4p_l->SetFillColor(kBlack);
  m_hist_sdc4p_t->SetFillColor(kGreen);

  m_hist_sdc4y_l = new TH2F("hist_sdc4y_l", "SDC4 Y (leading)", MaxWireSDC4, 0, MaxWireSDC4, 500, -500, 500);
  m_hist_sdc4y_t = new TH2F("hist_sdc4y_t", "SDC4 Y (trailing)", MaxWireSDC4, 0, MaxWireSDC4, 500, -500, 500);
  m_hist_sdc4y_t->SetFillColor(kRed);

  m_hist_sdc4yp_l = new TH2F("hist_sdc4yp_l", "SDC4 Yp(leading)", MaxWireSDC4, 0, MaxWireSDC4, 500, -500, 500);
  m_hist_sdc4yp_t = new TH2F("hist_sdc4yp_t", "SDC4 Yp(trailing)", MaxWireSDC4, 0, MaxWireSDC4, 500, -500, 500);
  m_hist_sdc4yp_l->SetFillColor(kBlack);
  m_hist_sdc4yp_t->SetFillColor(kGreen);

  m_canvas_hist4->cd(5)->SetGrid();
  m_hist_sdc3_l->Draw("box");
  m_hist_sdc3_t->Draw("samebox");
  m_hist_sdc3p_l->Draw("samebox");
  m_hist_sdc3p_t->Draw("samebox");

  m_canvas_hist4->cd(6)->SetGrid();
  m_hist_sdc3y_l->Draw("box");
  m_hist_sdc3y_t->Draw("samebox");
  m_hist_sdc3yp_l->Draw("samebox");
  m_hist_sdc3yp_t->Draw("samebox");

  m_canvas_hist4->cd(7)->SetGrid();
  m_hist_sdc4_l->Draw("box");
  m_hist_sdc4_t->Draw("samebox");
  m_hist_sdc4p_l->Draw("samebox");
  m_hist_sdc4p_t->Draw("samebox");

  m_canvas_hist4->cd(8)->SetGrid();
  m_hist_sdc4y_l->Draw("box");
  m_hist_sdc4y_t->Draw("samebox");
  m_hist_sdc4yp_l->Draw("samebox");
  m_hist_sdc4yp_t->Draw("samebox");


#endif

#if Hist_BcIn
  m_canvas_hist5 = new TCanvas("canvas_hist5", "EventDisplay Detector Timing (BcIn)",
                               800, 1000);
  m_canvas_hist5->Divide(2,2);
  m_hist_bh1  = new TH2F("hist_bh1", "BH1", NumOfSegBH1, 0., NumOfSegBH1, 500, -500, 500);
  m_hist_bh1->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bft  = new TH2F("hist_bft", "BFT", NumOfSegBFT, 0., NumOfSegBFT, 500, -500, 500);
  m_hist_bft_p  = new TH2F("hist_bft_p", "BFT prime", NumOfSegBFT, 0., NumOfSegBFT, 500, -500, 500);
  m_hist_bft->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bft_p->GetYaxis()->SetRangeUser(-100, 100);
  m_hist_bft_p->SetFillColor(kBlack);

  m_hist_bcIn  = new TH2F("hist_bcIn", "BcIn Tracking", 200, -100, 100, 200, -150, 50);

  m_hist_bcOut  = new TH2F("hist_bcOut", "BcOut Tracking", 200, -100, 100, 200, 0, 600);


  m_canvas_hist5->cd(1)->SetGrid();
  m_hist_bh1->Draw("box");
  m_canvas_hist5->cd(2)->SetGrid();
  m_hist_bft->Draw("box");
  m_hist_bft_p->Draw("samebox");
  m_canvas_hist5->cd(3)->SetGrid();
  m_hist_bcIn->Draw("box");

  Double_t Bh1SegX[NumOfSegBH1] = {30./2., 20./2., 16./2., 12./2., 8./2., 8./2., 8./2., 12./2., 16./2., 20./2., 30./2.};
  Double_t Bh1SegY[NumOfSegBH1] = {5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2.};

  //Double_t localPosBh1Z = -515.;
  Double_t localPosBh1Z = gGeom.GetLocalZ(IdBH1);
  Double_t localPosBh1X_dX = 0.;
  Double_t localPosBh1X[NumOfSegBH1] = {-70. + localPosBh1X_dX,
                                        -46. + localPosBh1X_dX,
                                        -29. + localPosBh1X_dX,
                                        -16. + localPosBh1X_dX,
                                        -7. + localPosBh1X_dX,
                                        0. + localPosBh1X_dX,
                                        7. + localPosBh1X_dX,
                                        16. + localPosBh1X_dX,
                                        29. + localPosBh1X_dX,
                                        46. + localPosBh1X_dX,
                                        70. + localPosBh1X_dX};
  Double_t localPosBh1_dZ[NumOfSegBH1] = {4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5};

  for (Int_t i=0; i<NumOfSegBH1; i++) {
    m_BH1box_cont.push_back(new TBox(localPosBh1X[i]-Bh1SegX[i],
                                     localPosBh1Z+localPosBh1_dZ[i]-Bh1SegY[i],
                                     localPosBh1X[i]+Bh1SegX[i],
                                     localPosBh1Z+localPosBh1_dZ[i]+Bh1SegY[i]));
  }
  for (Int_t i=0; i<NumOfSegBH1; i++) {
    m_BH1box_cont[i]->SetFillColor(kWhite);
    m_BH1box_cont[i]->SetLineColor(kBlack);
    m_BH1box_cont[i]->Draw("L");
  }

  m_canvas_hist5->cd(4)->SetGrid();
  m_hist_bcOut->Draw("box");

  Double_t Bh2SegX[NumOfSegBH2] = {35./2., 10./2., 7./2., 7./2., 7./2., 7./2., 10./2., 35./2.};
  Double_t Bh2SegY[NumOfSegBH2] = {5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2.};

  Double_t localPosBh2X[NumOfSegBH2] = {-41.5, -19.0, -10.5, -3.5, 3.5, 10.5, 19.0, 41.5};

  Double_t localPosBh2Z = gGeom.GetLocalZ(IdBH2);

  for (Int_t i=0; i<NumOfSegBH2; i++) {
    m_BH2box_cont.push_back(new TBox(localPosBh2X[i]-Bh2SegX[i],
                                     localPosBh2Z-Bh2SegY[i],
                                     localPosBh2X[i]+Bh2SegX[i],
                                     localPosBh2Z+Bh2SegY[i]));
  }
  for (Int_t i=0; i<NumOfSegBH2; i++) {
    m_BH2box_cont[i]->SetFillColor(kWhite);
    m_BH2box_cont[i]->SetLineColor(kBlack);
    m_BH2box_cont[i]->Draw("L");
  }

  gStyle->SetOptStat(0);


  m_canvas_hist6 = new TCanvas("canvas_hist6", "EventDisplay Detector Timing (BcOut SdcIn)",
                               800, 1000);

  m_hist_sdcIn_predict  = new TH2F("hist_sdcIn_predict", "BcOut-SdcIn Tracking", 200, -400, 400, 750, -3000, 0);
  m_hist_sdcIn_predict->SetFillColor(kRed);
  m_hist_sdcIn_predict->SetLineColor(kRed);
  m_hist_sdcIn_predict->Draw("box");

  m_hist_sdcIn_predict2  = new TH2F("hist_sdcIn_predict2", "BcOut-SdcIn Tracking", 200, -400, 400, 750, -3000, 0);
  m_hist_sdcIn_predict2->SetFillColor(kMagenta);
  m_hist_sdcIn_predict2->SetLineColor(kMagenta);
  m_hist_sdcIn_predict2->Draw("samebox");

  m_hist_bcOut_sdcIn  = new TH2F("hist_bcOut_sdcIn", "BcOut-SdcIn Tracking", 400, -400, 400, 1500, -3000, 0);
  m_canvas_hist6->SetGrid();
  m_hist_bcOut_sdcIn->Draw("samebox");

  Double_t globalPosBh2Z = gGeom.GetGlobalPosition(IdBH2).z();
  Double_t globalPosBh2X = gGeom.GetGlobalPosition(IdBH2).x();

  for (Int_t i=0; i<NumOfSegBH2; i++) {
    m_BH2box_cont2.push_back(new TBox(globalPosBh2X+localPosBh2X[i]-Bh2SegX[i],
                                      globalPosBh2Z-Bh2SegY[i]-20,
                                      globalPosBh2X+localPosBh2X[i]+Bh2SegX[i],
                                      globalPosBh2Z+Bh2SegY[i]+20));
  }
  for (Int_t i=0; i<NumOfSegBH2; i++) {
    m_BH2box_cont2[i]->SetFillColor(kWhite);
    m_BH2box_cont2[i]->SetLineColor(kBlack);
    m_BH2box_cont2[i]->Draw("L");
  }

  Double_t globalPosTarget_x = gGeom.GetGlobalPosition(IdTarget).x();
  Double_t globalPosTarget_y = gGeom.GetGlobalPosition(IdTarget).y();
  Double_t globalPosTarget_z = gGeom.GetGlobalPosition(IdTarget).z();
  Double_t target_r = 20.;
  Double_t target_z = 300./2;
  m_TargetXZ_box2 =  new TBox(globalPosTarget_x-target_r,
                              globalPosTarget_z-target_z,
                              globalPosTarget_x+target_r,
                              globalPosTarget_z+target_z);

  m_TargetXZ_box2->SetFillColor(kWhite);
  m_TargetXZ_box2->SetLineColor(kBlack);
  m_TargetXZ_box2->Draw("L");

  m_TargetYZ_box2 =  new TBox(globalPosTarget_y-target_r,
                              globalPosTarget_z-target_z,
                              globalPosTarget_y+target_r,
                              globalPosTarget_z+target_z);

  m_TargetYZ_box2->SetFillColor(kWhite);
  m_TargetYZ_box2->SetLineColor(kBlack);
  m_TargetYZ_box2->Draw("L");

#endif

#if CATCH_Timing
  m_canvas_hist7 = new TCanvas( "canvas_hist7", "EventDisplay Detector Timing (CATCH)", 400, 500 );
  m_canvas_hist7->Divide(4,3);


  for (Int_t layer=0; layer<8; layer++) {
    TH2 *hp_l=0, *hp_t=0;
    hp_l = new TH2F( Form( "hist_cft%d_l", layer+1 ),
		      Form( "CFT layer%d (Leading)", layer+1 ),
		      NumOfSegCFT[layer], 0, NumOfSegCFT[layer], 500, -500, 500 );
    hp_l->GetYaxis()->SetRangeUser(-100, 100);
    hp_t = new TH2F( Form( "hist_cft%d_t", layer+1 ),
		      Form( "CFT layer%d (Trailing)", layer+1 ),
		      NumOfSegCFT[layer], 0, NumOfSegCFT[layer], 500, -500, 500 );
    hp_t->GetYaxis()->SetRangeUser(-100, 100);
    hp_t->SetFillColor(kRed);

    m_canvas_hist7->cd(layer+1)->SetGrid();
    hp_l->Draw("box");
    hp_t->Draw("samebox");

    if (layer == 0) {
      m_hist_cft1_l = hp_l;
      m_hist_cft1_t = hp_t;
    } else if (layer == 1) {
      m_hist_cft2_l = hp_l;
      m_hist_cft2_t = hp_t;
    } else if (layer == 2) {
      m_hist_cft3_l = hp_l;
      m_hist_cft3_t = hp_t;
    } else if (layer == 3) {
      m_hist_cft4_l = hp_l;
      m_hist_cft4_t = hp_t;
    } else if (layer == 4) {
      m_hist_cft5_l = hp_l;
      m_hist_cft5_t = hp_t;
    } else if (layer == 5) {
      m_hist_cft6_l = hp_l;
      m_hist_cft6_t = hp_t;
    } else if (layer == 6) {
      m_hist_cft7_l = hp_l;
      m_hist_cft7_t = hp_t;
    } else if (layer == 7) {
      m_hist_cft8_l = hp_l;
      m_hist_cft8_t = hp_t;
    }

  }

  m_hist_bgo = new TH2F( "hist_bgo","BGO (Leading)",
		   NumOfSegBGO, 0, NumOfSegBGO, 500, -500, 500 );
  m_canvas_hist7->cd(9)->SetGrid();
  m_hist_bgo->Draw("box");

  m_hist_piid_l = new TH2F( "hist_piid_l","PiId (Leading)",
			   NumOfSegPiID, 0, NumOfSegPiID, 500, -500, 500 );
  m_hist_piid_t = new TH2F( "hist_piid_t","PiId (Leading)",
			   NumOfSegPiID, 0, NumOfSegPiID, 500, -500, 500 );
  m_canvas_hist7->cd(10)->SetGrid();
  m_hist_piid_l->Draw("box");
  m_hist_piid_t->Draw("samebox");
#endif

#if CATCH_ADC

  m_canvas_hist8 = new TCanvas( "canvas_hist8", "EventDisplay Detector ADC (CATCH)", 400, 500 );
  m_canvas_hist8->Divide(4,3);

  for (Int_t layer=0; layer<8; layer++) {
    TH2 *hp_hi=0, *hp_lo=0;
    hp_hi = new TH2F( Form( "hist_cft%d_hi", layer+1 ),
		      Form( "CFT layer%d (High Gain)", layer+1 ),
		      NumOfSegCFT[layer], 0, NumOfSegCFT[layer], 500, -50, 3450 );
    hp_lo = new TH2F( Form( "hist_cft%d_low", layer+1 ),
		      Form( "CFT layer%d (Low Gain)", layer+1 ),
		      NumOfSegCFT[layer], 0, NumOfSegCFT[layer], 500, -50, 3450 );
    hp_lo->SetFillColor(kRed);

    m_canvas_hist8->cd(layer+1)->SetGrid();
    hp_hi->Draw("box");
    hp_lo->Draw("samebox");

    if (layer == 0) {
      m_hist_cft1_hi = hp_hi;
      m_hist_cft1_lo = hp_lo;
    } else if (layer == 1) {
      m_hist_cft2_hi = hp_hi;
      m_hist_cft2_lo = hp_lo;
    } else if (layer == 2) {
      m_hist_cft3_hi = hp_hi;
      m_hist_cft3_lo = hp_lo;
    } else if (layer == 3) {
      m_hist_cft4_hi = hp_hi;
      m_hist_cft4_lo = hp_lo;
    } else if (layer == 4) {
      m_hist_cft5_hi = hp_hi;
      m_hist_cft5_lo = hp_lo;
    } else if (layer == 5) {
      m_hist_cft6_hi = hp_hi;
      m_hist_cft6_lo = hp_lo;
    } else if (layer == 6) {
      m_hist_cft7_hi = hp_hi;
      m_hist_cft7_lo = hp_lo;
    } else if (layer == 7) {
      m_hist_cft8_hi = hp_hi;
      m_hist_cft8_lo = hp_lo;
    }

  }

#endif

#if BGO_WF
  m_canvas_hist9 = new TCanvas( "canvas_hist9", "EventDisplay Detector BGO waveform (CATCH)", 500, 500 );
  m_canvas_hist9->Divide(4,6);
  for (int npad = 0; npad<24; npad++)
    m_canvas_hist9->cd(npad+1)->SetGrid();
#endif

#if TAG_WF
  m_canvas_hist10 = new TCanvas( "canvas_hist10", "EventDisplay Detector Tag-PL waveform", 500, 500 );
  m_canvas_hist10->Divide(2,5);
  for (int npad = 0; npad<10; npad++)
    m_canvas_hist10->cd(npad+1)->SetGrid();
#endif

#if CATCH
  ConstructCATCH();

  m_canvas_dE_E = new TCanvas( "canvas_dE_E", "EventDisplay Detector #DeltaE-E (CATCH)", 250, 250 );
  m_canvas_dE_E->SetGrid();
  m_hist_dE_E = new TH2F( "hist_dE_E","#Delta E - E",
			   100, 0, 200, 100, 0, 10 );
  m_hist_dE_E->Draw();

  const double p1[4] = {2.13544, -0.0357453, 0.000296504, -9.29258e-07};
  const double p2[4] = {4.47497, -0.0789024, 0.000748812, -2.62424e-06};

  TF1 *func1 = new TF1("func1","pol3", 0, 120);
  for (int i=0; i<4; i++)
    func1->SetParameter(i, p1[i]);
  m_CATCH_dE_E_line_cont.push_back(func1);

  TF1 *func2 = new TF1("func2","pol3", 0, 120);
  for (int i=0; i<4; i++)
    func2->SetParameter(i, p2[i]);
  m_CATCH_dE_E_line_cont.push_back(func2);

  TF1 *func3 = new TF1("func3","pol0", 100, 200);
  func3->SetParameter(0, 0.6);
  m_CATCH_dE_E_line_cont.push_back(func3);

  TF1 *func4 = new TF1("func4","pol0", 100, 200);
  func4->SetParameter(0, 1.2);
  m_CATCH_dE_E_line_cont.push_back(func4);

  for (int i=0; i<m_CATCH_dE_E_line_cont.size(); i++) {
    m_CATCH_dE_E_line_cont[i]->SetLineColor(kRed);
    m_CATCH_dE_E_line_cont[i]->Draw("same");
  }

#endif

#if CATCH3d
  ConstructCATCH3d();
#endif
  
  ResetVisibility();

  m_is_ready = true;
  return m_is_ready;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructBH2()
{
  const Int_t lid = gGeom.GetDetectorId("BH2");

  Double_t rotMatBH2[9] = {};
  Double_t BH2wallX = 120.0/2.; // X
  Double_t BH2wallY =   6.0/2.; // Z
  Double_t BH2wallZ =  40.0/2.; // Y
  Double_t BH2SizeX[NumOfSegBH2] = { 120./2. }; // X
  Double_t BH2SizeY[NumOfSegBH2] = {   6./2. }; // Z
  Double_t BH2SizeZ[NumOfSegBH2] = {  40./2. }; // Y
  Double_t BH2PosX[NumOfSegBH2]  = { 0./2. };
  Double_t BH2PosY[NumOfSegBH2]  = { 0./2. };
  Double_t BH2PosZ[NumOfSegBH2]  = { 0./2. };

  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatBH2);

  new TRotMatrix("rotBH2", "rotBH2", rotMatBH2);
  const ThreeVector& BH2wallPos = gGeom.GetGlobalPosition(lid);
  new TBRIK("BH2wall_brik", "BH2wall_brik", "void",
            BH2wallX, BH2wallY, BH2wallZ);
  m_BH2wall_node = new TNode("BH2wall_node", "BH2wall_node", "BH2wall_brik",
                             BH2wallPos.x(),
                             BH2wallPos.y(),
                             BH2wallPos.z(), "rotBH2", "void");
  m_BH2wall_node->SetVisibility(0);
  m_BH2wall_node->cd();

  for(Int_t i=0; i<NumOfSegBH2; ++i){
    new TBRIK(Form("BH2seg_brik_%d", i),
              Form("BH2seg_brik_%d", i),
              "void", BH2SizeX[i], BH2SizeY[i], BH2SizeZ[i]);
    m_BH2seg_node.push_back(new TNode(Form("BH2seg_node_%d", i),
                                      Form("BH2seg_node_%d", i),
                                      Form("BH2seg_brik_%d", i),
                                      BH2PosX[i], BH2PosY[i], BH2PosZ[i]));
  }
  m_node->cd();
  ConstructionDone(__func__);
  return true;
}

//______________________________________________________________________________
Bool_t
EventDisplay::ConstructKURAMA( void )
{
  double Matrix[9] = {};

  double inner_x = 1400.0/2.; // X
  double inner_y =  800.0/2.; // Z
  double inner_z =  800.0/2.; // Y

  double outer_x = 2200.0/2.; // X
  double outer_y =  800.0/2.; // Z
  double outer_z = 1540.0/2.; // Y

  double uguard_inner_x = 1600.0/2.; // X
  double uguard_inner_y =  100.0/2.; // Z
  double uguard_inner_z = 1940.0/2.; // Y
  double uguard_outer_x =  600.0/2.; // X
  double uguard_outer_y =  100.0/2.; // Z
  double uguard_outer_z =  300.0/2.; // Y

  double dguard_inner_x = 1600.0/2.; // X
  double dguard_inner_y =  100.0/2.; // Z
  double dguard_inner_z = 1940.0/2.; // Y
  double dguard_outer_x = 1100.0/2.; // X
  double dguard_outer_y =  100.0/2.; // Z
  double dguard_outer_z = 1100.0/2.; // Y

  CalcRotMatrix( 0., 0., 0., Matrix );

  new TRotMatrix( "rotKURAMA", "rotKURAMA", Matrix );

  new TBRIK( "kurama_inner_brik", "kurama_inner_brik",
	     "void", inner_x, inner_y, inner_z );

  new TBRIK( "kurama_outer_brik", "kurama_outer_brik",
	     "void", outer_x, outer_y, outer_z );

  new TBRIK( "uguard_inner_brik", "uguard_inner_brik",
	     "void", uguard_inner_x, uguard_inner_y, uguard_inner_z );

  new TBRIK( "uguard_outer_brik", "uguard_outer_brik",
	     "void", uguard_outer_x, uguard_outer_y, uguard_outer_z );

  new TBRIK( "dguard_inner_brik", "dguard_inner_brik",
	     "void", dguard_inner_x, dguard_inner_y, dguard_inner_z );

  new TBRIK( "dguard_outer_brik", "dguard_outer_brik",
	     "void", dguard_outer_x, dguard_outer_y, dguard_outer_z );


  m_kurama_inner_node = new TNode( "kurama_inner_node",
				   "kurama_inner_node",
				   "kurama_inner_brik",
				   0., 0., 0., "rotKURAMA", "void" );
  m_kurama_outer_node = new TNode( "kurama_outer_node",
				   "kurama_outer_node",
				   "kurama_outer_brik",
				   0., 0., 0., "rotKURAMA", "void" );
  /*
  TNode *uguard_inner = new TNode( "uguard_inner_node",
				   "uguard_inner_node",
				   "uguard_inner_brik",
				   0., 0., -820.,
				   "rotKURAMA", "void" );
  TNode *uguard_outer = new TNode( "uguard_outer_node",
				   "uguard_outer_node",
				   "uguard_outer_brik",
				   0., 0., -820.,
				   "rotKURAMA", "void" );

  TNode *dguard_inner = new TNode( "dguard_inner_node",
				   "dguard_inner_node",
				   "dguard_inner_brik",
				   0., 0., 820.,
				   "rotKURAMA", "void" );
  TNode *dguard_outer = new TNode( "dguard_outer_node",
				   "dguard_outer_node",
				   "dguard_outer_brik",
				   0., 0., 820.,
				   "rotKURAMA", "void" );
  */
  const Color_t color = kBlack;

  m_kurama_inner_node->SetLineColor( color );
  m_kurama_outer_node->SetLineColor( color );
  //uguard_inner->SetLineColor( color );
  //uguard_outer->SetLineColor( color );
  //dguard_inner->SetLineColor( color );
  //dguard_outer->SetLineColor( color );

  m_node->cd();
  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructS2S()
{
  Double_t Matrix[9] = {};
  CalcRotMatrix(0., 0., 0., Matrix);

  // Q1
  const TVector3 Q1InnerSize( 500./2.,  500./2., 880./2.);
  const TVector3 Q1OuterSize(2400./2., 2400./2., 880./2.);
  const TVector3 posS2SQ1(0., 0., -2775.6 - 1110.);
  new TRotMatrix("rotS2SQ1", "rotS2SQ1", Matrix);
  new TBRIK("q1i_brik", "q1i_brik", "void",
            Q1InnerSize.X(), Q1InnerSize.Z(), Q1InnerSize.Y());
  new TBRIK("q1o_brik", "q1o_brik", "void",
            Q1OuterSize.X(), Q1OuterSize.Z(), Q1OuterSize.Y());
  auto q1i = new TNode("q1i_node", "q1i_node", "q1i_brik",
                       posS2SQ1.X(), posS2SQ1.Y(), posS2SQ1.Z(),
                       "rotS2SQ1", "void");
  auto q1o = new TNode("q1o_node", "q1o_node", "q1o_brik",
                       posS2SQ1.X(), posS2SQ1.Y(), posS2SQ1.Z(),
                       "rotS2SQ1", "void");

  // Q2
  const TVector3 Q2InnerSize( 500./2.,  500./2., 540./2.);
  const TVector3 Q2OuterSize(2100./2., 1540./2., 540./2.);
  const TVector3 posS2SQ2(0., 0., -2775.);
  new TRotMatrix("rotS2SQ2", "rotS2SQ2", Matrix);
  new TBRIK("q2i_brik", "q2i_brik", "void",
            Q2InnerSize.X(), Q2InnerSize.Z(), Q2InnerSize.Y());
  new TBRIK("q2o_brik", "q2o_brik", "void",
            Q2OuterSize.X(), Q2OuterSize.Z(), Q2OuterSize.Y());
  auto q2i = new TNode("q2i_node", "q2i_node", "q2i_brik",
                       posS2SQ2.X(), posS2SQ2.Y(), posS2SQ2.Z(),
                       "rotS2SQ2", "void");
  auto q2o = new TNode("q2o_node", "q2o_node", "q2o_brik",
                       posS2SQ2.X(), posS2SQ2.Y(), posS2SQ2.Z(),
                       "rotS2SQ2", "void");

  // D1
  const Double_t rho       = 3000.;//gGeom.GetLocalZ("S2SD1");
  const Double_t bendAngle = 70.;//gGeom.GetRotAngle2("S2SD1");
  const TVector3 D1InnerSize(800./2., 320./2., 3665.19/2.);
  const TVector3 D1OuterSize(2453.1/2., 1200./2., 3665.19/2.);
  const TVector3 posS2SD1(3000., 0., -2100.623);
  new TRotMatrix("rotS2SD1", "rotS2SD1", Matrix);
  new TTUBS("d1i_tubs", "d1i_tubs",
            "void", rho-D1InnerSize.X(), rho+D1InnerSize.X(), D1InnerSize.Y(),
            180., 180.+bendAngle);
  new TTUBS("d1o_tubs", "d1o_tubs",
            "void", rho-D1OuterSize.X(), rho+D1OuterSize.X(), D1OuterSize.Y(),
            180., 180.+bendAngle);
  auto d1i = new TNode("d1i_node", "d1i_node", "d1i_tubs",
                       posS2SD1.X(), posS2SD1.Y(), posS2SD1.Z(),
                       "rotS2SD1", "void");
  auto d1o = new TNode("d1o_node", "d1o_node", "d1o_tubs",
                       posS2SD1.X(), posS2SD1.Y(), posS2SD1.Z(),
                       "rotS2SD1", "void");

  // D1 End Guard
  const TVector3 D1EGInnerSize(810./2., 330./2., 76./2.);
  const TVector3 D1EGOuterSize(1880./2., 1200./2., 76./2.);
  const auto& posS2SD1EG = gGeom.GetGlobalPosition("S2SD1EG");
  const Double_t RA2 = gGeom.GetRotAngle2("S2SD1EG");
  Double_t MatrixEG[9] = {};
  CalcRotMatrix(0., 0., RA2, MatrixEG);
  new TRotMatrix("rotS2SD1EG", "rotS2SD1EG", MatrixEG);
  new TBRIK("d1ei_brik", "d1ei_brik", "void",
            D1EGInnerSize.X(), D1EGInnerSize.Z(), D1EGInnerSize.Y());
  new TBRIK("d1eo_brik", "d1eo_brik", "void",
            D1EGOuterSize.X(), D1EGOuterSize.Z(), D1EGOuterSize.Y());
  auto d1ei = new TNode("d1ei_node", "d1ei_node", "d1ei_brik",
                        posS2SD1EG.X(), posS2SD1EG.Y(), posS2SD1EG.Z(),
                        "rotS2SD1EG", "void");
  auto d1eo = new TNode("d1eo_node", "d1eo_node", "d1eo_brik",
                        posS2SD1EG.X(), posS2SD1EG.Y(), posS2SD1EG.Z(),
                        "rotS2SD1EG", "void");

  const Color_t q1color = kBlack; // kCyan;
  const Color_t q2color = q1color;
  const Color_t d1color = kBlack; // kBlue;
  q1i->SetLineColor(q1color);
  q1o->SetLineColor(q1color);
  q2i->SetLineColor(q2color);
  q2o->SetLineColor(q2color);
  d1i->SetLineColor(d1color);
  d1o->SetLineColor(d1color);
  d1ei->SetLineColor(d1color);
  d1eo->SetLineColor(d1color);

  m_node->cd();
  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructBcOut()
{
  const Double_t wireL = 200.0;

  const Double_t offsetZ = -3068.4;

  // BC3 X1
  {
    const Int_t lid = gGeom.GetDetectorId("BC3-X1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC3X1", "rotBC3X1", Matrix);
    new TTUBE("BC3X1Tube", "BC3X1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC3x1_node.push_back(new TNode(Form("BC3x1_Node_%d", wire),
                                       Form("BC3x1_Node_%d", wire),
                                       "BC3X1Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC3X1", "void"));
    }
  }

  // BC3 X2
  {
    const Int_t lid = gGeom.GetDetectorId("BC3-X2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC3X2", "rotBC3X2", Matrix);
    new TTUBE("BC3X2Tube", "BC3X2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC3x2_node.push_back(new TNode(Form("BC3x2_Node_%d", wire),
                                       Form("BC3x2_Node_%d", wire),
                                       "BC3X2Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC3X2", "void"));
    }
  }

  // BC3 V1
  {
    const Int_t lid = gGeom.GetDetectorId("BC3-V1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC3V1", "rotBC3V1", Matrix);
    new TTUBE("BC3V1Tube", "BC3V1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC3v1_node.push_back(new TNode(Form("BC3v1_Node_%d", wire),
                                       Form("BC3v1_Node_%d", wire),
                                       "BC3V1Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC3V1", "void"));
    }
  }

  // BC3 V2
  {
    const Int_t lid = gGeom.GetDetectorId("BC3-V2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC3V2", "rotBC3V2", Matrix);
    new TTUBE("BC3V2Tube", "BC3V2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC3v2_node.push_back(new TNode(Form("BC3v2_Node_%d", wire),
                                       Form("BC3v2_Node_%d", wire),
                                       "BC3V2Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC3V2", "void"));
    }
  }

  // BC3 U1
  {
    const Int_t lid = gGeom.GetDetectorId("BC3-U1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC3U1", "rotBC3U1", Matrix);
    new TTUBE("BC3U1Tube", "BC3U1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC3u1_node.push_back(new TNode(Form("BC3u1_Node_%d", wire),
                                       Form("BC3u1_Node_%d", wire),
                                       "BC3U1Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC3U1", "void"));
    }
  }

  // BC3 U2
  {
    const Int_t lid = gGeom.GetDetectorId("BC3-U2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC3U2", "rotBC3U2", Matrix);
    new TTUBE("BC3U2Tube", "BC3U2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC3u2_node.push_back(new TNode(Form("BC3u2_Node_%d", wire),
                                       Form("BC3u2_Node_%d", wire),
                                       "BC3U2Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC3U2", "void"));
    }
  }

  // BC4 U1
  {
    const Int_t lid = gGeom.GetDetectorId("BC4-U1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC4U1", "rotBC4U1", Matrix);
    new TTUBE("BC4U1Tube", "BC4U1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC4; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC4u1_node.push_back(new TNode(Form("BC4u1_Node_%d", wire),
                                       Form("BC4u1_Node_%d", wire),
                                       "BC4U1Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC4U1", "void"));
    }
  }

  // BC4 U2
  {
    const Int_t lid = gGeom.GetDetectorId("BC4-U2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC4U2", "rotBC4U2", Matrix);
    new TTUBE("BC4U2Tube", "BC4U2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC4; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC4u2_node.push_back(new TNode(Form("BC4u2_Node_%d", wire),
                                       Form("BC4u2_Node_%d", wire),
                                       "BC4U2Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC4U2", "void"));
    }
  }

  // BC4 V1
  {
    const Int_t lid = gGeom.GetDetectorId("BC4-V1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC4V1", "rotBC4V1", Matrix);
    new TTUBE("BC4V1Tube", "BC4V1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC4; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC4v1_node.push_back(new TNode(Form("BC4v1_Node_%d", wire),
                                       Form("BC4v1_Node_%d", wire),
                                       "BC4V1Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC4V1", "void"));
    }
  }

  // BC4 V2
  {
    const Int_t lid = gGeom.GetDetectorId("BC4-V2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC4V2", "rotBC4V2", Matrix);
    new TTUBE("BC4V2Tube", "BC4V2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC4; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC4v2_node.push_back(new TNode(Form("BC4v2_Node_%d", wire),
                                       Form("BC4v2_Node_%d", wire),
                                       "BC4V2Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC4V2", "void"));
    }
  }

  // BC4 X1
  {
    const Int_t lid = gGeom.GetDetectorId("BC4-X1");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC4X1", "rotBC4X1", Matrix);
    new TTUBE("BC4X1Tube", "BC4X1Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC4; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC4x1_node.push_back(new TNode(Form("BC4x1_Node_%d", wire),
                                       Form("BC4x1_Node_%d", wire),
                                       "BC4X1Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC4X1", "void"));
    }
  }

  // BC4 X2
  {
    const Int_t lid = gGeom.GetDetectorId("BC4-X2");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotBC4X2", "rotBC4X2", Matrix);
    new TTUBE("BC4X2Tube", "BC4X2Tube", "void", Rmin, Rmax, L);
    for(Int_t wire=1; wire<=MaxWireBC4; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_BC4x2_node.push_back(new TNode(Form("BC4x2_Node_%d", wire),
                                       Form("BC4x2_Node_%d", wire),
                                       "BC4X2Tube",
                                       localPos+BeamAxis,
                                       wireGlobalPos.y(),
                                       wireGlobalPos.z()+offsetZ,
                                       "rotBC4X2", "void"));
    }
  }

  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructSdcIn()
{
  const Double_t wireL = 200.0;

  // SDC0 X1
  {
    const Int_t lid = gGeom.GetDetectorId("SDC0-X");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotSDC0X1", "rotSDC0X1", Matrix);
    new TTUBE("SDC0X1Tube", "SDC0X1Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<=MaxWireSDC0; ++wire){
    for(Int_t wire=0; wire<MaxWireSDC0; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_SDC0x1_node.push_back(new TNode(Form("SDC0x1_Node_%d", wire),
                                        Form("SDC0x1_Node_%d", wire),
                                        "SDC0X1Tube",
                                        wireGlobalPos.x()+localPos,
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotSDC0X1", "void"));
    }
  }

  // SDC0 XP (X2)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC0-XP");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t Z    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotX2", "rotX2", Matrix);
    new TTUBE("SDC0X2Tube", "SDC0X2Tube", "void", Rmin, Rmax, Z);
    //for(Int_t wire=1; wire<=MaxWireSDC0; ++wire){
    for(Int_t wire=0; wire<MaxWireSDC0; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_SDC0x2_node.push_back(new TNode(Form("SDC0x2_Node_%d", wire),
                                        Form("SDC0x2_Node_%d", wire),
                                        "SDC0X2Tube",
                                        wireGlobalPos.x()+localPos,
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotX2", "void"));
    }
  }

  // SDC0 U
  {
    const Int_t lid = gGeom.GetDetectorId("SDC0-U");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t Z    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotU1", "rotU1", Matrix);
    new TTUBE("SDC0U1Tube", "SDC0U1Tube", "void", Rmin, Rmax, Z);
    //for(Int_t wire=1; wire<=MaxWireSDC0; ++wire){
    for(Int_t wire=0; wire<MaxWireSDC0; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_SDC0u1_node.push_back(new TNode(Form("SDC0u1_Node_%d", wire),
                                        Form("SDC0u1_Node_%d", wire),
                                        "SDC0U1Tube",
                                        wireGlobalPos.x()+localPos,
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotU1", "void"));
    }
  }

  // SDC0 UP(U2)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC0-UP");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t Z    = wireL/cos(gGeom.GetTiltAngle(lid)*TMath::DegToRad())/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotU2", "rotU2", Matrix);
    new TTUBE("SDC0U2Tube", "SDC0U2Tube", "void", Rmin, Rmax, Z);
    //for(Int_t wire=1; wire<=MaxWireSDC0; ++wire){
    for(Int_t wire=0; wire<MaxWireSDC0; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireGlobalPos = gGeom.GetGlobalPosition(lid);
      m_SDC0u2_node.push_back(new TNode(Form("SDC0u2_Node_%d", wire),
                                        Form("SDC0u2_Node_%d", wire),
                                        "SDC0U2Tube",
                                        wireGlobalPos.x()+localPos,
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotU2", "void"));
    }
  }

  // SDC1 XPP (X1)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC1-XPP");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotSDC1X1", "rotSDC1X1", Matrix);
    new TTUBE("SDC1X1Tube", "SDC1X1Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<=MaxWireSDC1; ++wire){
    for(Int_t wire=0; wire<MaxWireSDC1; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0, 0);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC1x1_node.push_back(new TNode(Form("SDC1x1_Node_%d", wire),
                                        Form("SDC1x1_Node_%d", wire),
                                        "SDC1X1Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotSDC1X1", "void"));
    }
  }

  // SDC1 V1 (V)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC1-V");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotSDC1V1", "rotSDC1V1", Matrix);
    new TTUBE("SDC1V1Tube", "SDC1V1Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<=MaxWireSDC1; ++wire){
    for(Int_t wire=0; wire<MaxWireSDC1; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0, 0);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC1v1_node.push_back(new TNode(Form("SDC1v1_Node_%d", wire),
                                        Form("SDC1v1_Node_%d", wire),
                                        "SDC1V1Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotSDC1V1", "void"));
    }
  }

  // SDC1 U1 (UP)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC1-UP");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotSDC1U1", "rotSDC1U1", Matrix);
    new TTUBE("SDC1U1Tube", "SDC1U1Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<=MaxWireSDC1; ++wire){
    for(Int_t wire=0; wire<MaxWireSDC1; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0, 0);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC1u1_node.push_back(new TNode(Form("SDC1u1_Node_%d", wire),
                                        Form("SDC1u1_Node_%d", wire),
                                        "SDC1U1Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotSDC1U1", "void"));
    }
  }

  // SDC1 U2 (U)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC1-U");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotSDC1U2", "rotSDC1U2", Matrix);
    new TTUBE("SDC1U2Tube", "SDC1U2Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<=MaxWireSDC1; ++wire){
    for(Int_t wire=0; wire<MaxWireSDC1; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0, 0);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC1u2_node.push_back(new TNode(Form("SDC1u2_Node_%d", wire),
                                        Form("SDC1u2_Node_%d", wire),
                                        "SDC1U2Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotSDC1U2", "void"));
    }
  }

  // SDC1 X2 (XP)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC1-XP");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotSDC1X2", "rotSDC1X2", Matrix);
    new TTUBE("SDC1X2Tube", "SDC1X2Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<=MaxWireSDC1; ++wire){
    for(Int_t wire=0; wire<MaxWireSDC1; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0, 0);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC1x2_node.push_back(new TNode(Form("SDC1x2_Node_%d", wire),
                                        Form("SDC1x2_Node_%d", wire),
                                        "SDC1X2Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotSDC1X2", "void"));
    }
  }

  // SDC1 X3 (X)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC1-X");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotSDC1X3", "rotSDC1X3", Matrix);
    new TTUBE("SDC1X3Tube", "SDC1X3Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<=MaxWireSDC1; ++wire){
    for(Int_t wire=0; wire<MaxWireSDC1; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0, 0);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC1x3_node.push_back(new TNode(Form("SDC1x3_Node_%d", wire),
                                        Form("SDC1x3_Node_%d", wire),
                                        "SDC1X3Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotSDC1X3", "void"));
    }
  }


  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructSdcOut()
{
  const Double_t wireL = 1152.0;

  // SDC2 V1 (V)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC2-V");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotV1", "rotV1", Matrix);
    new TTUBE("SDC2V1Tube", "SDC2V1Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<= MaxWireSDC2; ++wire){
    for(Int_t wire=0; wire< MaxWireSDC2; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC2v1_node.push_back(new TNode(Form("SDC2v1_Node_%d", wire),
                                        Form("SDC2v1_Node_%d", wire),
                                        "SDC2V1Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotV1", "void"));
    }
  }

  // SDC2 U1 (UP)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC2-UP");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotU1", "rotU1", Matrix);
    new TTUBE("SDC2U1Tube", "SDC2U1Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<= MaxWireSDC2; ++wire){
    for(Int_t wire=0; wire< MaxWireSDC2; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC2u1_node.push_back(new TNode(Form("SDC2u1_Node_%d", wire),
                                        Form("SDC2u1_Node_%d", wire),
                                        "SDC2U1Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotU1", "void"));
    }
  }

  // SDC2 U2 (U)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC2-U");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotU2", "rotU2", Matrix);
    new TTUBE("SDC2U2Tube", "SDC2U2Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<= MaxWireSDC2; ++wire){
    for(Int_t wire=0; wire< MaxWireSDC2; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC2u2_node.push_back(new TNode(Form("SDC2u2_Node_%d", wire),
                                        Form("SDC2u2_Node_%d", wire),
                                        "SDC2U2Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotU2", "void"));
    }
  }

  // SDC2 X1 (XP)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC2-XP");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotX1", "rotX1", Matrix);
    new TTUBE("SDC2X1Tube", "SDC2X1Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<= MaxWireSDC2; ++wire){
    for(Int_t wire=0; wire< MaxWireSDC2; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC2x1_node.push_back(new TNode(Form("SDC2x1_Node_%d", wire),
                                        Form("SDC2x1_Node_%d", wire),
                                        "SDC2X1Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotX1", "void"));
    }
  }

  // SDC2 X2 (X)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC2-X");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotX2", "rotX2", Matrix);
    new TTUBE("SDC2X2Tube", "SDC2X2Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<= MaxWireSDC2; ++wire){
    for(Int_t wire=0; wire< MaxWireSDC2; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC2x2_node.push_back(new TNode(Form("SDC2x2_Node_%d", wire),
                                        Form("SDC2x2_Node_%d", wire),
                                        "SDC2X2Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotX2", "void"));
    }
  }


  // SDC3 V1 (V)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC3-V");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotV1", "rotV1", Matrix);
    new TTUBE("SDC3V1Tube", "SDC3V1Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<= MaxWireSDC3; ++wire){
    for(Int_t wire=0; wire< MaxWireSDC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC3v1_node.push_back(new TNode(Form("SDC3v1_Node_%d", wire),
                                        Form("SDC3v1_Node_%d", wire),
                                        "SDC3V1Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotV1", "void"));
    }
  }

  // SDC3 U1 (UP)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC3-UP");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotU1", "rotU1", Matrix);
    new TTUBE("SDC3U1Tube", "SDC3U1Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<= MaxWireSDC3; ++wire){
    for(Int_t wire=0; wire< MaxWireSDC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC3u1_node.push_back(new TNode(Form("SDC3u1_Node_%d", wire),
                                        Form("SDC3u1_Node_%d", wire),
                                        "SDC3U1Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotU1", "void"));
    }
  }

  // SDC3 U2 (U)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC3-U");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotU2", "rotU2", Matrix);
    new TTUBE("SDC3U2Tube", "SDC3U2Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<= MaxWireSDC3; ++wire){
    for(Int_t wire=0; wire< MaxWireSDC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC3u2_node.push_back(new TNode(Form("SDC3u2_Node_%d", wire),
                                        Form("SDC3u2_Node_%d", wire),
                                        "SDC3U2Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotU2", "void"));
    }
  }

  // SDC3 X1 (XP)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC3-XP");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotX1", "rotX1", Matrix);
    new TTUBE("SDC3X1Tube", "SDC3X1Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<= MaxWireSDC3; ++wire){
    for(Int_t wire=0; wire< MaxWireSDC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC3x1_node.push_back(new TNode(Form("SDC3x1_Node_%d", wire),
                                        Form("SDC3x1_Node_%d", wire),
                                        "SDC3X1Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotX1", "void"));
    }
  }

  // SDC3 X2 (X)
  {
    const Int_t lid = gGeom.GetDetectorId("SDC3-X");
    Double_t Rmin = 0.0;
    Double_t Rmax = 0.01;
    Double_t L    = wireL/2.;
    Double_t Matrix[9] = {};
    CalcRotMatrix(gGeom.GetTiltAngle(lid),
                  gGeom.GetRotAngle1(lid),
                  gGeom.GetRotAngle2(lid),
                  Matrix);
    new TRotMatrix("rotX2", "rotX2", Matrix);
    new TTUBE("SDC3X2Tube", "SDC3X2Tube", "void", Rmin, Rmax, L);
    //for(Int_t wire=1; wire<= MaxWireSDC3; ++wire){
    for(Int_t wire=0; wire< MaxWireSDC3; ++wire){
      Double_t localPos = gGeom.CalcWirePosition(lid, wire);
      ThreeVector wireLocalPos(localPos, 0., 0.);
      ThreeVector wireGlobalPos = gGeom.Local2GlobalPos(lid, wireLocalPos);
      m_SDC3x2_node.push_back(new TNode(Form("SDC3x2_Node_%d", wire),
                                        Form("SDC3x2_Node_%d", wire),
                                        "SDC3X2Tube",
                                        wireGlobalPos.x(),
                                        wireGlobalPos.y(),
                                        wireGlobalPos.z(),
                                        "rotX2", "void"));
    }
  }

  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructTarget()
{
  Double_t TargetX = 50.0/2.0;
  Double_t TargetY = 30.0/2.0; // Z
  Double_t TargetZ = 30.0/2.0; // Y
  Double_t rotMatTarget[9];

  new TBRIK("target_brik", "target_brik", "void",
            TargetX, TargetY, TargetZ);

  CalcRotMatrix(0., 0., 0., rotMatTarget);
  new TRotMatrix("rotTarget", "rotTarget", rotMatTarget);

  const Int_t lid = IdTarget;
  ThreeVector GlobalPos = gGeom.GetGlobalPosition(lid);
  m_target_node = new TNode("target_node", "target_node", "target_brik",
                            GlobalPos.x(), GlobalPos.y(), GlobalPos.z(),
                            "rotTarget", "void");

#if Vertex
  m_TargetXZ_box = new TBox(-TargetY, -TargetX, TargetY, TargetX);
  m_TargetXZ_box->SetFillColor(kWhite);
  m_TargetXZ_box->SetLineColor(kBlack);
  m_TargetYZ_box = new TBox(-TargetY, -TargetZ, TargetY, TargetZ);
  m_TargetYZ_box->SetFillColor(kWhite);
  m_TargetYZ_box->SetLineColor(kBlack);
#endif

  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
#if 0
Bool_t
EventDisplay::ConstructTOF()
{
  const Int_t lid = gGeom.GetDetectorId("TOF-X");

  Double_t rotMatTOF[9] = {};
  Double_t TOFwallX =  80.0/2.0*NumOfSegTOF; // X
  Double_t TOFwallY =  30.0/2.0*2.0; // Z
  Double_t TOFwallZ = 1800.0/2.0; // Y

  Double_t TOFSegX =  80.0/2.0; // X
  Double_t TOFSegY =  30.0/2.0; // Z
  Double_t TOFSegZ = 1800.0/2.0; // Y

  Double_t overlap = 5.0;

  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatTOF);

  new TRotMatrix("rotTOF", "rotTOF", rotMatTOF);
  ThreeVector  TOFwallPos = gGeom.GetGlobalPosition(lid);
  // Double_t offset = gGeom.CalcWirePosition(lid, (Double_t)NumOfSegTOF/2.-0.5);
  new TBRIK("TOFwall_brik", "TOFwall_brik", "void",
            TOFwallX, TOFwallY, TOFwallZ);
  m_TOFwall_node = new TNode("TOFwall_node", "TOFwall_node", "TOFwall_brik",
                             TOFwallPos.x(),// + offset,
                             TOFwallPos.y(),
                             TOFwallPos.z(),
                             "rotTOF", "void");

  m_TOFwall_node->SetVisibility(0);
  m_TOFwall_node->cd();

  new TBRIK("TOFseg_brik", "TOFseg_brik", "void",
            TOFSegX, TOFSegY, TOFSegZ);
  for(Int_t i=0; i<NumOfSegTOF; i++){
    ThreeVector tofSegLocalPos((-NumOfSegTOF/2.+i)*(TOFSegX*2.-overlap)+75./2.,
                               (-(i%2)*2+1)*TOFSegY-(i%2)*2+1,
                               0.);
    std::cout << i << " " << tofSegLocalPos << std::endl;
    m_TOFseg_node.push_back(new TNode(Form("TOFseg_node_%d", i),
                                      Form("TOFseg_node_%d", i),
                                      "TOFseg_brik",
                                      tofSegLocalPos.x(),
                                      -tofSegLocalPos.y(),
                                      tofSegLocalPos.z()));
  }

  m_node->cd();
  ConstructionDone(__func__);
  return true;
}
#endif

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructTOF()
{
  Double_t LEPS_TOFSegX =  120.0/2.0; // X
  Double_t LEPS_TOFSegY =  40.0/2.0; // Z
  Double_t LEPS_TOFSegZ = 1800.0/2.0; // Y
  Double_t LEPS_TOFPitch = 110.0; // X

  // segment 0, 2, 4, 6, 8
  Int_t lid = gGeom.GetDetectorId("LEPS-TOF-UX-Tilt-R");
  ThreeVector localPosSeg4(gGeom.CalcWirePosition(lid, 4), 0, 0);
  ThreeVector globalPosSeg4=gGeom.Local2GlobalPos(lid, localPosSeg4);
  std::cout << "globalPosSeg4 = " << globalPosSeg4 << std::endl;
  Double_t rotMatTOF[9] = {};
  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatTOF);
  new TRotMatrix("rotTOF_seg4", "rotTOF_seg4", rotMatTOF);

  Double_t TOFwallX_0_8 =  LEPS_TOFPitch * 8 + LEPS_TOFSegX*2; // X
  Double_t TOFwallY_0_8 =  LEPS_TOFSegY; // Z
  Double_t TOFwallZ_0_8 =  LEPS_TOFSegZ; // Y

  new TBRIK("TOFwall_0_8_brik", "TOFwall_0_8_brik", "void",
            TOFwallX_0_8, TOFwallY_0_8, TOFwallZ_0_8);

  m_TOFwall_0_8_node = new TNode("TOFwall_0_8_node", "TOFwall_0_8_node", "TOFwall_0_8_brik",
                             globalPosSeg4.x(),// + offset,
                             globalPosSeg4.y(),
                             globalPosSeg4.z(),
                             "rotTOF_seg4", "void");

  m_TOFwall_0_8_node->SetVisibility(0);
  m_TOFwall_0_8_node->cd();

  new TBRIK("LEPS_TOFseg_brik", "LEPS_TOFseg_brik", "void",
            LEPS_TOFSegX, LEPS_TOFSegY, LEPS_TOFSegZ);

  for(Int_t i=0; i<=8; i+=2){
    Int_t ref_seg = 4;
    Double_t ref_lpos = gGeom.CalcWirePosition(lid, ref_seg);
    Double_t lpos = gGeom.CalcWirePosition(lid, i);
    ThreeVector tofSegLocalPos(lpos - ref_lpos,
                               0., 0.);
    std::cout << i << " " << tofSegLocalPos << std::endl;
    m_TOFseg_node.push_back(new TNode(Form("TOFseg_node_%d", i),
                                      Form("TOFseg_node_%d", i),
                                      "LEPS_TOFseg_brik",
                                      tofSegLocalPos.x(),
                                      tofSegLocalPos.y(),
                                      tofSegLocalPos.z()));
  }

  m_node->cd();

  // segment 1, 3, 5, 7, 9
  lid = gGeom.GetDetectorId("LEPS-TOF-DX-Tilt-R");
  ThreeVector localPosSeg5(gGeom.CalcWirePosition(lid, 5), 0, 0);
  ThreeVector globalPosSeg5=gGeom.Local2GlobalPos(lid, localPosSeg5);
  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatTOF);
  new TRotMatrix("rotTOF_seg5", "rotTOF_seg5", rotMatTOF);

  Double_t TOFwallX_1_9 =  LEPS_TOFPitch * 8 + LEPS_TOFSegX*2; // X
  Double_t TOFwallY_1_9 =  LEPS_TOFSegY; // Z
  Double_t TOFwallZ_1_9 =  LEPS_TOFSegZ; // Y

  new TBRIK("TOFwall_1_9_brik", "TOFwall_1_9_brik", "void",
            TOFwallX_1_9, TOFwallY_1_9, TOFwallZ_1_9);

  m_TOFwall_1_9_node = new TNode("TOFwall_1_9_node", "TOFwall_1_9_node", "TOFwall_1_9_brik",
                             globalPosSeg5.x(),// + offset,
                             globalPosSeg5.y(),
                             globalPosSeg5.z(),
                             "rotTOF_seg5", "void");

  m_TOFwall_1_9_node->SetVisibility(0);
  m_TOFwall_1_9_node->cd();

  for(Int_t i=1; i<=9; i+=2){
    Int_t ref_seg = 5;
    Double_t ref_lpos = gGeom.CalcWirePosition(lid, ref_seg);
    Double_t lpos = gGeom.CalcWirePosition(lid, i);
    ThreeVector tofSegLocalPos(lpos - ref_lpos,
                               0., 0.);
    std::cout << i << " " << tofSegLocalPos << std::endl;
    m_TOFseg_node.push_back(new TNode(Form("TOFseg_node_%d", i),
                                      Form("TOFseg_node_%d", i),
                                      "LEPS_TOFseg_brik",
                                      tofSegLocalPos.x(),
                                      tofSegLocalPos.y(),
                                      tofSegLocalPos.z()));
  }

  m_node->cd();

  // segment 10
  lid = gGeom.GetDetectorId("LEPS-TOF-UX-R");
  ThreeVector localPosSeg10(gGeom.CalcWirePosition(lid, 10), 0, 0);
  ThreeVector globalPosSeg10=gGeom.Local2GlobalPos(lid, localPosSeg10);
  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatTOF);
  new TRotMatrix("rotTOF_seg10", "rotTOF_seg10", rotMatTOF);

  Double_t TOFwallX_10 =  LEPS_TOFSegX; // X
  Double_t TOFwallY_10 =  LEPS_TOFSegY; // Z
  Double_t TOFwallZ_10 =  LEPS_TOFSegZ; // Y

  new TBRIK("TOFwall_10_brik", "TOFwall_10_brik", "void",
            TOFwallX_10, TOFwallY_10, TOFwallZ_10);

  m_TOFwall_10_node = new TNode("TOFwall_10_node", "TOFwall_10_node", "TOFwall_10_brik",
                             globalPosSeg10.x(),// + offset,
                             globalPosSeg10.y(),
                             globalPosSeg10.z(),
                             "rotTOF_seg10", "void");

  m_TOFwall_10_node->SetVisibility(0);
  m_TOFwall_10_node->cd();

  for(Int_t i=10; i<=10; i++){
    Int_t ref_seg = 10;
    Double_t ref_lpos = gGeom.CalcWirePosition(lid, ref_seg);
    Double_t lpos = gGeom.CalcWirePosition(lid, i);
    ThreeVector tofSegLocalPos(lpos - ref_lpos,
                               0., 0.);
    std::cout << i << " " << tofSegLocalPos << std::endl;
    m_TOFseg_node.push_back(new TNode(Form("TOFseg_node_%d", i),
                                      Form("TOFseg_node_%d", i),
                                      "LEPS_TOFseg_brik",
                                      tofSegLocalPos.x(),
                                      tofSegLocalPos.y(),
                                      tofSegLocalPos.z()));
  }

  m_node->cd();

  // segment 11
  lid = gGeom.GetDetectorId("LEPS-TOF-DX-R");
  ThreeVector localPosSeg11(gGeom.CalcWirePosition(lid, 11), 0, 0);
  ThreeVector globalPosSeg11=gGeom.Local2GlobalPos(lid, localPosSeg11);
  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatTOF);
  new TRotMatrix("rotTOF_seg11", "rotTOF_seg11", rotMatTOF);

  Double_t TOFwallX_11 =  LEPS_TOFSegX; // X
  Double_t TOFwallY_11 =  LEPS_TOFSegY; // Z
  Double_t TOFwallZ_11 =  LEPS_TOFSegZ; // Y

  new TBRIK("TOFwall_11_brik", "TOFwall_11_brik", "void",
            TOFwallX_11, TOFwallY_11, TOFwallZ_11);

  m_TOFwall_11_node = new TNode("TOFwall_11_node", "TOFwall_11_node", "TOFwall_11_brik",
                             globalPosSeg11.x(),// + offset,
                             globalPosSeg11.y(),
                             globalPosSeg11.z(),
                             "rotTOF_seg11", "void");

  m_TOFwall_11_node->SetVisibility(0);
  m_TOFwall_11_node->cd();

  for(Int_t i=11; i<=11; i++){
    Int_t ref_seg = 11;
    Double_t ref_lpos = gGeom.CalcWirePosition(lid, ref_seg);
    Double_t lpos = gGeom.CalcWirePosition(lid, i);
    ThreeVector tofSegLocalPos(lpos - ref_lpos,
                               0., 0.);
    std::cout << i << " " << tofSegLocalPos << std::endl;
    m_TOFseg_node.push_back(new TNode(Form("TOFseg_node_%d", i),
                                      Form("TOFseg_node_%d", i),
                                      "LEPS_TOFseg_brik",
                                      tofSegLocalPos.x(),
                                      tofSegLocalPos.y(),
                                      tofSegLocalPos.z()));
  }

  m_node->cd();

  // segment 12, 14, 16, 18, 20, 22
  lid = gGeom.GetDetectorId("TOF-UX-R");
  ThreeVector localPosSeg12(gGeom.CalcWirePosition(lid, 12), 0, 0);
  ThreeVector globalPosSeg12=gGeom.Local2GlobalPos(lid, localPosSeg12);
  ThreeVector localPosSeg22(gGeom.CalcWirePosition(lid, 22), 0, 0);
  ThreeVector globalPosSeg22=gGeom.Local2GlobalPos(lid, localPosSeg22);
  ThreeVector globalPosSeg12_22 = (globalPosSeg12 + globalPosSeg22)*0.5;

  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatTOF);
  new TRotMatrix("rotTOF_seg12_22", "rotTOF_seg12_22", rotMatTOF);

  Double_t KURAMA_TOFSegX =  80.0/2.0; // X
  Double_t KURAMA_TOFSegY =  30.0/2.0; // Z
  Double_t KURAMA_TOFSegZ = 1800.0/2.0; // Y
  Double_t KURAMA_TOFPitch = 75.0;

  Double_t TOFwallX_12_22 =  KURAMA_TOFPitch * 10 + KURAMA_TOFSegX*2; // X
  Double_t TOFwallY_12_22 =  KURAMA_TOFSegY; // Z
  Double_t TOFwallZ_12_22 =  KURAMA_TOFSegZ; // Y

  new TBRIK("TOFwall_12_22_brik", "TOFwall_12_22_brik", "void",
            TOFwallX_12_22, TOFwallY_12_22, TOFwallZ_12_22);

  m_TOFwall_12_22_node = new TNode("TOFwall_12_22_node", "TOFwall_12_22_node", "TOFwall_12_22_brik",
                             globalPosSeg12_22.x(),// + offset,
                             globalPosSeg12_22.y(),
                             globalPosSeg12_22.z(),
                             "rotTOF_seg12_22", "void");

  m_TOFwall_12_22_node->SetVisibility(0);
  m_TOFwall_12_22_node->cd();

  new TBRIK("KURAMA_TOFseg_brik", "KURAMA_TOFseg_brik", "void",
            KURAMA_TOFSegX, KURAMA_TOFSegY, KURAMA_TOFSegZ);

  for(Int_t i=12; i<=22; i+=2){
    Int_t ref_seg = 17;
    Double_t ref_lpos = gGeom.CalcWirePosition(lid, ref_seg);
    Double_t lpos = gGeom.CalcWirePosition(lid, i);
    ThreeVector tofSegLocalPos(lpos - ref_lpos,
                               0., 0.);
    std::cout << i << " " << tofSegLocalPos << std::endl;
    m_TOFseg_node.push_back(new TNode(Form("TOFseg_node_%d", i),
                                      Form("TOFseg_node_%d", i),
                                      "KURAMA_TOFseg_brik",
                                      tofSegLocalPos.x(),
                                      tofSegLocalPos.y(),
                                      tofSegLocalPos.z()));
  }

  m_node->cd();

  // segment 13, 15, 17, 19, 21, 23
  lid = gGeom.GetDetectorId("TOF-DX-R");
  ThreeVector localPosSeg13(gGeom.CalcWirePosition(lid, 13), 0, 0);
  ThreeVector globalPosSeg13=gGeom.Local2GlobalPos(lid, localPosSeg13);
  ThreeVector localPosSeg23(gGeom.CalcWirePosition(lid, 23), 0, 0);
  ThreeVector globalPosSeg23=gGeom.Local2GlobalPos(lid, localPosSeg23);
  ThreeVector globalPosSeg13_23 = (globalPosSeg13 + globalPosSeg23)*0.5;

  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatTOF);
  new TRotMatrix("rotTOF_seg13_23", "rotTOF_seg13_23", rotMatTOF);

  Double_t TOFwallX_13_23 =  KURAMA_TOFPitch * 10 + KURAMA_TOFSegX*2; // X
  Double_t TOFwallY_13_23 =  KURAMA_TOFSegY; // Z
  Double_t TOFwallZ_13_23 =  KURAMA_TOFSegZ; // Y

  new TBRIK("TOFwall_13_23_brik", "TOFwall_13_23_brik", "void",
            TOFwallX_13_23, TOFwallY_13_23, TOFwallZ_13_23);

  m_TOFwall_13_23_node = new TNode("TOFwall_13_23_node", "TOFwall_13_23_node", "TOFwall_13_23_brik",
                             globalPosSeg13_23.x(),// + offset,
                             globalPosSeg13_23.y(),
                             globalPosSeg13_23.z(),
                             "rotTOF_seg13_23", "void");

  m_TOFwall_13_23_node->SetVisibility(0);
  m_TOFwall_13_23_node->cd();

  for(Int_t i=13; i<=23; i+=2){
    Int_t ref_seg = 18;
    Double_t ref_lpos = gGeom.CalcWirePosition(lid, ref_seg);
    Double_t lpos = gGeom.CalcWirePosition(lid, i);
    ThreeVector tofSegLocalPos(lpos - ref_lpos,
                               0., 0.);
    std::cout << i << " " << tofSegLocalPos << std::endl;
    m_TOFseg_node.push_back(new TNode(Form("TOFseg_node_%d", i),
                                      Form("TOFseg_node_%d", i),
                                      "KURAMA_TOFseg_brik",
                                      tofSegLocalPos.x(),
                                      tofSegLocalPos.y(),
                                      tofSegLocalPos.z()));
  }

  m_node->cd();

  // segment 24, 26, 28, 30, 32, 34
  lid = gGeom.GetDetectorId("TOF-DX-L");
  ThreeVector localPosSeg24(gGeom.CalcWirePosition(lid, 24), 0, 0);
  ThreeVector globalPosSeg24=gGeom.Local2GlobalPos(lid, localPosSeg24);
  ThreeVector localPosSeg34(gGeom.CalcWirePosition(lid, 34), 0, 0);
  ThreeVector globalPosSeg34=gGeom.Local2GlobalPos(lid, localPosSeg34);
  ThreeVector globalPosSeg24_34 = (globalPosSeg24 + globalPosSeg34)*0.5;

  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatTOF);
  new TRotMatrix("rotTOF_seg24_34", "rotTOF_seg24_34", rotMatTOF);

  Double_t TOFwallX_24_34 =  KURAMA_TOFPitch * 10 + KURAMA_TOFSegX*2; // X
  Double_t TOFwallY_24_34 =  KURAMA_TOFSegY; // Z
  Double_t TOFwallZ_24_34 =  KURAMA_TOFSegZ; // Y

  new TBRIK("TOFwall_24_34_brik", "TOFwall_24_34_brik", "void",
            TOFwallX_24_34, TOFwallY_24_34, TOFwallZ_24_34);

  m_TOFwall_24_34_node = new TNode("TOFwall_24_34_node", "TOFwall_24_34_node", "TOFwall_24_34_brik",
                             globalPosSeg24_34.x(),// + offset,
                             globalPosSeg24_34.y(),
                             globalPosSeg24_34.z(),
                             "rotTOF_seg24_34", "void");

  m_TOFwall_24_34_node->SetVisibility(0);
  m_TOFwall_24_34_node->cd();

  for(Int_t i=24; i<=34; i+=2){
    Int_t ref_seg = 29;
    Double_t ref_lpos = gGeom.CalcWirePosition(lid, ref_seg);
    Double_t lpos = gGeom.CalcWirePosition(lid, i);
    ThreeVector tofSegLocalPos(lpos - ref_lpos,
                               0., 0.);
    std::cout << i << " " << tofSegLocalPos << std::endl;
    m_TOFseg_node.push_back(new TNode(Form("TOFseg_node_%d", i),
                                      Form("TOFseg_node_%d", i),
                                      "KURAMA_TOFseg_brik",
                                      tofSegLocalPos.x(),
                                      tofSegLocalPos.y(),
                                      tofSegLocalPos.z()));
  }

  m_node->cd();

  // segment 25, 27, 29, 31, 33, 35
  lid = gGeom.GetDetectorId("TOF-UX-L");
  ThreeVector localPosSeg25(gGeom.CalcWirePosition(lid, 25), 0, 0);
  ThreeVector globalPosSeg25=gGeom.Local2GlobalPos(lid, localPosSeg25);
  ThreeVector localPosSeg35(gGeom.CalcWirePosition(lid, 35), 0, 0);
  ThreeVector globalPosSeg35=gGeom.Local2GlobalPos(lid, localPosSeg35);
  ThreeVector globalPosSeg25_35 = (globalPosSeg25 + globalPosSeg35)*0.5;

  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatTOF);
  new TRotMatrix("rotTOF_seg25_35", "rotTOF_seg25_35", rotMatTOF);

  Double_t TOFwallX_25_35 =  KURAMA_TOFPitch * 10 + KURAMA_TOFSegX*2; // X
  Double_t TOFwallY_25_35 =  KURAMA_TOFSegY; // Z
  Double_t TOFwallZ_25_35 =  KURAMA_TOFSegZ; // Y

  new TBRIK("TOFwall_25_35_brik", "TOFwall_25_35_brik", "void",
            TOFwallX_25_35, TOFwallY_25_35, TOFwallZ_25_35);

  m_TOFwall_25_35_node = new TNode("TOFwall_25_35_node", "TOFwall_25_35_node", "TOFwall_25_35_brik",
                             globalPosSeg25_35.x(),// + offset,
                             globalPosSeg25_35.y(),
                             globalPosSeg25_35.z(),
                             "rotTOF_seg25_35", "void");

  m_TOFwall_25_35_node->SetVisibility(0);
  m_TOFwall_25_35_node->cd();

  for(Int_t i=25; i<=35; i+=2){
    Int_t ref_seg = 30;
    Double_t ref_lpos = gGeom.CalcWirePosition(lid, ref_seg);
    Double_t lpos = gGeom.CalcWirePosition(lid, i);
    ThreeVector tofSegLocalPos(lpos - ref_lpos,
                               0., 0.);
    std::cout << i << " " << tofSegLocalPos << std::endl;
    m_TOFseg_node.push_back(new TNode(Form("TOFseg_node_%d", i),
                                      Form("TOFseg_node_%d", i),
                                      "KURAMA_TOFseg_brik",
                                      tofSegLocalPos.x(),
                                      tofSegLocalPos.y(),
                                      tofSegLocalPos.z()));
  }

  m_node->cd();

  // segment 36
  lid = gGeom.GetDetectorId("LEPS-TOF-DX-L");
  ThreeVector localPosSeg36(gGeom.CalcWirePosition(lid, 36), 0, 0);
  ThreeVector globalPosSeg36=gGeom.Local2GlobalPos(lid, localPosSeg36);
  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatTOF);
  new TRotMatrix("rotTOF_seg36", "rotTOF_seg36", rotMatTOF);

  Double_t TOFwallX_36 =  LEPS_TOFSegX; // X
  Double_t TOFwallY_36 =  LEPS_TOFSegY; // Z
  Double_t TOFwallZ_36 =  LEPS_TOFSegZ; // Y

  new TBRIK("TOFwall_36_brik", "TOFwall_36_brik", "void",
            TOFwallX_36, TOFwallY_36, TOFwallZ_36);

  m_TOFwall_36_node = new TNode("TOFwall_36_node", "TOFwall_36_node", "TOFwall_36_brik",
                             globalPosSeg36.x(),// + offset,
                             globalPosSeg36.y(),
                             globalPosSeg36.z(),
                             "rotTOF_seg36", "void");

  m_TOFwall_36_node->SetVisibility(0);
  m_TOFwall_36_node->cd();

  for(Int_t i=36; i<=36; i++){
    Int_t ref_seg = 36;
    Double_t ref_lpos = gGeom.CalcWirePosition(lid, ref_seg);
    Double_t lpos = gGeom.CalcWirePosition(lid, i);
    ThreeVector tofSegLocalPos(lpos - ref_lpos,
                               0., 0.);
    std::cout << i << " " << tofSegLocalPos << std::endl;
    m_TOFseg_node.push_back(new TNode(Form("TOFseg_node_%d", i),
                                      Form("TOFseg_node_%d", i),
                                      "LEPS_TOFseg_brik",
                                      tofSegLocalPos.x(),
                                      tofSegLocalPos.y(),
                                      tofSegLocalPos.z()));
  }

  m_node->cd();

  // segment 37
  lid = gGeom.GetDetectorId("LEPS-TOF-UX-L");
  ThreeVector localPosSeg37(gGeom.CalcWirePosition(lid, 37), 0, 0);
  ThreeVector globalPosSeg37=gGeom.Local2GlobalPos(lid, localPosSeg37);
  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatTOF);
  new TRotMatrix("rotTOF_seg37", "rotTOF_seg37", rotMatTOF);

  Double_t TOFwallX_37 =  LEPS_TOFSegX; // X
  Double_t TOFwallY_37 =  LEPS_TOFSegY; // Z
  Double_t TOFwallZ_37 =  LEPS_TOFSegZ; // Y

  new TBRIK("TOFwall_37_brik", "TOFwall_37_brik", "void",
            TOFwallX_37, TOFwallY_37, TOFwallZ_37);

  m_TOFwall_37_node = new TNode("TOFwall_37_node", "TOFwall_37_node", "TOFwall_37_brik",
                             globalPosSeg37.x(),// + offset,
                             globalPosSeg37.y(),
                             globalPosSeg37.z(),
                             "rotTOF_seg37", "void");

  m_TOFwall_37_node->SetVisibility(0);
  m_TOFwall_37_node->cd();

  for(Int_t i=37; i<=37; i++){
    Int_t ref_seg = 37;
    Double_t ref_lpos = gGeom.CalcWirePosition(lid, ref_seg);
    Double_t lpos = gGeom.CalcWirePosition(lid, i);
    ThreeVector tofSegLocalPos(lpos - ref_lpos,
                               0., 0.);
    std::cout << i << " " << tofSegLocalPos << std::endl;
    m_TOFseg_node.push_back(new TNode(Form("TOFseg_node_%d", i),
                                      Form("TOFseg_node_%d", i),
                                      "LEPS_TOFseg_brik",
                                      tofSegLocalPos.x(),
                                      tofSegLocalPos.y(),
                                      tofSegLocalPos.z()));
  }

  m_node->cd();

  // segment 38, 40, 42, 44, 46
  lid = gGeom.GetDetectorId("LEPS-TOF-DX-Tilt-L");
  ThreeVector localPosSeg42(gGeom.CalcWirePosition(lid, 42), 0, 0);
  ThreeVector globalPosSeg42=gGeom.Local2GlobalPos(lid, localPosSeg42);
  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatTOF);
  new TRotMatrix("rotTOF_seg42", "rotTOF_seg42", rotMatTOF);

  Double_t TOFwallX_38_46 =  LEPS_TOFPitch * 8 + LEPS_TOFSegX*2; // X
  Double_t TOFwallY_38_46 =  LEPS_TOFSegY; // Z
  Double_t TOFwallZ_38_46 =  LEPS_TOFSegZ; // Y

  new TBRIK("TOFwall_38_46_brik", "TOFwall_38_46_brik", "void",
            TOFwallX_38_46, TOFwallY_38_46, TOFwallZ_38_46);

  m_TOFwall_38_46_node = new TNode("TOFwall_38_46_node", "TOFwall_38_46_node", "TOFwall_38_46_brik",
                             globalPosSeg42.x(),// + offset,
                             globalPosSeg42.y(),
                             globalPosSeg42.z(),
                             "rotTOF_seg42", "void");

  m_TOFwall_38_46_node->SetVisibility(0);
  m_TOFwall_38_46_node->cd();

  for(Int_t i=38; i<=46; i+=2){
    Int_t ref_seg = 42;
    Double_t ref_lpos = gGeom.CalcWirePosition(lid, ref_seg);
    Double_t lpos = gGeom.CalcWirePosition(lid, i);
    ThreeVector tofSegLocalPos(lpos - ref_lpos,
                               0., 0.);
    std::cout << i << " " << tofSegLocalPos << std::endl;
    m_TOFseg_node.push_back(new TNode(Form("TOFseg_node_%d", i),
                                      Form("TOFseg_node_%d", i),
                                      "LEPS_TOFseg_brik",
                                      tofSegLocalPos.x(),
                                      tofSegLocalPos.y(),
                                      tofSegLocalPos.z()));
  }

  m_node->cd();

  // segment 39, 41, 43, 45, 47
  lid = gGeom.GetDetectorId("LEPS-TOF-UX-Tilt-L");
  ThreeVector localPosSeg43(gGeom.CalcWirePosition(lid, 43), 0, 0);
  ThreeVector globalPosSeg43=gGeom.Local2GlobalPos(lid, localPosSeg43);
  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatTOF);
  new TRotMatrix("rotTOF_seg43", "rotTOF_seg43", rotMatTOF);

  Double_t TOFwallX_39_47 =  LEPS_TOFPitch * 8 + LEPS_TOFSegX*2; // X
  Double_t TOFwallY_39_47 =  LEPS_TOFSegY; // Z
  Double_t TOFwallZ_39_47 =  LEPS_TOFSegZ; // Y

  new TBRIK("TOFwall_39_47_brik", "TOFwall_39_47_brik", "void",
            TOFwallX_39_47, TOFwallY_39_47, TOFwallZ_39_47);

  m_TOFwall_39_47_node = new TNode("TOFwall_39_47_node", "TOFwall_39_47_node", "TOFwall_39_47_brik",
                             globalPosSeg43.x(),// + offset,
                             globalPosSeg43.y(),
                             globalPosSeg43.z(),
                             "rotTOF_seg43", "void");

  m_TOFwall_39_47_node->SetVisibility(0);
  m_TOFwall_39_47_node->cd();

  for(Int_t i=39; i<=47; i+=2){
    Int_t ref_seg = 43;
    Double_t ref_lpos = gGeom.CalcWirePosition(lid, ref_seg);
    Double_t lpos = gGeom.CalcWirePosition(lid, i);
    ThreeVector tofSegLocalPos(lpos - ref_lpos,
                               0., 0.);
    std::cout << i << " " << tofSegLocalPos << std::endl;
    m_TOFseg_node.push_back(new TNode(Form("TOFseg_node_%d", i),
                                      Form("TOFseg_node_%d", i),
                                      "LEPS_TOFseg_brik",
                                      tofSegLocalPos.x(),
                                      tofSegLocalPos.y(),
                                      tofSegLocalPos.z()));
  }



  m_node->cd();
  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructAC1()
{
  const Int_t lid = gGeom.GetDetectorId("AC1");

  Double_t rot[9] = {};
  Double_t AC1boxX = 1450.0/2.; // X
  Double_t AC1boxY = 480.0/2.; // Z
  Double_t AC1boxZ = 1250.0/2.; // Y

  Double_t AC1radX = 1450.0/2.; // X
  Double_t AC1radY = 113.0/2.; // Z
  Double_t AC1radZ = 995.0/2.; // Y

  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rot);

  new TRotMatrix("rotAC1", "rotAC1", rot);
  ThreeVector rad_pos = gGeom.GetGlobalPosition(lid);
  ThreeVector offset(0., 0., -113./2.-0.5+480./2); // ?
  offset.RotateY(70.);
  ThreeVector box_pos = rad_pos + offset;
  new TBRIK("AC1box_brik", "AC1box_brik", "void",
            AC1boxX, AC1boxY, AC1boxZ);
  new TBRIK("AC1rad_brik", "AC1rad_brik", "void",
            AC1radX, AC1radY, AC1radZ);
  // new TNode("AC1box_node", "AC1box_node", "AC1box_brik",
  //           box_pos.x(),
  //           box_pos.y(),
  //           box_pos.z(),
  //           "rotAC1", "void");
  m_AC1_node = new TNode("AC1rad_node", "AC1rad_node", "AC1rad_brik",
                         rad_pos.x(),
                         rad_pos.y(),
                         rad_pos.z(),
                        "rotAC1", "void");
  // m_WCwall_node->SetVisibility(0);
  // m_WCwall_node->cd();

  // new TBRIK("WCseg_brik", "WCseg_brik", "void",
  //           WCSegX, WCSegY, WCSegZ);
  // for(Int_t i=0; i<NumOfSegWC; i++){
  //   ThreeVector wcSegLocalPos((-NumOfSegWC/2.+i)*(WCSegX*2.-overlap)+255./4.,
  //                             (-(i%2)*2+1)*WCSegY-(i%2)*2+1,
  //                             0.);
  //   std::cout << i << " " << wcSegLocalPos << std::endl;
  //   m_WCseg_node.push_back(new TNode(Form("WCseg_node_%d", i),
  //                                    Form("WCseg_node_%d", i),
  //                                    "WCseg_brik",
  //                                    wcSegLocalPos.x(),
  //                                    wcSegLocalPos.y(),
  //                                    wcSegLocalPos.z()));
  // }

  m_node->cd();
  ConstructionDone(__func__);
  return true;
}

//_____________________________________________________________________________
Bool_t
EventDisplay::ConstructWC()
{
  const Int_t lid = gGeom.GetDetectorId("WC");

  Double_t rotMatWC[9] = {};
  Double_t WCwallX = 180.0/2.0*NumOfSegWC; // X
  Double_t WCwallY = 230.0/2.0*2.0; // Z
  Double_t WCwallZ = 730.0/2.0; // Y

  Double_t WCSegX = 180.0/2.0; // X
  Double_t WCSegY = 230.0/2.0; // Z
  Double_t WCSegZ = 730.0/2.0; // Y

  Double_t overlap = 180.0/2.0;

  CalcRotMatrix(gGeom.GetTiltAngle(lid),
                gGeom.GetRotAngle1(lid),
                gGeom.GetRotAngle2(lid),
                rotMatWC);

  new TRotMatrix("rotWC", "rotWC", rotMatWC);
  ThreeVector  WCwallPos = gGeom.GetGlobalPosition(lid);
  // Double_t offset = gGeom.CalcWirePosition(lid, (Double_t)NumOfSegWC/2.-0.5);
  new TBRIK("WCwall_brik", "WCwall_brik", "void",
            WCwallX, WCwallY, WCwallZ);
  m_WCwall_node = new TNode("WCwall_node", "WCwall_node", "WCwall_brik",
                             WCwallPos.x(),// + offset,
                             WCwallPos.y(),
                             WCwallPos.z(),
                             "rotWC", "void");

  m_WCwall_node->SetVisibility(0);
  m_WCwall_node->cd();

  new TBRIK("WCseg_brik", "WCseg_brik", "void",
            WCSegX, WCSegY, WCSegZ);
  for(Int_t i=0; i<NumOfSegWC; i++){
    ThreeVector wcSegLocalPos((-NumOfSegWC/2.+i)*(WCSegX*2.-overlap)+255./4.,
                              (-(i%2)*2+1)*WCSegY-(i%2)*2+1,
                              0.);
    std::cout << i << " " << wcSegLocalPos << std::endl;
    m_WCseg_node.push_back(new TNode(Form("WCseg_node_%d", i),
                                     Form("WCseg_node_%d", i),
                                     "WCseg_brik",
                                     wcSegLocalPos.x(),
                                     wcSegLocalPos.y(),
                                     wcSegLocalPos.z()));
  }

  m_node->cd();
  ConstructionDone(__func__);
  return true;
}

//______________________________________________________________________________
Bool_t EventDisplay::ConstructCATCH(void)
{
  m_canvas->cd(2)->cd(1)->Divide(1,3);
  m_canvas->cd(2)->cd(1)->cd(1)->SetPad( 0.001,0.001,0.499,0.999);
  m_canvas->cd(2)->cd(1)->cd(2)->SetPad( 0.501,0.501,0.999,0.999);
  m_canvas->cd(2)->cd(1)->cd(3)->SetPad( 0.501,0.001,0.999,0.499);
  m_canvas->cd(2)->cd(1)->cd(1)->SetGrid();

  m_hbase_catch = new TH2F("hbase_catch","Event Display XY plane", 180, -180, 180, 180, -180, 180);
  m_hbase_catch->SetMaximum(200);
  m_hbase_catch->SetMinimum(-1);
  m_hbase_catch->Draw();

  m_Tgt_Arc = new TArc(0, 0, 40./2);
  m_Tgt_Arc->SetLineColor(kCyan);
  m_Tgt_Arc->SetFillStyle(0);
  m_Tgt_Arc->Draw("same");

  m_CFRP_Arc = new TArc(0, 0, 80./2);
  m_CFRP_Arc->SetLineColor(kBlack);
  m_CFRP_Arc->SetFillStyle(0);
  m_CFRP_Arc->SetLineWidth(2);
  m_CFRP_Arc->Draw("same");

  ConstructCFT();
  ConstructBGO();
  ConstructPiID();

  //m_canvas_catch->cd(2)->SetGrid();
  m_canvas->cd(2)->cd(1)->cd(2)->SetGrid();
  m_hbase_catch_zx = new TH2F("hbase_catch_zx","Event Display ZX plane", 600, -200, 400, 180, -180, 180);
  m_hbase_catch_zx->Draw();

  //m_canvas_catch->cd(3)->SetGrid();
  m_canvas->cd(2)->cd(1)->cd(3)->SetGrid();
  m_hbase_catch_zy = new TH2F("hbase_catch_zy","Event Display ZY plane", 600, -200, 400, 180, -180, 180);
  m_hbase_catch_zy->Draw();


  return true;
}

//______________________________________________________________________________
Bool_t EventDisplay::ConstructCFT(void)
{
  //m_canvas_catch->cd(1);
  m_canvas->cd(2)->cd(1)->cd(1);

  for (Int_t i=0; i<NumOfPlaneCFT; ++i) {
    m_CFT_Arc_cont[i].reserve(NumOfSegCFT[i]);
    for (Int_t seg=0; seg<NumOfSegCFT[i]; ++seg) {
      Double_t x, y;
      FiberPosPhi(i, seg, &x, &y);
      TArc *arc = new TArc(x, y, 0.75/2);
      arc->SetLineColor(kBlack);
      if (seg%32 == 0)
	arc->SetLineColor(kOrange);
      arc->SetFillStyle(0);
      arc->Draw("same");
      m_CFT_Arc_cont[i].push_back(arc);
    }
  }

  return true;
}

//______________________________________________________________________________
Bool_t EventDisplay::ConstructBGO(void)
{
  //m_canvas_catch->cd(1);
  m_canvas->cd(2)->cd(1)->cd(1);

  int unit=0;

  for (int i=0; i<NumOfBGOUnit; i++) {
    double theta = (double)i*45.;

    for (int j=0; j<NumOfBGOInOneUnit; j++) {
      double x0 = RadiusOfBGOSurface+BGO_Y/2;
      double y0 = (double)(j-0.5)*BGO_X;

      double x1 = x0+BGO_Y/2;
      double y1 = y0+BGO_X/2;

      double x2 = x0-BGO_Y/2;
      double y2 = y0+BGO_X/2;

      double x3 = x0-BGO_Y/2;
      double y3 = y0-BGO_X/2;

      double x4 = x0+BGO_Y/2;
      double y4 = y0-BGO_X/2;

      ThreeVector pos1((x1*cos(theta*TMath::DegToRad()) - y1*sin(theta*TMath::DegToRad())),
		       (x1*sin(theta*TMath::DegToRad()) + y1*cos(theta*TMath::DegToRad())),
		       0);
      ThreeVector pos2((x2*cos(theta*TMath::DegToRad()) - y2*sin(theta*TMath::DegToRad())),
		       (x2*sin(theta*TMath::DegToRad()) + y2*cos(theta*TMath::DegToRad())),
		       0);
      ThreeVector pos3((x3*cos(theta*TMath::DegToRad()) - y3*sin(theta*TMath::DegToRad())),
		       (x3*sin(theta*TMath::DegToRad()) + y3*cos(theta*TMath::DegToRad())),
		       0);
      ThreeVector pos4((x4*cos(theta*TMath::DegToRad()) - y4*sin(theta*TMath::DegToRad())),
		       (x4*sin(theta*TMath::DegToRad()) + y4*cos(theta*TMath::DegToRad())),
		       0);

      TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
      TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
      TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
      TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
      l1->Draw("same");
      l2->Draw("same");
      l3->Draw("same");
      l4->Draw("same");

      unit = j+3*i;
      m_BGO_Line_cont[unit].push_back(l1);
      m_BGO_Line_cont[unit].push_back(l2);
      m_BGO_Line_cont[unit].push_back(l3);
      m_BGO_Line_cont[unit].push_back(l4);

      //std::cout << unit << std::endl;
    }
  }

  for (int i=0; i<NumOfBGOUnit; i++) {
    double theta = 22.5 + (double)i*45.;

    for (int j=0; j<NumOfBGOInOneUnit2; j++) {
      double x0 = RadiusOfBGOSurface2+BGO_Y/2;
      double y0 = (double)(j)*BGO_X;

      double x1 = x0+BGO_Y/2;
      double y1 = y0+BGO_X/2;

      double x2 = x0-BGO_Y/2;
      double y2 = y0+BGO_X/2;

      double x3 = x0-BGO_Y/2;
      double y3 = y0-BGO_X/2;

      double x4 = x0+BGO_Y/2;
      double y4 = y0-BGO_X/2;

      ThreeVector pos1((x1*cos(theta*TMath::DegToRad()) - y1*sin(theta*TMath::DegToRad())),
		       (x1*sin(theta*TMath::DegToRad()) + y1*cos(theta*TMath::DegToRad())),
		       0);
      ThreeVector pos2((x2*cos(theta*TMath::DegToRad()) - y2*sin(theta*TMath::DegToRad())),
		       (x2*sin(theta*TMath::DegToRad()) + y2*cos(theta*TMath::DegToRad())),
		       0);
      ThreeVector pos3((x3*cos(theta*TMath::DegToRad()) - y3*sin(theta*TMath::DegToRad())),
		       (x3*sin(theta*TMath::DegToRad()) + y3*cos(theta*TMath::DegToRad())),
		       0);
      ThreeVector pos4((x4*cos(theta*TMath::DegToRad()) - y4*sin(theta*TMath::DegToRad())),
		       (x4*sin(theta*TMath::DegToRad()) + y4*cos(theta*TMath::DegToRad())),
		       0);

      TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
      TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
      TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
      TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
      l1->Draw("same");
      l2->Draw("same");
      l3->Draw("same");
      l4->Draw("same");

      unit = j+NumOfBGOInOneUnit+3*i;
      m_BGO_Line_cont[unit].push_back(l1);
      m_BGO_Line_cont[unit].push_back(l2);
      m_BGO_Line_cont[unit].push_back(l3);
      m_BGO_Line_cont[unit].push_back(l4);

      //std::cout << unit << std::endl;
    }
  }

  return true;
}

//______________________________________________________________________________
Bool_t EventDisplay::ConstructPiID(void)
{
  //m_canvas_catch->cd(1);
  m_canvas->cd(2)->cd(1)->cd(1);

  int unit=0;

  for (int i=0; i<NumOfPiIDUnit; i++) {
    double theta = (double)i*45.;

    for (int j=0; j<NumOfPiIDInOneUnit; j++) {
      double x0 = RadiusOfPiIDSurface+PiID_Y/2;
      double y0 = (double)(j-1)*PiID_X;

      double x1 = x0+PiID_Y/2;
      double y1 = y0+PiID_X/2;

      double x2 = x0-PiID_Y/2;
      double y2 = y0+PiID_X/2;

      double x3 = x0-PiID_Y/2;
      double y3 = y0-PiID_X/2;

      double x4 = x0+PiID_Y/2;
      double y4 = y0-PiID_X/2;

      ThreeVector pos1((x1*cos(theta*TMath::DegToRad()) - y1*sin(theta*TMath::DegToRad())),
		       (x1*sin(theta*TMath::DegToRad()) + y1*cos(theta*TMath::DegToRad())),
		       0);
      ThreeVector pos2((x2*cos(theta*TMath::DegToRad()) - y2*sin(theta*TMath::DegToRad())),
		       (x2*sin(theta*TMath::DegToRad()) + y2*cos(theta*TMath::DegToRad())),
		       0);
      ThreeVector pos3((x3*cos(theta*TMath::DegToRad()) - y3*sin(theta*TMath::DegToRad())),
		       (x3*sin(theta*TMath::DegToRad()) + y3*cos(theta*TMath::DegToRad())),
		       0);
      ThreeVector pos4((x4*cos(theta*TMath::DegToRad()) - y4*sin(theta*TMath::DegToRad())),
		       (x4*sin(theta*TMath::DegToRad()) + y4*cos(theta*TMath::DegToRad())),
		       0);

      TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
      TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
      TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
      TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
      l1->Draw("same");
      l2->Draw("same");
      l3->Draw("same");
      l4->Draw("same");

      unit = i*NumOfPiIDInOneUnit + j + i;
      m_PiID_Line_cont[unit].push_back(l1);
      m_PiID_Line_cont[unit].push_back(l2);
      m_PiID_Line_cont[unit].push_back(l3);
      m_PiID_Line_cont[unit].push_back(l4);

      //std::cout << unit << std::endl;
    }
  }

  for (int i=0; i<NumOfPiIDUnit; i++) {
    double theta = 22.5 + (double)i*45.;

    for (int j=0; j<NumOfPiIDInOneUnit2; j++) {
      double x0 = RadiusOfPiID2Surface+PiID2_Y/2;
      double y0 = (double)(j)*PiID2_X;

      double x1 = x0+PiID2_Y/2;
      double y1 = y0+PiID2_X/2;

      double x2 = x0-PiID2_Y/2;
      double y2 = y0+PiID2_X/2;

      double x3 = x0-PiID2_Y/2;
      double y3 = y0-PiID2_X/2;

      double x4 = x0+PiID2_Y/2;
      double y4 = y0-PiID2_X/2;

      ThreeVector pos1((x1*cos(theta*TMath::DegToRad()) - y1*sin(theta*TMath::DegToRad())),
		       (x1*sin(theta*TMath::DegToRad()) + y1*cos(theta*TMath::DegToRad())),
		       0);
      ThreeVector pos2((x2*cos(theta*TMath::DegToRad()) - y2*sin(theta*TMath::DegToRad())),
		       (x2*sin(theta*TMath::DegToRad()) + y2*cos(theta*TMath::DegToRad())),
		       0);
      ThreeVector pos3((x3*cos(theta*TMath::DegToRad()) - y3*sin(theta*TMath::DegToRad())),
		       (x3*sin(theta*TMath::DegToRad()) + y3*cos(theta*TMath::DegToRad())),
		       0);
      ThreeVector pos4((x4*cos(theta*TMath::DegToRad()) - y4*sin(theta*TMath::DegToRad())),
		       (x4*sin(theta*TMath::DegToRad()) + y4*cos(theta*TMath::DegToRad())),
		       0);

      TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
      TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
      TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
      TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
      l1->Draw("same");
      l2->Draw("same");
      l3->Draw("same");
      l4->Draw("same");

      unit = (i+1)*NumOfPiIDInOneUnit + i*NumOfPiIDInOneUnit2 + j;
      m_PiID_Line_cont[unit].push_back(l1);
      m_PiID_Line_cont[unit].push_back(l2);
      m_PiID_Line_cont[unit].push_back(l3);
      m_PiID_Line_cont[unit].push_back(l4);

      //std::cout << unit << std::endl;
    }
  }

  return true;
}


//______________________________________________________________________________

Bool_t EventDisplay::ConstructCATCH3d(void)
{
  static const std::string func_name("[EventDisplay::ConstructCATCH3d()]");

  m_canvas_catch3d = new TCanvas( "canvas_catch3d", "CATCH Event Display",
				400, 400 );


  m_geometry_catch = new TGeometry( "evdisp_catch","CATCH Event Display" );

  ThreeVector worldSize( 200., 200., 400. ); /*mm*/
  new TBRIK( "world_catch", "world_catch", "void",
	     worldSize.x(), worldSize.y(), worldSize.z() );

  m_node_catch = new TNode( "node_catch", "node_catch", "world_catch", 0., 0., 0. );
  m_geometry_catch->GetNode("node_catch")->SetVisibility(0);


  double Rmin = 0.0;
  double Rmax = 0.75/2;
  double L    = 400./2.;

  new TTUBE( "CFTFiberTube", "CFTFiberTube", "void", Rmin, Rmax, L );

  for (int i=0; i<NumOfPlaneCFT/2; i++) {
    int layer = 2*i+1;
    m_CFT_node_cont[i].reserve(NumOfSegCFT[layer]);
    for (int seg=0; seg<NumOfSegCFT[layer]; seg++) {
      double x, y;
      FiberPosPhi(layer, seg, &x, &y);

      m_CFT_node_cont[i].push_back( new TNode( Form( "CFT%d_Node_%d", i, seg ),
					       Form( "CFT%d_Node_%d", i, seg ),
					       "CFTFiberTube",
					       x, y, L - offsetCATCH));
    }
  }

  new TBRIK( "BGOBrik", "BGOBrik", "void", BGO_Y/2., BGO_Z/2., BGO_X/2. );
  for (int i=0; i<NumOfBGOUnit; i++) {
    double theta = (double)i*45.;

    double rotMat[9] = {};
    CalcRotMatrix( theta, 0., 0., rotMat );

    new TRotMatrix( Form( "rotBGO_%d", i ),
		    Form( "rotBGO_%d", i ),
		    rotMat);

    for (int j=0; j<NumOfBGOInOneUnit; j++) {
      double x0 = RadiusOfBGOSurface+BGO_Y/2;
      double y0 = (double)(j-0.5)*BGO_X;

      ThreeVector pos0((x0*cos(theta*TMath::DegToRad()) - y0*sin(theta*TMath::DegToRad())),
		       (x0*sin(theta*TMath::DegToRad()) + y0*cos(theta*TMath::DegToRad())),
		       L+offsetBGO - offsetCATCH);

      m_BGOseg_node_cont.push_back( new TNode( Form( "BGOseg_node_%d", i*NumOfBGOInOneUnit + j + i ),
					       Form( "BGOseg_node_%d", i*NumOfBGOInOneUnit + j + i ),
					       "BGOBrik",
					       pos0.x(),
					       pos0.y(),
					       pos0.z(),
					       Form( "rotBGO_%d", i ),
					       "void") );
      //std::cout << unit << std::endl;
    }
  }

  for (int i=0; i<NumOfBGOUnit; i++) {
    double theta = 22.5 + (double)i*45.;

    double rotMat[9] = {};
    CalcRotMatrix( theta, 0., 0., rotMat );
    new TRotMatrix( Form( "rotBGO2_%d", i ),
		    Form( "rotBGO2_%d", i ),
		    rotMat);

    for (int j=0; j<NumOfBGOInOneUnit2; j++) {
      double x0 = RadiusOfBGOSurface2+BGO_Y/2;
      double y0 = (double)(j)*BGO_X;

      ThreeVector pos0((x0*cos(theta*TMath::DegToRad()) - y0*sin(theta*TMath::DegToRad())),
		       (x0*sin(theta*TMath::DegToRad()) + y0*cos(theta*TMath::DegToRad())),
		       L+offsetBGO - offsetCATCH);

      m_BGOseg_node_cont.push_back( new TNode( Form( "BGOseg_node_%d", (i+1)*NumOfBGOInOneUnit + i*NumOfBGOInOneUnit2 + j ),
					       Form( "BGOseg_node_%d", (i+1)*NumOfBGOInOneUnit + i*NumOfBGOInOneUnit2 + j ),
					       "BGOBrik",
					       pos0.x(),
					       pos0.y(),
					       pos0.z(),
					       Form( "rotBGO2_%d", i ),
					       "void") );
      //std::cout << unit << std::endl;
    }
  }

  std::string node_name;
  for (int seg=0; seg<2; seg++) {
    node_name = Form( "BGOseg_node_%d", seg );
    
    TNode *node = m_geometry_catch->GetNode( node_name.c_str() );
    if( !node ){
      hddaq::cout << "#E " << func_name << " "
		  << "no such node : " << node_name << std::endl;
      return false;
    }
    node->SetVisibility(0);    
  }

  new TBRIK( "PiIDBrik", "PiIDBrik", "void", PiID_Y/2., PiID_Z/2., PiID_X/2. );
  for (int i=0; i<NumOfPiIDUnit; i++) {
    double theta = (double)i*45.;

    double rotMat[9] = {};
    CalcRotMatrix( theta, 0., 0., rotMat );

    new TRotMatrix( Form( "rotPiID_%d", i ),
		    Form( "rotPiID_%d", i ),
		    rotMat);

    for (int j=0; j<NumOfPiIDInOneUnit; j++) {
      double x0 = RadiusOfPiIDSurface+PiID_Y/2;
      double y0 = (double)(j-1)*PiID_X;

      ThreeVector pos0((x0*cos(theta*TMath::DegToRad()) - y0*sin(theta*TMath::DegToRad())),
		       (x0*sin(theta*TMath::DegToRad()) + y0*cos(theta*TMath::DegToRad())),
		       L+offsetBGO - offsetCATCH);

      m_PiIDseg_node_cont.push_back( new TNode( Form( "PiIDseg_node_%d", i*NumOfPiIDInOneUnit + j + i ),
					       Form( "PiIDseg_node_%d", i*NumOfPiIDInOneUnit + j + i ),
					       "PiIDBrik",
					       pos0.x(),
					       pos0.y(),
					       pos0.z(),
					       Form( "rotPiID_%d", i ),
					       "void") );



      //std::cout << unit << std::endl;
    }
  }

  new TBRIK( "PiIDBrik2", "PiIDBrik2", "void", PiID2_Y/2., PiID2_Z/2., PiID2_X/2. );
  for (int i=0; i<NumOfPiIDUnit; i++) {
    double theta = 22.5 + (double)i*45.;

    double rotMat[9] = {};
    CalcRotMatrix( theta, 0., 0., rotMat );
    new TRotMatrix( Form( "rotPiID2_%d", i ),
		    Form( "rotPiID2_%d", i ),
		    rotMat);

    for (int j=0; j<NumOfPiIDInOneUnit2; j++) {
      double x0 = RadiusOfPiID2Surface+PiID_Y/2;
      double y0 = (double)(j)*PiID_X;

      ThreeVector pos0((x0*cos(theta*TMath::DegToRad()) - y0*sin(theta*TMath::DegToRad())),
		       (x0*sin(theta*TMath::DegToRad()) + y0*cos(theta*TMath::DegToRad())),
		       L+offsetBGO - offsetCATCH);

      m_PiIDseg_node_cont.push_back( new TNode( Form( "PiIDseg_node_%d", (i+1)*NumOfPiIDInOneUnit + i*NumOfPiIDInOneUnit2 + j ),
					       Form( "PiIDseg_node_%d", (i+1)*NumOfPiIDInOneUnit + i*NumOfPiIDInOneUnit2 + j ),
					       "PiIDBrik2",
					       pos0.x(),
					       pos0.y(),
					       pos0.z(),
					       Form( "rotPiID2_%d", i ),
					       "void") );
      //std::cout << unit << std::endl;
    }
  }

  for (int seg=0; seg<3; seg++) {
    node_name = Form( "PiIDseg_node_%d", seg );
    
    TNode *node = m_geometry_catch->GetNode( node_name.c_str() );
    if( !node ){
      hddaq::cout << "#E " << func_name << " "
		  << "no such node : " << node_name << std::endl;
      return false;
    }
    node->SetVisibility(0);    
  }

  m_geometry_catch->Draw();
  m_canvas_catch3d->Update();

  return true;
}

//_____________________________________________________________________________
void EventDisplay::FiberPosPhi(int layer, int seg, double *x, double *y) const
{
  int lnum;

  if (layer == CFT_PHI1) {
    lnum = gGeom.GetDetectorId("CFT-PHI1");
  } else if (layer == CFT_PHI2) {
    lnum = gGeom.GetDetectorId("CFT-PHI2");
  } else if (layer == CFT_PHI3) {
    lnum = gGeom.GetDetectorId("CFT-PHI3");
  } else if (layer == CFT_PHI4) {
    lnum = gGeom.GetDetectorId("CFT-PHI4");
  } else if (layer == CFT_U1) {
    lnum = gGeom.GetDetectorId("CFT-UV1");
  } else if (layer == CFT_V2) {
    lnum = gGeom.GetDetectorId("CFT-UV2");
  } else if (layer == CFT_U3) {
    lnum = gGeom.GetDetectorId("CFT-UV3");
  } else if (layer == CFT_V4) {
    lnum = gGeom.GetDetectorId("CFT-UV4");
  } else {
    std::cout <<  "EventDisplay::FiberPosPhi : No PHI Layer " <<  layer << std::endl;
    return;
  }

  if (layer == CFT_PHI1 || layer == CFT_PHI2 || layer == CFT_PHI3 || layer == CFT_PHI4) {
    double phi = gGeom.CalcWirePosition(lnum, seg);
    double r   = gGeom.GetLocalZ(lnum);
    if (seg%2 == 0)
      r -= 0.4;
    else
      r += 0.4;

    *x = r * cos(phi*TMath::DegToRad());
    *y = r * sin(phi*TMath::DegToRad());
  } else if (layer == CFT_U1 || layer == CFT_U3) {
    double z = 200.;

    double SegNumU=NumOfSegCFT[layer];
    double phi = -(360./SegNumU)*(double)seg + 90.;
    double offset = 360./400.*z;

    phi += offset;

    double r   = gGeom.GetLocalZ(lnum);
    if (seg%2 == 0)
      r -= 0.4755/2;
    else
      r += 0.4755/2;

    *x = r * cos(phi*TMath::DegToRad());
    *y = r * sin(phi*TMath::DegToRad());

  } else if (layer == CFT_V2 || layer == CFT_V4) {
    double z = 200.;

    double SegNumV=NumOfSegCFT[layer];
    double phi = (360./SegNumV)*(double)seg + 90.;
    double offset = -360./400.*z;

    phi += offset;

    double r   = gGeom.GetLocalZ(lnum);
    if (seg%2 == 0)
      r -= 0.4755/2;
    else
      r += 0.4755/2;

    *x = r * cos(phi*TMath::DegToRad());
    *y = r * sin(phi*TMath::DegToRad());

  }
  /*
  int lnum=301+layer;

  double R     = gGeom.GetLocalZ(lnum);
  double Phi = gGeom.CalcWirePosition(lnum, seg);

  *x = R * cos(Phi*TMath::DegToRad());
  *y = R * sin(Phi*TMath::DegToRad());
  */
  /*
  std::cout << "FiberPos layer" << layer << ", seg" << seg
	    << "(" << R * cos(Theta*Deg2Rad)
	    << "," << R * sin(Theta*Deg2Rad)
	    << ")" << std::endl;
  */

}

//______________________________________________________________________________

void EventDisplay::BGOPos(int seg, double *x, double *y) const
{

  int UnitNum = seg/(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);
  int SegInUnit = seg%(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);

  double theta;
  double x0;
  double y0;

  if (SegInUnit==0 || SegInUnit==1 ) {
    theta = (double)UnitNum*45.;
    x0 = RadiusOfBGOSurface+BGO_Y/2;
    y0 = (double)((double)SegInUnit-0.5)*BGO_X;
  } else {
    theta = 22.5+(double)UnitNum*45.;
    x0 = RadiusOfBGOSurface2+BGO_Y/2;
    y0 = 0.;

  }


  *x = x0*cos(theta*TMath::DegToRad()) - y0*sin(theta*TMath::DegToRad());
  *y = x0*sin(theta*TMath::DegToRad()) + y0*cos(theta*TMath::DegToRad());

}

//______________________________________________________________________________
Bool_t EventDisplay::ConstructTagger(void)
{
  m_canvas->cd(2)->cd(2);

  m_hbase_tagger = new TH2F("m_hbase_tagger", "Tagger", 100, -40, 60, 100, -30, 70);
  m_hbase_tagger->SetMinimum(-1);
  m_hbase_tagger->Draw("");

  double TagSF_X = 6.;
  double TagSF_Y = 1.;

  for (int i=0; i<NumOfSegTagSF; i++) {
    double x0 = 0.0;
    double y0 = i*TagSF_Y;

    double x1 = x0+TagSF_X/2;
    double y1 = y0+TagSF_Y/2;

    double x2 = x0-TagSF_X/2;
    double y2 = y0+TagSF_Y/2;

    double x3 = x0-TagSF_X/2;
    double y3 = y0-TagSF_Y/2;

    double x4 = x0+TagSF_X/2;
    double y4 = y0-TagSF_Y/2;


    double theta = 0.;

    ThreeVector pos1((x1*cos(theta*TMath::DegToRad()) - y1*sin(theta*TMath::DegToRad())),
		     (x1*sin(theta*TMath::DegToRad()) + y1*cos(theta*TMath::DegToRad())),
		     0);
    ThreeVector pos2((x2*cos(theta*TMath::DegToRad()) - y2*sin(theta*TMath::DegToRad())),
		     (x2*sin(theta*TMath::DegToRad()) + y2*cos(theta*TMath::DegToRad())),
		     0);
    ThreeVector pos3((x3*cos(theta*TMath::DegToRad()) - y3*sin(theta*TMath::DegToRad())),
		     (x3*sin(theta*TMath::DegToRad()) + y3*cos(theta*TMath::DegToRad())),
		     0);
    ThreeVector pos4((x4*cos(theta*TMath::DegToRad()) - y4*sin(theta*TMath::DegToRad())),
		     (x4*sin(theta*TMath::DegToRad()) + y4*cos(theta*TMath::DegToRad())),
		     0);

    TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
    TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
    TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
    TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
    l1->Draw("same");
    l2->Draw("same");
    l3->Draw("same");
    l4->Draw("same");

    m_TagSFF_Line_cont[i].push_back(l1);
    m_TagSFF_Line_cont[i].push_back(l2);
    m_TagSFF_Line_cont[i].push_back(l3);
    m_TagSFF_Line_cont[i].push_back(l4);
  }

  for (int i=0; i<NumOfSegTagSF; i++) {
    double x0 = 6.0;
    double y0 = 0.5 + i*TagSF_Y;

    double x1 = x0+TagSF_X/2;
    double y1 = y0+TagSF_Y/2;

    double x2 = x0-TagSF_X/2;
    double y2 = y0+TagSF_Y/2;

    double x3 = x0-TagSF_X/2;
    double y3 = y0-TagSF_Y/2;

    double x4 = x0+TagSF_X/2;
    double y4 = y0-TagSF_Y/2;


    double theta = 0.;

    ThreeVector pos1((x1*cos(theta*TMath::DegToRad()) - y1*sin(theta*TMath::DegToRad())),
		     (x1*sin(theta*TMath::DegToRad()) + y1*cos(theta*TMath::DegToRad())),
		     0);
    ThreeVector pos2((x2*cos(theta*TMath::DegToRad()) - y2*sin(theta*TMath::DegToRad())),
		     (x2*sin(theta*TMath::DegToRad()) + y2*cos(theta*TMath::DegToRad())),
		     0);
    ThreeVector pos3((x3*cos(theta*TMath::DegToRad()) - y3*sin(theta*TMath::DegToRad())),
		     (x3*sin(theta*TMath::DegToRad()) + y3*cos(theta*TMath::DegToRad())),
		     0);
    ThreeVector pos4((x4*cos(theta*TMath::DegToRad()) - y4*sin(theta*TMath::DegToRad())),
		     (x4*sin(theta*TMath::DegToRad()) + y4*cos(theta*TMath::DegToRad())),
		     0);

    TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
    TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
    TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
    TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
    l1->Draw("same");
    l2->Draw("same");
    l3->Draw("same");
    l4->Draw("same");

    m_TagSFB_Line_cont[i].push_back(l1);
    m_TagSFB_Line_cont[i].push_back(l2);
    m_TagSFB_Line_cont[i].push_back(l3);
    m_TagSFB_Line_cont[i].push_back(l4);
  }

  double TagPL_X = 3.;
  double TagPL_Y = 7.4;
  double overlap = 2.7;
  for (int i=0; i<NumOfSegTagPL; i++) {
    double x0 = 24.0 + 10*(i%2);
    double y0 = i*(TagPL_Y - overlap);

    double x1 = x0+TagPL_X/2;
    double y1 = y0+TagPL_Y/2;

    double x2 = x0-TagPL_X/2;
    double y2 = y0+TagPL_Y/2;

    double x3 = x0-TagPL_X/2;
    double y3 = y0-TagPL_Y/2;

    double x4 = x0+TagPL_X/2;
    double y4 = y0-TagPL_Y/2;


    double theta = 0.;

    ThreeVector pos1((x1*cos(theta*TMath::DegToRad()) - y1*sin(theta*TMath::DegToRad())),
		     (x1*sin(theta*TMath::DegToRad()) + y1*cos(theta*TMath::DegToRad())),
		     0);
    ThreeVector pos2((x2*cos(theta*TMath::DegToRad()) - y2*sin(theta*TMath::DegToRad())),
		     (x2*sin(theta*TMath::DegToRad()) + y2*cos(theta*TMath::DegToRad())),
		     0);
    ThreeVector pos3((x3*cos(theta*TMath::DegToRad()) - y3*sin(theta*TMath::DegToRad())),
		     (x3*sin(theta*TMath::DegToRad()) + y3*cos(theta*TMath::DegToRad())),
		     0);
    ThreeVector pos4((x4*cos(theta*TMath::DegToRad()) - y4*sin(theta*TMath::DegToRad())),
		     (x4*sin(theta*TMath::DegToRad()) + y4*cos(theta*TMath::DegToRad())),
		     0);

    TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
    TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
    TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
    TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
    l1->Draw("same");
    l2->Draw("same");
    l3->Draw("same");
    l4->Draw("same");

    m_TagPL_Line_cont[i].push_back(l1);
    m_TagPL_Line_cont[i].push_back(l2);
    m_TagPL_Line_cont[i].push_back(l3);
    m_TagPL_Line_cont[i].push_back(l4);
  }

  return true;
}


//_____________________________________________________________________________
void
EventDisplay::DrawInitTrack()
{
  m_canvas->cd(1)->cd(2);
  if(m_init_step_mark) m_init_step_mark->Draw();

  m_canvas->Update();

}

//_____________________________________________________________________________
void
EventDisplay::DrawHitWire(Int_t lid, Int_t hit_wire,
                          Bool_t range_check, Bool_t tdc_check)
{
  if(hit_wire<=0) return;

  TString node_name;

  const TString bcout_node_name[NumOfLayersBcOut]
    = { Form("BC3x1_Node_%d", hit_wire),
        Form("BC3x2_Node_%d", hit_wire),
        Form("BC3u1_Node_%d", hit_wire),
        Form("BC3u2_Node_%d", hit_wire),
        Form("BC3v1_Node_%d", hit_wire),
        Form("BC3v2_Node_%d", hit_wire),
        Form("BC4x1_Node_%d", hit_wire),
        Form("BC4x2_Node_%d", hit_wire),
        Form("BC4u1_Node_%d", hit_wire),
        Form("BC4u2_Node_%d", hit_wire),
        Form("BC4v1_Node_%d", hit_wire),
        Form("BC4v2_Node_%d", hit_wire) };

  const TString sdcin_node_name[NumOfLayersSdcIn]
    = { Form("SDC0x1_Node_%d", hit_wire),
        Form("SDC0x2_Node_%d", hit_wire),
        Form("SDC0u1_Node_%d", hit_wire),
        Form("SDC0u2_Node_%d", hit_wire),
        Form("SDC1x1_Node_%d", hit_wire),
        Form("SDC1v1_Node_%d", hit_wire),
        Form("SDC1u1_Node_%d", hit_wire),
        Form("SDC1u2_Node_%d", hit_wire),
        Form("SDC1x2_Node_%d", hit_wire),
        Form("SDC1x3_Node_%d", hit_wire)};

  const TString sdcout_node_name[]
    = { Form("SDC2v1_Node_%d", hit_wire),
        Form("SDC2u1_Node_%d", hit_wire),
        Form("SDC2u2_Node_%d", hit_wire),
        Form("SDC2x1_Node_%d", hit_wire),
        Form("SDC2x2_Node_%d", hit_wire),
        Form("SDC3v1_Node_%d", hit_wire),
        Form("SDC3u1_Node_%d", hit_wire),
        Form("SDC3u2_Node_%d", hit_wire),
        Form("SDC3x1_Node_%d", hit_wire),
        Form("SDC3x2_Node_%d", hit_wire)
  };

  switch (lid) {

    // SDC0
  case 1: case 2: case 3: case 4:
    if(hit_wire>MaxWireSDC0) return;
    node_name = sdcin_node_name[lid-1];
    break;

    // SDC1
  case 5: case 6: case 7: case 8: case 9: case 10:
    if(hit_wire>MaxWireSDC1) return;
    node_name = sdcin_node_name[lid-1];
    break;

    // SDC2
  case 31: case 32: case 33: case 34: case 35:
    if(hit_wire>MaxWireSDC2) return;
    node_name = sdcout_node_name[lid-31];
    break;

    // SDC3
  case 36: case 37: case 38: case 39: case 40:
    if(hit_wire>MaxWireSDC3) return;
    node_name = sdcout_node_name[lid-31];
    break;

    // BC3
  case 113: case 114: case 115: case 116: case 117: case 118:
    if(hit_wire>MaxWireBC4) return;
    node_name = bcout_node_name[lid-113];
    break;

    // BC4
  case 119: case 120: case 121: case 122: case 123: case 124:
    if(hit_wire>MaxWireBC3) return;
    node_name = bcout_node_name[lid-113];
    break;

  default:
    throw Exception(FUNC_NAME + Form(" no such plane : %d", lid));
  }

  auto node = m_geometry->GetNode(node_name);
  if(!node){
    throw Exception(FUNC_NAME + Form(" no such node : %s", node_name.Data()));
  }

  node->SetVisibility(1);
  if(range_check && tdc_check)
    node->SetLineColor(kBlue);
  else if(range_check && !tdc_check)
    node->SetLineColor(28);
  else
    node->SetLineColor(kBlack);

  // m_canvas->cd(2);
  // m_geometry->Draw();
  m_canvas->Update();
}

//_____________________________________________________________________________
void
EventDisplay::DrawHitHodoscope(Int_t lid, Int_t seg, Int_t Tu, Int_t Td)
{
  TString node_name;
  if(seg<0) return;
  /*
  if(lid == IdBH2){
    node_name = Form("BH2seg_node_%d", seg);
  }else if(lid == IdTOF){
    node_name = Form("TOFseg_node_%d", seg);
  }else if(lid == IdWC){
    node_name = Form("WCseg_node_%d", seg);
  }else{
    throw Exception(FUNC_NAME + Form(" no such plane : %d", lid));
  }
  */
  if(lid == IdTOF){
    node_name = Form("TOFseg_node_%d", seg);
  }else{
    throw Exception(FUNC_NAME + Form(" no such plane : %d", lid));
  }

  auto node = m_geometry->GetNode(node_name);
  if(!node){
    throw Exception(FUNC_NAME + Form(" no such node : %s",
                                     node_name.Data()));
  }

  node->SetVisibility(1);

  if(Tu>0 && Td>0){
    node->SetLineWidth(2);
    node->SetLineColor(kBlue);
  }
  else {
    node->SetLineWidth(2);
    node->SetLineColor(kGreen);
  }

  m_canvas->cd(1)->cd(2);
  m_geometry->Draw();
  m_canvas->Update();
}

//_____________________________________________________________________________
#if 0
void
EventDisplay::DrawBcOutLocalTrack(const DCLocalTrack *tp)
{
  const Double_t offsetZ = -6526.5; // 5326.5(S2S-FF) + 1200.0(FF-V0)

  Double_t z0 = gGeom.GetLocalZ("BC3-X1") - 100.;
  Double_t x0 = tp->GetX(z0) + BeamAxis;
  Double_t y0 = tp->GetY(z0);
  z0 += offsetZ;

  Double_t z1 = zK18Target + 100.;
  Double_t x1 = tp->GetX(z1) + BeamAxis;
  Double_t y1 = tp->GetY(z1);
  z1 += offsetZ;

  ThreeVector gPos0(x0, y0, z0);
  ThreeVector gPos1(x1, y1, z1);

  auto p = new TPolyLine3D(2);
  p->SetLineColor(kRed);
  p->SetLineWidth(1);
  p->SetPoint(0, gPos0.x(), gPos0.y(), gPos0.z());
  p->SetPoint(1, gPos1.x(), gPos1.y(), gPos1.z());
  m_BcOutTrack.push_back(p);
  m_canvas->cd(1)->cd(2);
  p->Draw();

#if Vertex
  z0 = zK18Target + MinZ;
  z1 = zK18Target;
  x0 = tp->GetX(z0); y0 = tp->GetY(z0);
  x1 = tp->GetX(z1); y1 = tp->GetY(z1);
  z0 -= zK18Target;
  z1 -= zK18Target;
  {
    m_canvas_vertex->cd(1);
    TPolyLine *line = new TPolyLine(2);
    line->SetPoint(0, z0, x0);
    line->SetPoint(1, z1, x1);
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    line->Draw();
    m_BcOutXZ_line.push_back(line);
  }
  {
    m_canvas_vertex->cd(2);
    TPolyLine *line = new TPolyLine(2);
    line->SetPoint(0, z0, y0);
    line->SetPoint(1, z1, y1);
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    line->Draw();
    m_BcOutYZ_line.push_back(line);
  }
#endif
}
#endif
//_____________________________________________________________________________
// void
// EventDisplay::DrawBcOutLocalTrack(Double_t x0, Double_t y0, Double_t u0, Double_t v0)
// {
//   const Double_t offsetZ = -3108.4;

//   Double_t z0 = offsetZ;

//   Double_t z1 = zK18Target + 100.;
//   Double_t x1 = x0+u0*z1;
//   Double_t y1 = y0+v0*z1;
//   z1 += offsetZ;

//   ThreeVector gPos0(x0+BeamAxis, y0, z0);
//   ThreeVector gPos1(x1+BeamAxis, y1, z1);

//   TPolyLine3D *p = new TPolyLine3D(2);
//   p->SetLineColor(kRed);
//   p->SetLineWidth(1);
//   p->SetPoint(0, gPos0.x(), gPos0.y(), gPos0.z());
//   p->SetPoint(1, gPos1.x(), gPos1.y(), gPos1.z());
//   m_BcOutTrack.push_back(p);
//   m_canvas->cd(1)->cd(2);
//   p->Draw();
//   gPad->Update();

// #if Vertex
//   Double_t z2 = zK18Target + MinZ;
//   z1 = zK18Target;
//   Double_t x2 = x0+u0*z0, y2 = y0+v0*z0;
//   x1 = x0+u0*z1; y1 = y0+v0*z1;
//   z2 -= zK18Target;
//   z1 -= zK18Target;
//   {
//     TPolyLine *line = new TPolyLine(2);
//     line->SetPoint(0, z2, x2);
//     line->SetPoint(1, z1, x1);
//     line->SetLineColor(kRed);
//     line->SetLineWidth(1);
//     m_BcOutXZ_line.push_back(line);
//     m_canvas_vertex->cd(1);
//     line->Draw();
//   }
//   {
//     m_canvas_vertex->cd(2);
//     TPolyLine *line = new TPolyLine(2);
//     line->SetPoint(0, z2, y2);
//     line->SetPoint(1, z1, y1);
//     line->SetLineColor(kRed);
//     line->SetLineWidth(1);
//     line->Draw();
//     m_BcOutYZ_line.push_back(line);
//   }
//   m_canvas_vertex->Update();
// #endif

// }

//_____________________________________________________________________________
void
EventDisplay::DrawSdcInLocalTrack(const DCLocalTrack *tp)
{
#if SdcIn
  Double_t z0 = gGeom.GetLocalZ(IdTarget);
  Double_t x0 = tp->GetX(z0), y0 = tp->GetY(z0);
  ThreeVector gPos0 = gGeom.Local2GlobalPos(IdTarget, TVector3(x0, y0, 0.));

  Double_t z1 = gGeom.GetLocalZ(IdSDC1X3);
  Double_t x1 = tp->GetX(z1), y1 = tp->GetY(z1);
  ThreeVector gPos1 = gGeom.Local2GlobalPos(IdSDC1X3,  TVector3(x1, y1, 0.));

  auto p = new TPolyLine3D(2);
  p->SetLineColor(kRed);
  p->SetLineWidth(1);
  p->SetPoint(0, gPos0.x(), gPos0.y(), gPos0.z());
  p->SetPoint(1, gPos1.x(), gPos1.y(), gPos1.z());
  m_SdcInTrack.push_back(p);
  m_canvas->cd(1)->cd(2);
  p->Draw();
  gPad->Update();
#endif

#if Vertex
  Double_t z0 = zTarget;
  Double_t z1 = zTarget + MaxZ;
  x0 = tp->GetX(z0); y0 = tp->GetY(z0);
  x1 = tp->GetX(z1); y1 = tp->GetY(z1);
  z0 -= zTarget;
  z1 -= zTarget;
  {
    m_canvas_vertex->cd(1);
    TPolyLine *line = new TPolyLine(2);
    line->SetPoint(0, z0, x0);
    line->SetPoint(1, z1, x1);
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    line->Draw();
    m_SdcInXZ_line.push_back(line);
  }
  {
    m_canvas_vertex->cd(2);
    TPolyLine *line = new TPolyLine(2);
    line->SetPoint(0, z0, y0);
    line->SetPoint(1, z1, y1);
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    line->Draw();
    m_SdcInYZ_line.push_back(line);
  }
  m_canvas_vertex->Update();
#endif

#if CATCH
  double zc0 = -150, zc1 = 400;
  double xc0 = tp->GetX( zc0 ), yc0 = tp->GetY( zc0 );
  double xc1 = tp->GetX( zc1 ), yc1 = tp->GetY( zc1 );
  {
    TPolyLine3D *p = new TPolyLine3D(2);
    p->SetLineColor(kRed);
    p->SetLineWidth(1);
    p->SetPoint( 0, xc0, yc0, zc0 );
    p->SetPoint( 1, xc1, yc1, zc1 );
    m_SdcInTrack_Catch_cont.push_back(p);

    TPolyLine *lxy = new TPolyLine(2);
    lxy->SetPoint( 0, xc0, yc0 );
    lxy->SetPoint( 1, xc1, yc1 );
    lxy->SetLineColor(kRed);
    lxy->SetLineWidth(1);
    m_SdcInTrack_Catch_xy_cont.push_back(lxy);

    TPolyLine *lzx = new TPolyLine(2);
    lzx->SetPoint( 0, zc0, xc0 );
    lzx->SetPoint( 1, zc1, xc1 );
    lzx->SetLineColor(kRed);
    lzx->SetLineWidth(1);
    m_SdcInTrack_Catch_zx_cont.push_back(lzx);

    TPolyLine *lzy = new TPolyLine(2);
    lzy->SetPoint( 0, zc0, yc0 );
    lzy->SetPoint( 1, zc1, yc1 );
    lzy->SetLineColor(kRed);
    lzy->SetLineWidth(1);
    m_SdcInTrack_Catch_zy_cont.push_back(lzy);


  }
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawSdcOutLocalTrack(const DCLocalTrack *tp)
{
#if SdcOut
  Double_t z0 = gGeom.GetLocalZ(IdSDC2V);
  Double_t x0 = tp->GetX(z0), y0 = tp->GetY(z0);
  Double_t z0_X = gGeom.GetLocalZ(IdSDC2X);

  ThreeVector gPos0 = gGeom.Local2GlobalPos(IdSDC2X , TVector3(x0, y0, z0-z0_X));

  //Double_t z1 = gGeom.GetLocalZ(IdRKINIT);
  //Double_t z1 = gGeom.GetLocalZ(IdTOF);
  Double_t z_offset = 1500.;
  Double_t z1 = gGeom.GetLocalZ(IdSDC3X)+z_offset;
  Double_t x1 = tp->GetX(z1), y1 = tp->GetY(z1);
  //ThreeVector gPos1 = gGeom.Local2GlobalPos(IdTOF,  TVector3(x1, y1, 0.));
  ThreeVector gPos1 = gGeom.Local2GlobalPos(IdSDC3X,  TVector3(x1, y1, z_offset));

  auto p = new TPolyLine3D(2);
  p->SetLineColor(kRed);
  p->SetLineWidth(1);
  p->SetPoint(0, gPos0.x(), gPos0.y(), gPos0.z());
  p->SetPoint(1, gPos1.x(), gPos1.y(), gPos1.z());
  m_SdcOutTrack.push_back(p);
  m_canvas->cd(1)->cd(2);
  p->Draw();
  gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawVertex(const ThreeVector& vertex)
{
#if Vertex
  del::DeleteObject(m_VertexPointXZ);
  del::DeleteObject(m_VertexPointYZ);

  Double_t x = vertex.x();
  Double_t y = vertex.y();
  Double_t z = vertex.z();

  m_VertexPointXZ = new TMarker(z, x, 34);
  m_VertexPointXZ->SetMarkerSize(1);
  m_VertexPointXZ->SetMarkerColor(kBlue);
  m_VertexPointYZ = new TMarker(z, y, 34);
  m_VertexPointYZ->SetMarkerSize(1);
  m_VertexPointYZ->SetMarkerColor(kBlue);

  m_canvas_vertex->cd(1);
  m_VertexPointXZ->Draw();
  m_canvas_vertex->cd(2);
  m_VertexPointYZ->Draw();
  m_canvas_vertex->Update();
  hddaq::cout <<"Draw Vertex! " << x << " " << y << " " << z << std::endl;
#endif

}

//_____________________________________________________________________________
void
EventDisplay::DrawMissingMomentum(const ThreeVector& mom, const ThreeVector& pos)
{
#if Vertex
  if(std::abs(pos.x())>40.) return;
  if(std::abs(pos.y())>20.) return;
  if(std::abs(pos.z())>20.) return;

  Double_t z0 = pos.z();
  Double_t x0 = pos.x();
  Double_t y0 = pos.y();
  Double_t u0 = mom.x()/mom.z();
  Double_t v0 = mom.y()/mom.z();
  Double_t z1 = MaxZ;
  Double_t x1 = x0 + u0*z1;
  Double_t y1 = y0 + v0*z1;
  m_canvas_vertex->cd(1);
  m_MissMomXZ_line = new TPolyLine(2);
  m_MissMomXZ_line->SetPoint(0, z0, x0);
  m_MissMomXZ_line->SetPoint(1, z1, x1);
  m_MissMomXZ_line->SetLineColor(kBlue);
  m_MissMomXZ_line->SetLineWidth(1);
  m_MissMomXZ_line->Draw();
  m_canvas_vertex->cd(2);
  m_MissMomYZ_line = new TPolyLine(2);
  m_MissMomYZ_line->SetPoint(0, z0, y0);
  m_MissMomYZ_line->SetPoint(1, z1, y1);
  m_MissMomYZ_line->SetLineColor(kBlue);
  m_MissMomYZ_line->SetLineWidth(1);
  m_MissMomYZ_line->Draw();
  m_canvas_vertex->Update();
  ::sleep(3);
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawHypsTrack(Int_t nStep, const std::vector<TVector3>& StepPoint,
			   Double_t q)
{

#if DrawOneHypsTrack
  del::DeleteObject(m_Hyps_step_mark);
#endif

  auto step_mark = new TPolyMarker3D(nStep);
  for(Int_t i=0; i<nStep; ++i){
    step_mark->SetPoint(i,
			StepPoint[i].x(),
			StepPoint[i].y(),
			StepPoint[i].z());
  }

  Color_t color = (q > 0) ? kRed : kBlue;

  if(!m_is_save_mode){
    step_mark->SetMarkerSize(1);
    step_mark->SetMarkerStyle(6);
  }
  step_mark->SetMarkerColor(color);

  m_Hyps_step_mark.push_back(step_mark);

  m_canvas->cd(1)->cd(2);
  step_mark->Draw();
  m_canvas->Update();

#if Vertex
  del::DeleteObject(m_HypsMarkVertexX);
  m_HypsMarkVertexX = new TPolyMarker(nStep);
  for(Int_t i=0; i<nStep; ++i){
    Double_t x = StepPoint[i].x()-BeamAxis;
    Double_t z = StepPoint[i].z()-zTarget;
    m_HypsMarkVertexX->SetPoint(i, z, x);
  }
  m_HypsMarkVertexX->SetMarkerSize(0.4);
  m_HypsMarkVertexX->SetMarkerColor(color);
  m_HypsMarkVertexX->SetMarkerStyle(6);
  m_canvas_vertex->cd(1);
  m_HypsMarkVertexX->Draw();
  del::DeleteObject(m_HypsMarkVertexY);
  m_HypsMarkVertexY = new TPolyMarker(nStep);
  for(Int_t i=0; i<nStep; ++i){
    Double_t y = StepPoint[i].y();
    Double_t z = StepPoint[i].z()-zTarget;
    m_HypsMarkVertexY->SetPoint(i, z, y);
  }
  m_HypsMarkVertexY->SetMarkerSize(0.4);
  m_HypsMarkVertexY->SetMarkerColor(color);
  m_HypsMarkVertexY->SetMarkerStyle(6);
  m_canvas_vertex->cd(2);
  m_HypsMarkVertexY->Draw();
  m_canvas_vertex->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawHypsTrackToLast(Int_t nStep, const std::vector<TVector3>& StepPoint,
				 Double_t q)
{

#if DrawOneHypsTrack
  del::DeleteObject(m_Hyps_step_mark_tolast);
#endif

  auto step_mark = new TPolyMarker3D(nStep);
  for(Int_t i=0; i<nStep; ++i){
    step_mark->SetPoint(i,
			StepPoint[i].x(),
			StepPoint[i].y(),
			StepPoint[i].z());
  }

  Color_t color = (q > 0) ? kRed : kBlue;

  if(!m_is_save_mode){
    step_mark->SetMarkerSize(1);
    step_mark->SetMarkerStyle(6);
  }
  step_mark->SetMarkerColor(color);

  m_Hyps_step_mark_tolast.push_back(step_mark);

  m_canvas->cd(1)->cd(2);
  step_mark->Draw();
  m_canvas->Update();

#if Vertex
  del::DeleteObject(m_HypsMarkVertexX);
  m_HypsMarkVertexX = new TPolyMarker(nStep);
  for(Int_t i=0; i<nStep; ++i){
    Double_t x = StepPoint[i].x()-BeamAxis;
    Double_t z = StepPoint[i].z()-zTarget;
    m_HypsMarkVertexX->SetPoint(i, z, x);
  }
  m_HypsMarkVertexX->SetMarkerSize(0.4);
  m_HypsMarkVertexX->SetMarkerColor(color);
  m_HypsMarkVertexX->SetMarkerStyle(6);
  m_canvas_vertex->cd(1);
  m_HypsMarkVertexX->Draw();
  del::DeleteObject(m_HypsMarkVertexY);
  m_HypsMarkVertexY = new TPolyMarker(nStep);
  for(Int_t i=0; i<nStep; ++i){
    Double_t y = StepPoint[i].y();
    Double_t z = StepPoint[i].z()-zTarget;
    m_HypsMarkVertexY->SetPoint(i, z, y);
  }
  m_HypsMarkVertexY->SetMarkerSize(0.4);
  m_HypsMarkVertexY->SetMarkerColor(color);
  m_HypsMarkVertexY->SetMarkerStyle(6);
  m_canvas_vertex->cd(2);
  m_HypsMarkVertexY->Draw();
  m_canvas_vertex->Update();
#endif
}


//_____________________________________________________________________________
void EventDisplay::ShowHitFiber(Int_t layer, Int_t segment, Double_t pe, Double_t ctime)// const
{
  Color_t colorPallet[5] = {kAzure, kTeal, kSpring, kOrange, kPink};
  Color_t color = kBlack;

  Double_t p0[8] = {-0.371816, 0.996831, -1.594047, -0.366108, -1.736902, -0.362409, -1.585484, -0.196685};
  Double_t p1[8] = {0.007238, 0.004976, 0.009632, 0.007535, 0.010529, 0.007181, 0.010438, 0.006284};

  Double_t z_time = (ctime-p0[layer])/p1[layer];

  if (!(pe > 0))
    pe = -1;

  if (pe <= 0)
    color = kGray;
  else {
    Int_t color_unit = 200;
    Int_t color_index = ((Int_t)pe/color_unit);
    Int_t sub_color = ((Int_t)pe%color_unit)/(color_unit/10);
    if (color_index>=5) {
      color_index = 4;
      sub_color = 10;
    }
    color = colorPallet[color_index] + sub_color ;

  }


  if (segment>=0 && segment<NumOfSegCFT[layer]) {
#if CATCH
    Double_t x, y;
    FiberPosPhi(layer, segment, &x, &y);
    //hddaq::cout << x << ", " << y << ", adcLow = " << pe << std::endl;
    m_hbase_catch->Fill(x, y, pe);

    m_CFT_Arc_cont[layer][segment]->SetLineColor(kRed);
#endif

#if CATCH3d    

    if (layer%2 == 1) {
      // phi layer
      std::string node_name;
      node_name = Form( "CFT%d_Node_%d", layer/2, segment );

      TNode *node = m_geometry_catch->GetNode( node_name.c_str() );
      if( !node ){
	hddaq::cout << "#E EventDisplay::ShowHitFiber() " 
		    << "no such node : " << node_name << std::endl;
	return;
      }

      node->SetVisibility(1);
      node->SetLineColor(color);

      if (z_time>-500 && z_time<1000) {
	Int_t z_min =(Int_t) (z_time-100);
	if (z_min<-50)
	  z_min = -50;
	Int_t z_max = (Int_t) (z_time+100);
	if (z_max>450)
	  z_max = 450;

	Int_t n = ((Int_t)z_max - (Int_t)z_min)/5 + 1;
	TPolyMarker3D *z_pos = new TPolyMarker3D( n );

	for (Int_t in=0; in<n; in++) {
	  Double_t iz = z_min + (Double_t)in*5;
	  z_pos->SetPoint( in, x, y, (Double_t)iz - offsetCATCH);

	  z_pos->SetMarkerSize(0.5);
	  z_pos->SetMarkerColor(kPink);
	  z_pos->SetMarkerStyle(20);
	}
	m_CFT_UV_cont.push_back(z_pos);

      }
    } else {
      // UV layer
      const DCGeomMan & geomMan=DCGeomMan::GetInstance();

      Int_t lnum=301+layer;

      Double_t R     = geomMan.GetLocalZ(lnum);
      if (segment%2 == 0)
	R -= 0.4755/2;
      else
	R += 0.4755/2;
      Double_t SegNumUV=NumOfSegCFT[layer];
      Double_t Phi0 = -(360./SegNumUV)*(Double_t)segment + 90.;
      if (layer == CFT_V2 || layer == CFT_V4) {
	Phi0 = (360./SegNumUV)*(Double_t)segment + 90.;
      }

      Double_t slope = 0.;
      Double_t d_phi = 5.0;
      Int_t nStep = (Int_t)(360/d_phi);
      if(layer==0||layer==4){
	slope = 400. /360.;
      }else if(layer==2||layer==6){
	slope = -400. /360.;
	d_phi  *= -1.;
      }

      TPolyMarker3D *uv_fiber = new TPolyMarker3D( nStep );
      for (Int_t i=0; i<nStep; ++i) {
	Double_t phi = Phi0 + d_phi*(Double_t)i;
	Double_t z = slope*(phi-Phi0);
	Double_t x = R * cos(phi * TMath::DegToRad());
	Double_t y = R * sin(phi * TMath::DegToRad());
	uv_fiber->SetPoint( i, x, y, z - offsetCATCH);
      }
      uv_fiber->SetMarkerSize(10);
      uv_fiber->SetMarkerColor(color);
      uv_fiber->SetMarkerStyle(6);
      m_CFT_UV_cont.push_back(uv_fiber);


      if (z_time>-500 && z_time<1000) {
	Int_t z_min =(Int_t) (z_time-100);
	if (z_min<-50)
	  z_min = -50;
	Int_t z_max = (Int_t) (z_time+100);
	if (z_max>450)
	  z_max = 450;

	Int_t n = ((Int_t)z_max - (Int_t)z_min)/5 + 1;

	TPolyMarker3D *z_pos = new TPolyMarker3D( n );
	for (Int_t in=0; in<n; in++) {
	  Double_t iz = z_min + (Double_t)in*5;
	  Double_t phi = (Double_t)iz/slope + Phi0;
	  Double_t x = R * cos(phi * TMath::DegToRad());
	  Double_t y = R * sin(phi * TMath::DegToRad());

	  z_pos->SetPoint( in, x, y, (Double_t)iz - offsetCATCH);
	  z_pos->SetMarkerSize(0.5);
	  z_pos->SetMarkerColor(kPink);
	  z_pos->SetMarkerStyle(20);
	}
	m_CFT_UV_cont.push_back(z_pos);
      }
    }
#endif
  }
}

//______________________________________________________________________________
void EventDisplay::ShowHitFiberTracked(Int_t layer, Int_t segment, Double_t z, Bool_t flagProton)// const
{
  static const std::string func_name("[EventDisplay::ShowHitFiberTracked()]");

  //printf("ShowHitFiber : layer %d, seg %d, pe %f\n", layer, segment, pe);

  Color_t colorPallet[5] = {kAzure, kTeal, kSpring, kOrange, kPink};
  Color_t color = kGreen;
  if (m_CFTTrack_cont.size()==2)
    color = kYellow;
  else if (m_CFTTrack_cont.size()==3)
    color = kOrange;
  else if (m_CFTTrack_cont.size()==4)
    color = kPink;


#if CATCH3d
  if (segment>=0 && segment<NumOfSegCFT[layer]) {
    double x, y;
    FiberPosPhi(layer, segment, &x, &y);

    if (layer%2 == 1) {
      // phi layer
      if (z>-1000 && z<1000) {
	int z_min =(int) (z-10);
	if (z_min<-200)
	  z_min = -50;
	int z_max = (int) (z+10);
	if (z_max>300)
	  z_max = 300;

	int n = (int)z_max - (int)z_min + 1;
	TPolyMarker3D *z_pos = new TPolyMarker3D( n );

	for (int in=0; in<n; in++) {
	  double iz = z_min + (double)in;
	  z_pos->SetPoint( in, x, y, (double)iz);      

	  z_pos->SetMarkerSize(0.5);
	  z_pos->SetMarkerColor(color);
	  if (flagProton)
	    z_pos->SetMarkerColor(kBlue);

	  z_pos->SetMarkerStyle(20);
	}
	m_CFT_UV_cont.push_back(z_pos);

      }
    } else {
      // UV layer
      const DCGeomMan & geomMan=DCGeomMan::GetInstance();

      int lnum=301+layer;

      double R     = geomMan.GetLocalZ(lnum);
      if (segment%2 == 0)
	R -= 0.4755/2;
      else
	R += 0.4755/2;
      double SegNumUV=NumOfSegCFT[layer];
      double Phi0 = -(360./SegNumUV)*(double)segment + 90.;
      if (layer == CFT_V2 || layer == CFT_V4) {
	Phi0 = (360./SegNumUV)*(double)segment + 90.;
      }

      double slope = 0.;
      double d_phi = 5.0;
      int nStep = (int)(360/d_phi);
      if(layer==0||layer==4){
	slope = 400. /360.;
      }else if(layer==2||layer==6){
	slope = -400. /360.;
	d_phi  *= -1.;
      }

      if (z>-200 && z<300) {
	int z_min =(int) (z-10);
	if (z_min<-200)
	  z_min = -50;
	int z_max = (int) (z+10);
	if (z_max>300)
	  z_max = 300;

	int n = (int)z_max - (int)z_min + 1;

	TPolyMarker3D *z_pos = new TPolyMarker3D( n + 10);
	for (int in=0; in<n; in++) {
	  double iz = z_min + (double)in;
	  double phi = (double)(iz + offsetCATCH)/slope + Phi0;
	  double x = R * cos(phi * TMath::DegToRad());
	  double y = R * sin(phi * TMath::DegToRad());

	  z_pos->SetPoint( in, x, y, (double)iz);      
	  z_pos->SetMarkerSize(0.5);
	  z_pos->SetMarkerColor(color);
	  if (flagProton)
	    z_pos->SetMarkerColor(kBlue);

	  z_pos->SetMarkerStyle(20);
	}

	for (int in=0; in<10; in++) {
	  double iz = (double)in * 40;
	  double phi = (double)iz/slope + Phi0;
	  double x = R * cos(phi * TMath::DegToRad());
	  double y = R * sin(phi * TMath::DegToRad());

	  z_pos->SetPoint( in, x, y, (double)iz - offsetCATCH);      
	  z_pos->SetMarkerSize(0.5);
	  z_pos->SetMarkerColor(color);
	  if (flagProton)
	    z_pos->SetMarkerColor(kBlue);

	  z_pos->SetMarkerStyle(20);
	}
	m_CFT_UV_cont.push_back(z_pos);
      }
    }
  }
#endif
}

//______________________________________________________________________________
void EventDisplay::ShowHitBGO(Int_t segment, Double_t de) const
{
  Color_t colorPallet[5] = {kAzure, kTeal, kSpring, kOrange, kPink};
  Color_t color = kBlack;

  if (de <= 0)
    color = kGray;
  else {
    Int_t color_unit = 5000;
    Int_t color_index = ((Int_t)de/color_unit);
    Int_t sub_color = ((Int_t)de%color_unit)/(color_unit/10);
    if (color_index>=5) {
      color_index = 4;
      sub_color = 10;
    }
    color = colorPallet[color_index] + sub_color ;

  }

#if CATCH
  Int_t size = m_BGO_Line_cont[segment].size();
  for (Int_t i=0; i<size; i++)
    m_BGO_Line_cont[segment][i]->SetLineColor(kRed);

  Double_t x, y;
  BGOPos(segment, &x, &y);
  m_hbase_catch->Fill(x, y, de);
#endif
  
#if CATCH3d
  std::string node_name;
  node_name = Form( "BGOseg_node_%d", segment );

  TNode *node = m_geometry_catch->GetNode( node_name.c_str() );
  if( !node ){
    hddaq::cout << "#E EventDisplay::ShowHitBGO() "
		<< "no such node : " << node_name << std::endl;
    return;
  }

  node->SetVisibility(1);
  node->SetLineColor(color);
#endif
}

//______________________________________________________________________________
void
EventDisplay::SetBGOWaveformCanvas(Int_t nhit )
{
#if BGO_WF
  m_canvas_hist9->Clear();

  if (nhit == 2)
    m_canvas_hist9->Divide(2, 1);
  else if (nhit>= 3 && nhit<=4)
    m_canvas_hist9->Divide(2, 2);
  else if (nhit>= 5 && nhit<=6)
    m_canvas_hist9->Divide(2, 3);
  else if (nhit>= 7 && nhit<=8)
    m_canvas_hist9->Divide(2, 4);
  else if (nhit>= 9 && nhit<=12)
    m_canvas_hist9->Divide(3, 4);
  else if (nhit>= 13 && nhit<=16)
    m_canvas_hist9->Divide(4, 4);
  else if (nhit>= 17 && nhit<=20)
    m_canvas_hist9->Divide(5, 4);
  else if (nhit>= 21 && nhit<=24)
    m_canvas_hist9->Divide(6, 4);


#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawBGOWaveform(Int_t nc, Int_t ngraph, Int_t seg, TGraphErrors* gr )
{
#if BGO_WF
  m_canvas_hist9->cd(nc);

  gr->SetNameTitle(Form("gra_%d_%d", seg, ngraph),
		  Form("Segment %d", seg));

  if (ngraph == 0) {
    gr->SetMarkerSize(1.);
    gr->SetMarkerStyle(20);
    gr->Draw("ap");
  } else {
    gr->SetMarkerSize(1.);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kRed);
    gr->Draw("p");

  }

#endif
}


//______________________________________________________________________________
void
EventDisplay::DrawBGOFitFunc(Int_t nc, Int_t seg, TF1* func )
{
#if BGO_WF
  m_canvas_hist9->cd(nc);

  func->Draw("same");

#endif
}

//______________________________________________________________________________
void
EventDisplay::SetTagWaveformCanvas(Int_t nhit )
{
#if TAG_WF
  m_canvas_hist10->Clear();

  if (nhit == 2)
    m_canvas_hist10->Divide(2, 1);
  else if (nhit>= 3 && nhit<=4)
    m_canvas_hist10->Divide(2, 2);
  else if (nhit>= 5 && nhit<=6)
    m_canvas_hist10->Divide(2, 3);
  else if (nhit>= 7 && nhit<=8)
    m_canvas_hist10->Divide(2, 4);
  else if (nhit>= 9 && nhit<=12)
    m_canvas_hist10->Divide(3, 4);

#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawTagWaveform(Int_t nc, Int_t ngraph, Int_t seg, TGraphErrors* gr )
{
#if TAG_WF
  m_canvas_hist10->cd(nc);

  gr->SetNameTitle(Form("gra_%d_%d", seg, ngraph),
		  Form("Segment %d", seg));

  if (ngraph == 0) {
    gr->SetMarkerSize(1.);
    gr->SetMarkerStyle(20);
    gr->Draw("ap");
  } else {
    gr->SetMarkerSize(1.);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kRed);
    gr->Draw("p");

  }

#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawTagFitFunc(Int_t nc, Int_t seg, TF1* func )
{
#if TAG_WF
  m_canvas_hist10->cd(nc);

  func->Draw("same");

#endif
}

//______________________________________________________________________________
void
EventDisplay::DrawCFT_Time( Int_t layer, Int_t seg, Int_t LorT, Double_t time )
{
#if CATCH_Timing

  TH2 *hp=0;

  if (layer == 0 && LorT == 0) {
    hp = m_hist_cft1_l;
  } else if (layer == 0 && LorT == 1) {
    hp = m_hist_cft1_t;
  } else if (layer == 1 && LorT == 0) {
    hp = m_hist_cft2_l;
  } else if (layer == 1 && LorT == 1) {
    hp = m_hist_cft2_t;
  } else if (layer == 2 && LorT == 0) {
    hp = m_hist_cft3_l;
  } else if (layer == 2 && LorT == 1) {
    hp = m_hist_cft3_t;
  } else if (layer == 3 && LorT == 0) {
    hp = m_hist_cft4_l;
  } else if (layer == 3 && LorT == 1) {
    hp = m_hist_cft4_t;
  } else if (layer == 4 && LorT == 0) {
    hp = m_hist_cft5_l;
  } else if (layer == 4 && LorT == 1) {
    hp = m_hist_cft5_t;
  } else if (layer == 5 && LorT == 0) {
    hp = m_hist_cft6_l;
  } else if (layer == 5 && LorT == 1) {
    hp = m_hist_cft6_t;
  } else if (layer == 6 && LorT == 0) {
    hp = m_hist_cft7_l;
  } else if (layer == 6 && LorT == 1) {
    hp = m_hist_cft7_t;
  } else if (layer == 7 && LorT == 0) {
    hp = m_hist_cft8_l;
  } else if (layer == 7 && LorT == 1) {
    hp = m_hist_cft8_t;
  }

  hp->Fill( seg, time );
#endif

}

//______________________________________________________________________________
void
EventDisplay::DrawCFT_AdcCor( Int_t layer, Int_t seg, Int_t HorL, Double_t adccor )
{
#if CATCH_ADC

  TH2 *hp=0;

  if (layer == 0 && HorL == 0) {
    hp = m_hist_cft1_hi;
  } else if (layer == 0 && HorL == 1) {
    hp = m_hist_cft1_lo;
  } else if (layer == 1 && HorL == 0) {
    hp = m_hist_cft2_hi;
  } else if (layer == 1 && HorL == 1) {
    hp = m_hist_cft2_lo;
  } else if (layer == 2 && HorL == 0) {
    hp = m_hist_cft3_hi;
  } else if (layer == 2 && HorL == 1) {
    hp = m_hist_cft3_lo;
  } else if (layer == 3 && HorL == 0) {
    hp = m_hist_cft4_hi;
  } else if (layer == 3 && HorL == 1) {
    hp = m_hist_cft4_lo;
  } else if (layer == 4 && HorL == 0) {
    hp = m_hist_cft5_hi;
  } else if (layer == 4 && HorL == 1) {
    hp = m_hist_cft5_lo;
  } else if (layer == 5 && HorL == 0) {
    hp = m_hist_cft6_hi;
  } else if (layer == 5 && HorL == 1) {
    hp = m_hist_cft6_lo;
  } else if (layer == 6 && HorL == 0) {
    hp = m_hist_cft7_hi;
  } else if (layer == 6 && HorL == 1) {
    hp = m_hist_cft7_lo;
  } else if (layer == 7 && HorL == 0) {
    hp = m_hist_cft8_hi;
  } else if (layer == 7 && HorL == 1) {
    hp = m_hist_cft8_lo;
  }


  hp->Fill( seg, adccor );

#endif
}

//______________________________________________________________________________
void EventDisplay::ShowHitPiID(Int_t segment)
{

#if CATCH
  int size = m_PiID_Line_cont[segment].size();
  for (int i=0; i<size; i++)
    m_PiID_Line_cont[segment][i]->SetLineColor(kRed);
#endif

#if CATCH3d  
  std::string node_name;
  node_name = Form( "PiIDseg_node_%d", segment );

  TNode *node = m_geometry_catch->GetNode( node_name.c_str() );
  if( !node ){
    hddaq::cout << "#E EventDisplay::ShowHitPiID() "
		<< "no such node : " << node_name << std::endl;
    return;
  }

  node->SetVisibility(1);
  node->SetLineColor(kRed);
#endif

}


//______________________________________________________________________________
void
EventDisplay::DrawCFTLocalTrack( const CFTLocalTrack *tp, Bool_t flagP , Int_t k_color)
{
#if CATCH
  ThreeVector Pos0 = tp->GetPos0();
  ThreeVector Dir = tp->GetDir();

  {
    Int_t color = kGreen;
    if (m_CFTTrack_cont.size()==1)
      color = kYellow;
    else if (m_CFTTrack_cont.size()==2)
      color = kOrange;
    else if (m_CFTTrack_cont.size()==3)
      color = kPink;


    if (flagP)
      color = kBlue;

    if (k_color == 1) {
      color = kCyan;
    } else if (k_color == 2) {
      color = kRed; // Lambda
    }

    TPolyLine3D *p = new TPolyLine3D(2);
    p->SetLineColor(color);
    p->SetLineWidth(2);


    ThreeVector pos1 = Pos0 - 3.0*Dir;
    p->SetPoint( 0, pos1.x(), pos1.y(), pos1.z() );
    ThreeVector pos2 = Pos0 + 2.5*Dir;
    p->SetPoint( 1, pos2.x(), pos2.y(), pos2.z() );
    m_CFTTrack_cont.push_back(p);

    TPolyLine *lxy = new TPolyLine(2);
    lxy->SetPoint( 0, pos1.x(), pos1.y() );
    lxy->SetPoint( 1, pos2.x(), pos2.y() );
    lxy->SetLineColor(color);
    lxy->SetLineWidth(1);
    m_CFTTrack_xy_cont.push_back(lxy);

    TPolyLine *lzx = new TPolyLine(2);
    lzx->SetPoint( 0, pos1.z(), pos1.x() );
    lzx->SetPoint( 1, pos2.z(), pos2.x() );
    lzx->SetLineColor(color);
    lzx->SetLineWidth(1);
    m_CFTTrack_zx_cont.push_back(lzx);

    TPolyLine *lzy = new TPolyLine(2);
    lzy->SetPoint( 0, pos1.z(), pos1.y() );
    lzy->SetPoint( 1, pos2.z(), pos2.y() );
    lzy->SetLineColor(color);
    lzy->SetLineWidth(1);
    m_CFTTrack_zy_cont.push_back(lzy);

    /*
    Double_t slope_xy = tp->GetAxy();
    Double_t theta_xy = TMath::Atan(slope_xy)*TMath::DegToRad();
    for (Int_t i=0; i<m_CFTTrackCand_zx_cont.size(); i++) {
      Double_t theta_xy1 = TMath::Atan(m_CFTTrackCandSlope_xy_cont1[i])*TMath::DegToRad();
      if (std::abs(theta_xy-theta_xy1)<5)
	m_CFTTrackCand_zx_cont[i]->SetMarkerColor(color);
    }

    for (Int_t i=0; i<m_CFTTrackCand_zy_cont.size(); i++) {
      Double_t theta_xy2 = TMath::Atan(m_CFTTrackCandSlope_xy_cont2[i])*TMath::DegToRad();
      if (std::abs(theta_xy-theta_xy2)<5)
	m_CFTTrackCand_zy_cont[i]->SetMarkerColor(color);
    }
    */
  }

#endif

}

//______________________________________________________________________________
void
EventDisplay::DrawCFTLocalTrack_dE_E( CFTParticle *CFTPart, bool flagP )
{


  const CFTLocalTrack *tp = CFTPart->GetTrack();

  ThreeVector Pos0 = tp->GetPos0();
  ThreeVector Dir = tp->GetDir();

  {
    int color = kGreen;
    if (m_CFTTrack_cont.size()==1)
      color = kYellow;
    else if (m_CFTTrack_cont.size()==2)
      color = kOrange;
    else if (m_CFTTrack_cont.size()==3)
      color = kPink;


    if (flagP)
      color = kBlue;

    ThreeVector pos1 = Pos0 - 3.0*Dir;
    ThreeVector pos2 = Pos0 + 2.5*Dir;
    
#if CATCH3d    
    TPolyLine3D *p = new TPolyLine3D(2);
    p->SetLineColor(color);
    p->SetLineWidth(2);

    p->SetPoint( 0, pos1.x(), pos1.y(), pos1.z() );
    p->SetPoint( 1, pos2.x(), pos2.y(), pos2.z() );
    m_CFTTrack_cont.push_back(p);
#endif

#if CATCH    
    TPolyLine *lxy = new TPolyLine(2);
    lxy->SetPoint( 0, pos1.x(), pos1.y() );
    lxy->SetPoint( 1, pos2.x(), pos2.y() );
    lxy->SetLineColor(color);
    lxy->SetLineWidth(1);
    m_CFTTrack_xy_cont.push_back(lxy);

    TPolyLine *lzx = new TPolyLine(2);
    lzx->SetPoint( 0, pos1.z(), pos1.x() );
    lzx->SetPoint( 1, pos2.z(), pos2.x() );
    lzx->SetLineColor(color);
    lzx->SetLineWidth(1);
    m_CFTTrack_zx_cont.push_back(lzx);

    TPolyLine *lzy = new TPolyLine(2);
    lzy->SetPoint( 0, pos1.z(), pos1.y() );
    lzy->SetPoint( 1, pos2.z(), pos2.y() );
    lzy->SetLineColor(color);
    lzy->SetLineWidth(1);
    m_CFTTrack_zy_cont.push_back(lzy);

    /*
    double slope_xy = tp->GetAxy();
    double theta_xy = atan(slope_xy)*math::Rad2Deg();
    for (int i=0; i<m_CFTTrackCand_zx_cont.size(); i++) {
      double theta_xy1 = atan(m_CFTTrackCandSlope_xy_cont1[i])*math::Rad2Deg();
      if (fabs(theta_xy-theta_xy1)<5)
	m_CFTTrackCand_zx_cont[i]->SetMarkerColor(color);
    }

    for (int i=0; i<m_CFTTrackCand_zy_cont.size(); i++) {
      double theta_xy2 = atan(m_CFTTrackCandSlope_xy_cont2[i])*math::Rad2Deg();
      if (fabs(theta_xy-theta_xy2)<5)
	m_CFTTrackCand_zy_cont[i]->SetMarkerColor(color);
    }
    */
    int nhPhi   = tp->GetNHit();
    int nhUV = tp->GetNHitUV();
    double Total_dEphi_max = tp->GetTotalMaxdEphi();
    double Total_dEuv_max  = tp->GetTotalMaxdEuv ();
    double theta = tp->GetThetaCFT();

    double dE = Total_dEphi_max*sin(theta*TMath::DegToRad())/(double)(nhPhi)
      + Total_dEuv_max/(double)(nhUV);
    double bgo_energy = CFTPart->GetBGOEnergy();
    int    bgo_seg = CFTPart->GetTrackBGOSeg();
    int    piid_seg = CFTPart->GetTrackPiIDSeg();
    std::cout << "Total_dEphi_max : " << Total_dEphi_max
	      << ", nhPhi : " << nhPhi
	      << ", Total_dEuv_max : " << Total_dEuv_max
	      << ", nhUV : " << nhUV
	      << std::endl;
    std::cout << "BGO ( " << bgo_seg << " ) : " << bgo_energy << ", dE : " << dE << std::endl;
    std::cout << "PiID ( " << piid_seg << " ) " << std::endl;

    m_hist_dE_E->Fill(bgo_energy, dE);

    TGraph *p_gr = new TGraph(1, &bgo_energy, &dE);
    p_gr->SetMarkerColor(color);
    p_gr->SetMarkerStyle(20);

    m_CATCH_dE_E_cont.push_back(p_gr);
#endif

  }


}

//______________________________________________________________________________
void EventDisplay::ShowHitTagger(const TString& name, Int_t segment, Double_t de) const
{
  Color_t colorPallet[5] = {kAzure, kTeal, kSpring, kOrange, kPink};
  Color_t color = kBlack;

  if (de <= 0)
    color = kGray;
  else {
    Int_t color_unit = 5000;
    Int_t color_index = ((Int_t)de/color_unit);
    Int_t sub_color = ((Int_t)de%color_unit)/(color_unit/10);
    if (color_index>=5) {
      color_index = 4;
      sub_color = 10;
    }
    color = colorPallet[color_index] + sub_color ;

  }

  double TagSF_X = 6.;
  double TagSF_Y = 1.;

  double TagPL_X = 3.;
  double TagPL_Y = 7.4;
  double overlap = 2.7;

  if (name == "TAG-PL") {
    Int_t size = m_TagPL_Line_cont[segment].size();
    for (Int_t i=0; i<size; i++)
      m_TagPL_Line_cont[segment][i]->SetLineColor(kRed);

    double x0 = 24.0 + 10*(segment%2);
    double y0 = segment*(TagPL_Y - overlap);
    m_hbase_tagger->Fill(x0, y0, de);

  } else if (name == "TAG-SFF" ) {
    Int_t size = m_TagSFF_Line_cont[segment].size();
    for (Int_t i=0; i<size; i++)
      m_TagSFF_Line_cont[segment][i]->SetLineColor(kRed);

    double x0 = 0.0;
    double y0 = segment*TagSF_Y;
    m_hbase_tagger->Fill(x0, y0, de);

  } else if (name == "TAG-SFB" ) {
    Int_t size = m_TagSFB_Line_cont[segment].size();
    for (Int_t i=0; i<size; i++)
      m_TagSFB_Line_cont[segment][i]->SetLineColor(kRed);

    double x0 = 6.0;
    double y0 = 0.5 + segment*TagSF_Y;
    m_hbase_tagger->Fill(x0, y0, de);

  }

}

//_____________________________________________________________________________
void
EventDisplay::DrawTarget()
{
  if(!m_target_node){
    hddaq::cout << FUNC_NAME << " " << "node is null" << std::endl;
    return;
  }
  m_target_node->SetLineColor(kMagenta);
  m_canvas->cd(1)->cd(2);
  m_canvas->Update();
}


//_____________________________________________________________________________
void
EventDisplay::DrawRunEvent(Double_t xpos, Double_t ypos, const TString& arg)
{
  if(arg.Contains("Run")){
    // std::cout << arg << " find " << std::endl;
    m_canvas->cd(1)->cd(1)->Clear();
    m_canvas->cd(1)->cd(2)->Clear();
  }
  m_canvas->cd(1)->cd(1);
  TLatex tex;
  tex.SetTextAlign(12);
  tex.SetTextSize(0.5);
  // tex.SetTextSize(0.15);
  tex.SetNDC();
  tex.DrawLatex(xpos, ypos, arg);
  m_canvas->Update();
}

//_____________________________________________________________________________
void
EventDisplay::DrawText(Double_t xpos, Double_t ypos, const TString& arg)
{
  m_canvas->cd(1)->cd(2);
  TLatex tex;
  tex.SetTextAlign(12);
  tex.SetTextFont(42);
  tex.SetTextSize(0.025);
  tex.SetNDC();
  tex.DrawLatex(xpos, ypos, arg);
  m_canvas->Update();
}

//_____________________________________________________________________________
void
EventDisplay::EndOfEvent()
{
  del::DeleteObject(m_init_step_mark);
  del::DeleteObject(m_BcOutXZ_line);
  del::DeleteObject(m_BcOutYZ_line);
  del::DeleteObject(m_SdcInXZ_line);
  del::DeleteObject(m_SdcInYZ_line);
  del::DeleteObject(m_BcInTrack);
  del::DeleteObject(m_BcOutTrack);
  del::DeleteObject(m_BcOutTrack2);
  del::DeleteObject(m_BcOutTrack3);
  del::DeleteObject(m_SdcInTrack);
  del::DeleteObject(m_SdcInTrack2);
  del::DeleteObject(m_SdcOutTrack);
  del::DeleteObject(m_hs_step_mark);
  del::DeleteObject(m_Hyps_step_mark);
  del::DeleteObject(m_Hyps_step_mark_tolast);
  del::DeleteObject(m_VertexPointXZ);
  del::DeleteObject(m_VertexPointYZ);
  del::DeleteObject(m_MissMomXZ_line);
  del::DeleteObject(m_MissMomYZ_line);
  del::DeleteObject(m_HSMarkVertexXShs);
  del::DeleteObject(m_HypsMarkVertexXShs);
  del::DeleteObject(m_HypsMarkVertexX);
  del::DeleteObject(m_HypsMarkVertexY);
  ResetVisibility();
  ResetHist();

  ResetCATCH();
}

//_____________________________________________________________________________
void
EventDisplay::ResetVisibility(TNode *& node, Color_t c)
{
  if(!node) return;
  node->SetLineWidth(1);
  if(c==kWhite)
    node->SetVisibility(kFALSE);
  else
    node->SetLineColor(c);
}

//_____________________________________________________________________________
void
EventDisplay::ResetVisibility(std::vector<TNode*>& node, Color_t c)
{
  const std::size_t n = node.size();
  for(std::size_t i=0; i<n; ++i){
    ResetVisibility(node[i], c);
  }
}

//_____________________________________________________________________________
void
EventDisplay::ResetVisibility()
{
  //ResetVisibility(m_BC3x1_node);
  //ResetVisibility(m_BC3x2_node);
  //ResetVisibility(m_BC3v1_node);
  //ResetVisibility(m_BC3v2_node);
  //ResetVisibility(m_BC3u1_node);
  //ResetVisibility(m_BC3u2_node);
  //ResetVisibility(m_BC4u1_node);
  //ResetVisibility(m_BC4u2_node);
  //ResetVisibility(m_BC4v1_node);
  //ResetVisibility(m_BC4v2_node);
  //ResetVisibility(m_BC4x1_node);
  //ResetVisibility(m_BC4x2_node);
  ResetVisibility(m_SDC0x1_node);
  ResetVisibility(m_SDC0x2_node);
  ResetVisibility(m_SDC0u1_node);
  ResetVisibility(m_SDC0u2_node);
  ResetVisibility(m_SDC1x1_node);
  ResetVisibility(m_SDC1v1_node);
  ResetVisibility(m_SDC1u1_node);
  ResetVisibility(m_SDC1u2_node);
  ResetVisibility(m_SDC1x2_node);
  ResetVisibility(m_SDC1x3_node);
  ResetVisibility(m_SDC2v1_node);
  ResetVisibility(m_SDC2u1_node);
  ResetVisibility(m_SDC2u2_node);
  ResetVisibility(m_SDC2x1_node);
  ResetVisibility(m_SDC2x2_node);
  ResetVisibility(m_SDC3v1_node);
  ResetVisibility(m_SDC3u1_node);
  ResetVisibility(m_SDC3u2_node);
  ResetVisibility(m_SDC3x1_node);
  ResetVisibility(m_SDC3x2_node);
  //ResetVisibility(m_SDC4y1_node);
  //ResetVisibility(m_SDC4y2_node);
  //ResetVisibility(m_SDC4x1_node);
  //ResetVisibility(m_SDC4x2_node);
  //ResetVisibility(m_SDC5y1_node);
  //ResetVisibility(m_SDC5y2_node);
  //ResetVisibility(m_SDC5x1_node);
  //ResetVisibility(m_SDC5x2_node);
  //ResetVisibility(m_BH2seg_node, kBlack);
  ResetVisibility(m_TOFseg_node, kBlack);
  //ResetVisibility(m_WCseg_node, kBlack);
  //ResetVisibility(m_target_node, kBlack);

  for (int layer=0; layer<NumOfPlaneCFT/2; layer++)
    ResetVisibility( m_CFT_node_cont[layer] );

  ResetVisibility( m_BGOseg_node_cont );
  ResetVisibility( m_PiIDseg_node_cont );
  
}

//_____________________________________________________________________________
void
EventDisplay::ResetHist()
{

  for (int seg=0; seg<NumOfSegTagSF; seg++) {
    int size = m_TagSFF_Line_cont[seg].size();
    for (int i=0; i<size; i++)
      m_TagSFF_Line_cont[seg][i]->SetLineColor(kBlack);
  }
  for (int seg=0; seg<NumOfSegTagSF; seg++) {
    int size = m_TagSFB_Line_cont[seg].size();
    for (int i=0; i<size; i++)
      m_TagSFB_Line_cont[seg][i]->SetLineColor(kBlack);
  }
  for (int seg=0; seg<NumOfSegTagPL; seg++) {
    int size = m_TagPL_Line_cont[seg].size();
    for (int i=0; i<size; i++)
      m_TagPL_Line_cont[seg][i]->SetLineColor(kBlack);
  }
  m_hbase_tagger->Reset("ICES");

  //m_hist_aft_x->Reset("");
  //m_hist_aft_y->Reset("");

#if Hist_Timing
  m_hist_bh2->Reset();
  m_hist_sch->Reset();
  m_hist_tof->Reset();
  m_hist_sdc1->Reset();
  m_hist_sdc1p->Reset();

  m_hist_bc3->Reset();
  m_hist_bc3p->Reset();
  m_hist_bc3u->Reset();
  m_hist_bc3up->Reset();
  m_hist_bc3v->Reset();
  m_hist_bc3vp->Reset();

  m_hist_bc4->Reset();
  m_hist_bc4p->Reset();
  m_hist_bc4u->Reset();
  m_hist_bc4up->Reset();
  m_hist_bc4v->Reset();
  m_hist_bc4vp->Reset();

  m_hist_bc3_time->Reset();
  m_hist_bc3_time->SetMaximum(-1111);
  m_hist_bc3p_time->Reset();

  m_hist_bc4_time->Reset();
  m_hist_bc4_time->SetMaximum(-1111);
  m_hist_bc4p_time->Reset();
#endif

#if Hist_SdcOut
  m_hist_sdc3_l->Reset();
  m_hist_sdc3_t->Reset();
  m_hist_sdc3p_l->Reset();
  m_hist_sdc3p_t->Reset();

  m_hist_sdc3y_l->Reset();
  m_hist_sdc3y_t->Reset();
  m_hist_sdc3yp_l->Reset();
  m_hist_sdc3yp_t->Reset();

  m_hist_sdc4_l->Reset();
  m_hist_sdc4_t->Reset();
  m_hist_sdc4p_l->Reset();
  m_hist_sdc4p_t->Reset();

  m_hist_sdc4y_l->Reset();
  m_hist_sdc4y_t->Reset();
  m_hist_sdc4yp_l->Reset();
  m_hist_sdc4yp_t->Reset();

#endif

#if Hist_Timing
  m_hist_bh1->Reset();
  m_hist_bft->Reset();
  m_hist_bft_p->Reset();
  m_hist_bcIn->Reset();
  m_hist_bcOut->Reset();

  Int_t nc=m_BH1box_cont.size();
  for (Int_t i=0; i<nc; i++)
    m_BH1box_cont[i]->SetFillColor(kWhite);

  Int_t ncBh2=m_BH2box_cont.size();
  for (Int_t i=0; i<ncBh2; i++) {
    m_BH2box_cont[i]->SetFillColor(kWhite);
    m_BH2box_cont2[i]->SetFillColor(kWhite);
  }

  m_hist_bcOut_sdcIn->Reset();
  m_hist_sdcIn_predict->Reset();
  m_hist_sdcIn_predict2->Reset();
#endif

#if CATCH_Timing
  TH2 *hp_l=0, *hp_t=0;
  for (int layer=0; layer<8; layer++) {
    if (layer == 0) {
      hp_l = m_hist_cft1_l;
      hp_t = m_hist_cft1_t;
    } else if (layer == 1) {
      hp_l = m_hist_cft2_l;
      hp_t = m_hist_cft2_t;
    } else if (layer == 2) {
      hp_l = m_hist_cft3_l;
      hp_t = m_hist_cft3_t;
    } else if (layer == 3) {
      hp_l = m_hist_cft4_l;
      hp_t = m_hist_cft4_t;
    } else if (layer == 4) {
      hp_l = m_hist_cft5_l;
      hp_t = m_hist_cft5_t;
    } else if (layer == 5) {
      hp_l = m_hist_cft6_l;
      hp_t = m_hist_cft6_t;
    } else if (layer == 6) {
      hp_l = m_hist_cft7_l;
      hp_t = m_hist_cft7_t;
    } else if (layer == 7) {
      hp_l = m_hist_cft8_l;
      hp_t = m_hist_cft8_t;
    }
    hp_l->Reset();
    hp_t->Reset();
  }

  m_hist_piid_l->Reset();
  m_hist_piid_t->Reset();
  m_hist_bgo->Reset();

#endif

#if CATCH_ADC

  TH2 *hp_hi=0, *hp_lo=0;
  for (int layer=0; layer<8; layer++) {
    if (layer == 0) {
      hp_hi = m_hist_cft1_hi;
      hp_lo = m_hist_cft1_lo;
    } else if (layer == 1) {
      hp_hi = m_hist_cft2_hi;
      hp_lo = m_hist_cft2_lo;
    } else if (layer == 2) {
      hp_hi = m_hist_cft3_hi;
      hp_lo = m_hist_cft3_lo;
    } else if (layer == 3) {
      hp_hi = m_hist_cft4_hi;
      hp_lo = m_hist_cft4_lo;
    } else if (layer == 4) {
      hp_hi = m_hist_cft5_hi;
      hp_lo = m_hist_cft5_lo;
    } else if (layer == 5) {
      hp_hi = m_hist_cft6_hi;
      hp_lo = m_hist_cft6_lo;
    } else if (layer == 6) {
      hp_hi = m_hist_cft7_hi;
      hp_lo = m_hist_cft7_lo;
    } else if (layer == 7) {
      hp_hi = m_hist_cft8_hi;
      hp_lo = m_hist_cft8_lo;
    }
    hp_hi->Reset();
    hp_lo->Reset();
  }
#endif

}

//______________________________________________________________________________

void EventDisplay::ResetCATCH( void )
{
#if CATCH
  for (int layer=0; layer<NumOfPlaneCFT; layer++) {
    for (int seg=0; seg<NumOfSegCFT[layer]; seg++) {
      if (seg%32==0)
	m_CFT_Arc_cont[layer][seg]->SetLineColor(kOrange);
      else
	m_CFT_Arc_cont[layer][seg]->SetLineColor(kBlack);
    }
  }

  for (int seg=0; seg<NumOfSegBGO; seg++) {
    int size = m_BGO_Line_cont[seg].size();
    for (int i=0; i<size; i++)
      m_BGO_Line_cont[seg][i]->SetLineColor(kBlack);
  }

  for (int seg=0; seg<NumOfSegPiID; seg++) {
    int size = m_PiID_Line_cont[seg].size();
    for (int i=0; i<size; i++)
      m_PiID_Line_cont[seg][i]->SetLineColor(kBlack);
  }

  m_hbase_catch->Reset("ICES");

  del::DeleteObject( m_CFTTrack_cont);
  del::DeleteObject( m_CFTTrack_xy_cont);
  del::DeleteObject( m_CFTTrack_zx_cont);
  del::DeleteObject( m_CFTTrack_zy_cont);

  del::DeleteObject( m_SdcInTrack_Catch_cont);
  del::DeleteObject( m_SdcInTrack_Catch_xy_cont);
  del::DeleteObject( m_SdcInTrack_Catch_zx_cont);
  del::DeleteObject( m_SdcInTrack_Catch_zy_cont);

  //del::DeleteObject( m_CFTTrackCand_zx_cont);
  //del::DeleteObject( m_CFTTrackCand_zy_cont);
  del::DeleteObject( m_CATCH_dE_E_cont);
#endif

  del::DeleteObject( m_CFT_UV_cont);

}


//_____________________________________________________________________________
void
EventDisplay::FillMomentum(Double_t momentum)
{
#if Hist
  m_hist_p->Fill(momentum);
  m_canvas_hist->cd(1);
  gPad->Modified();
  gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillMassSquare(Double_t mass_square)
{
#if Hist
  m_hist_m2->Fill(mass_square);
  m_canvas_hist->cd(2);
  gPad->Modified();
  gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillMissMass(Double_t missmass)
{
#if Hist
  m_hist_missmass->Fill(missmass);
  m_canvas_hist->cd(3);
  gPad->Modified();
  gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillBH1(Int_t seg, Int_t tdc)
{
#if Hist_Timing
  Double_t p0 = gHodo.GetOffset(DetIdBH1, 0, seg, 0);
  Double_t p1 = gHodo.GetGain(DetIdBH1, 0, seg, 0);
  m_hist_bh1->Fill(seg, p1*((Double_t)tdc-p0));
#endif
}

//_____________________________________________________________________________
void
EventDisplay::SetCorrectTimeBH1(Int_t seg, Double_t de)
{
#if Hist_BcIn
  Int_t color;
  if (de<0.5)
    color = kBlue;
  else if (de >= 0.5 && de <= 1.5)
    color = kOrange;
  else
    color = kRed;
  m_BH1box_cont[seg]->SetFillColor(color);
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillBFT(Int_t layer, Int_t seg, Int_t tdc)
{
#if Hist_Timing
  Double_t p0 = gHodo.GetOffset(DetIdBFT, layer, seg, 0);
  Double_t p1 = gHodo.GetGain(DetIdBFT, layer, seg, 0);

  TH2 *hp=0;

  if (layer == 0)
    hp = m_hist_bft;
  else
    hp = m_hist_bft_p;

  hp->Fill(seg, p1*((Double_t)tdc-p0));
  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawBcInTrack(Double_t x0, Double_t u0)
{
#if Hist_Timing

  Double_t z1 = -150, z2 = 50;
  TLine *l = new TLine(x0+u0*(z1-zBFT), z1, x0+u0*(z2-zBFT), z2);
  m_BcInTrack.push_back(l);
#endif
}

//_____________________________________________________________________________
void
EventDisplay::SetCorrectTimeBFT(Double_t pos)
{
#if Hist_BcIn
  m_hist_bcIn->Fill(pos, zBFT);
#endif
}

//_____________________________________________________________________________
#if 0
void
EventDisplay::FillAFT(Int_t plane, Int_t seg, Double_t de_high)
{

  double posx = gAftHelper.GetX( plane, seg );
  double posz = gAftHelper.GetZ( plane, seg );
  if( plane%4 == 0 || plane%4 == 1 ) m_hist_aft_x->Fill(posz, posx, de_high);
  if( plane%4 == 2 || plane%4 == 3 ) m_hist_aft_y->Fill(posz, posx, de_high);
  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
}
#endif
//_____________________________________________________________________________
void
EventDisplay::FillBH2(Int_t seg, Int_t tdc)
{
#if Hist_Timing
  Double_t p0 = gHodo.GetOffset(DetIdBH2, 0, seg, 0);
  Double_t p1 = gHodo.GetGain(DetIdBH2, 0, seg, 0);
  m_hist_bh2->Fill(seg, p1*((Double_t)tdc-p0));
  //m_canvas_hist2->cd(1);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::SetCorrectTimeBH2(Int_t seg, Double_t de)
{
#if Hist_BcIn
  Int_t color;
  if (de<0.5)
    color = kBlue;
  else if (de >= 0.5 && de <= 1.5)
    color = kOrange;
  else
    color = kRed;
  m_BH2box_cont[seg]->SetFillColor(color);
  m_BH2box_cont2[seg]->SetFillColor(color);
#endif
}

//_____________________________________________________________________________
void
EventDisplay::SetCorrectTimeBcOut(Int_t layer, Double_t pos)
{
#if Hist_BcIn
  Double_t z = gGeom.GetLocalZ(layer);
  m_hist_bcOut->Fill(pos, z);

  m_hist_bcOut_sdcIn->Fill(pos + gGeom.GetGlobalPosition(layer).x(),
                           gGeom.GetGlobalPosition(layer).z());
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawBcOutTrack(Double_t x0, Double_t u0, Double_t y0, Double_t v0, Bool_t flagGoodForTracking )
{
#if Hist_BcIn

  Double_t z1 = 0, z2 = 600;
  TLine *l = new TLine(x0+u0*z1, z1, x0+u0*z2, z2);
  if (! flagGoodForTracking)
    l->SetLineColor(kOrange);

  m_BcOutTrack2.push_back(l);

  Double_t z3 = -50, z4 = 2950;
  TLine *l2 = new TLine(x0+u0*z3 + gxK18Target, z3 - (zK18Target-gzK18Target),
                        x0+u0*z4 + gxK18Target, z4 - (zK18Target-gzK18Target));
  if (! flagGoodForTracking)
    l2->SetLineColor(kOrange);

  m_BcOutTrack3.push_back(l2);

  TLine *l3 = new TLine(y0+v0*z3, z3 - (zK18Target-gzK18Target),
                        y0+v0*z4, z4 - (zK18Target-gzK18Target));
  if (! flagGoodForTracking)
    l3->SetLineColor(kOrange);

  m_BcOutTrack3.push_back(l3);

  for (Int_t layer=1; layer<=9; layer++) {
    Double_t z_bcOut = gGeom.GetLocalZ(layer) - gzK18Target + zK18Target;
    Double_t x = x0 + u0*z_bcOut + gxK18Target;
    Double_t y = y0 + v0*z_bcOut;
    Double_t z = gGeom.GetGlobalPosition(layer).z();

    ThreeVector gloPos(x, y, z);
    ThreeVector localPos = gGeom.Global2LocalPos(layer, gloPos);
    m_hist_sdcIn_predict->Fill(localPos.x(), z);

  }



#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawSdcInTrack(Double_t x0, Double_t u0, Double_t y0, Double_t v0, Bool_t flagHyps, Bool_t flagBeam)
{
#if Hist_BcIn

  Double_t z1 = -1500, z2 = 0;
  TLine *l = new TLine(x0+u0*z1, z1, x0+u0*z2, z2);
  if (flagHyps)
    l->SetLineColor(kRed);
  if (flagBeam)
    l->SetLineColor(kOrange);
  m_SdcInTrack2.push_back(l);

  TLine *l2 = new TLine(y0+v0*z1, z1, y0+v0*z2, z2);
  if (flagHyps)
    l2->SetLineColor(kRed);
  if (flagBeam)
    l2->SetLineColor(kOrange);

  m_SdcInTrack2.push_back(l2);


  for (Int_t layer=1; layer<=9; layer++) {
    Double_t z = gGeom.GetLocalZ(layer);
    Double_t x = x0 + u0*z;
    Double_t y = y0 + v0*z;
    Double_t gz = gGeom.GetGlobalPosition(layer).z();

    ThreeVector gloPos(x, y, gz);
    ThreeVector localPos = gGeom.Global2LocalPos(layer, gloPos);
    m_hist_sdcIn_predict2->Fill(localPos.x(), z);

  }

#endif
}

//_____________________________________________________________________________
void
EventDisplay::SetCorrectTimeSdcIn(Int_t layer, Double_t pos)
{
#if Hist_BcIn
  Double_t z = gGeom.GetLocalZ(layer);
  m_hist_bcOut_sdcIn->Fill(pos, z);

#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillTOF(Int_t seg, Int_t tdc)
{
#if Hist_Timing
  Double_t p0 = gHodo.GetOffset(DetIdTOF, 0, seg, 0);
  Double_t p1 = gHodo.GetGain(DetIdTOF, 0, seg, 0);

  m_hist_tof->Fill(seg, p1*((Double_t)tdc-p0));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillBcOutHit(Int_t layer,  Int_t wire, Int_t tdc)
{
#if Hist_Timing
  TH2 *hp=0;

  if (layer == 1)
    hp = m_hist_bc3;
  else if (layer == 2)
    hp = m_hist_bc3p;
  else if (layer == 3)
    hp = m_hist_bc3v;
  else if (layer == 4)
    hp = m_hist_bc3vp;
  else if (layer == 5)
    hp = m_hist_bc3u;
  else if (layer == 6)
    hp = m_hist_bc3up;
  else if (layer == 7)
    hp = m_hist_bc4u;
  else if (layer == 8)
    hp = m_hist_bc4up;
  else if (layer == 9)
    hp = m_hist_bc4v;
  else if (layer == 10)
    hp = m_hist_bc4vp;
  else if (layer == 11)
    hp = m_hist_bc4;
  else if (layer == 12)
    hp = m_hist_bc4p;

  Double_t p0=0.0, p1=-0.;
  gTdc.GetParameter(layer + 112, wire, p0, p1);

  hp->Fill(wire, -p1*(tdc+p0));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawBC3(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_bc3->Fill(wire, -1.*(tdc-351));
  m_hist_bc3_time->Fill(-1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawBC3p(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_bc3p->Fill(wire, -1.*(tdc-351));
  m_hist_bc3p_time->Fill(-1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//_____________________________________________________________________________
void
EventDisplay::DrawBC4(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_bc4->Fill(wire, -1.*(tdc-351));
  m_hist_bc4_time->Fill(-1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::DrawBC4p(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_bc4p->Fill(wire, -1.*(tdc-351));
  m_hist_bc4p_time->Fill(-1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillSDC1(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc1->Fill(wire, -1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//_____________________________________________________________________________
void
EventDisplay::FillSDC1p(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc1p->Fill(wire, -1.*(tdc-351));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//_____________________________________________________________________________
void
EventDisplay::FillSdcOutHit(Int_t layer,  Int_t wire, Int_t LorT, Int_t tdc)
{
#if Hist_SdcOut
  TH2 *hp=0;

  if (layer == 1 && LorT == 0)
    hp = m_hist_sdc3_l;
  else if (layer == 1 && LorT == 1)
    hp = m_hist_sdc3_t;
  else if (layer == 2 && LorT == 0)
    hp = m_hist_sdc3p_l;
  else if (layer == 2 && LorT == 1)
    hp = m_hist_sdc3p_t;
  else if (layer == 3 && LorT == 0)
    hp = m_hist_sdc3y_l;
  else if (layer == 3 && LorT == 1)
    hp = m_hist_sdc3y_t;
  else if (layer == 4 && LorT == 0)
    hp = m_hist_sdc3yp_l;
  else if (layer == 4 && LorT == 1)
    hp = m_hist_sdc3yp_t;
  else if (layer == 5 && LorT == 0)
    hp = m_hist_sdc4y_l;
  else if (layer == 5 && LorT == 1)
    hp = m_hist_sdc4y_t;
  else if (layer == 6 && LorT == 0)
    hp = m_hist_sdc4yp_l;
  else if (layer == 6 && LorT == 1)
    hp = m_hist_sdc4yp_t;
  else if (layer == 7 && LorT == 0)
    hp = m_hist_sdc4_l;
  else if (layer == 7 && LorT == 1)
    hp = m_hist_sdc4_t;
  else if (layer == 8 && LorT == 0)
    hp = m_hist_sdc4p_l;
  else if (layer == 8 && LorT == 1)
    hp = m_hist_sdc4p_t;

  Double_t p0=0.0, p1=-0.;
  gTdc.GetParameter(layer + 30, wire, p0, p1);

  hp->Fill(wire, -p1*(tdc+p0));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//_____________________________________________________________________________
void
EventDisplay::FillSDC3_Leading(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc3_l->Fill(wire, -0.833*(tdc-890));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//_____________________________________________________________________________
void
EventDisplay::FillSDC3_Trailing(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc3_t->Fill(wire, -0.833*(tdc-890));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillSDC3p_Leading(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc3p_l->Fill(wire, -0.833*(tdc-890));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillSDC3p_Trailing(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc3p_t->Fill(wire, -0.833*(tdc-890));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}


//_____________________________________________________________________________
void
EventDisplay::FillSDC4_Leading(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc4_l->Fill(wire, -0.833*(tdc-885));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillSDC4_Trailing(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc4_t->Fill(wire, -0.833*(tdc-885));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillSDC4p_Leading(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc4p_l->Fill(wire, -0.833*(tdc-885));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::FillSDC4p_Trailing(Int_t wire, Int_t tdc)
{
#if Hist_Timing
  m_hist_sdc4p_t->Fill(wire, -0.833*(tdc-885));
  //m_canvas_hist2->cd(3);
  //gPad->Modified();
  //gPad->Update();
#endif
}

//_____________________________________________________________________________
void
EventDisplay::Update()
{
  TIter canvas_iterator(gROOT->GetListOfCanvases());
  while (true) {
    auto canvas = dynamic_cast<TCanvas*>(canvas_iterator.Next());
    if (!canvas) break;
    canvas->UseCurrentStyle();
    canvas->cd(1)->SetLogz();
    canvas->Modified();
    canvas->Update();
  }

  UpdateCATCH();
// #if Hist_Timing
//   for (Int_t i=0; i<9; i++) {
//     m_canvas_hist2->cd(i+1);
//     gPad->Modified();
//     gPad->Update();
//   }

//   Int_t max1 = m_hist_bc3_time->GetMaximum();
//   Int_t max2 = m_hist_bc3p_time->GetMaximum();
//   if (max2 > max1)
//     m_hist_bc3_time->SetMaximum(max2*1.1);

//   max1 = m_hist_bc4_time->GetMaximum();
//   max2 = m_hist_bc4p_time->GetMaximum();
//   if (max2 > max1)
//     m_hist_bc4_time->SetMaximum(max2*1.1);

//   for (Int_t i=0; i<6; i++) {
//     m_canvas_hist3->cd(i+1);
//     gPad->Modified();
//     gPad->Update();
//   }
// #endif

// #if Hist_SdcOut
//   for (Int_t i=0; i<8; i++) {
//     m_canvas_hist4->cd(i+1);
//     gPad->Modified();
//     gPad->Update();
//   }
// #endif

// #if Hist_BcIn
//   for (Int_t i=0; i<4; i++) {
//     m_canvas_hist5->cd(i+1);


//     if (i==2) {
//       Int_t nc=m_BcInTrack.size();
//       for (Int_t i=0; i<nc; i++)
//  m_BcInTrack[i]->Draw("same");
//     } else if (i==3) {
//       Int_t nc=m_BcOutTrack2.size();
//       for (Int_t i=0; i<nc; i++)
//  m_BcOutTrack2[i]->Draw("same");
//     }


//     gPad->Modified();
//     gPad->Update();
//   }

//   {
//     m_canvas_hist6->cd();

//     Int_t nc=m_BcOutTrack3.size();
//     for (Int_t i=0; i<nc; i++)
//       m_BcOutTrack3[i]->Draw("same");

//     Int_t ncSdcIn=m_SdcInTrack2.size();
//     for (Int_t i=0; i<ncSdcIn; i++)
//  m_SdcInTrack2[i]->Draw("same");

//     gPad->Modified();
//     gPad->Update();
//   }

// #endif

}

//______________________________________________________________________________
void EventDisplay::UpdateCATCH( void )
{
#if CATCH
  //m_canvas_catch->cd(1);
  //m_canvas_catch->cd(1)->SetLogz(0);
  m_canvas->cd(2)->cd(1)->cd(1);
  m_canvas->cd(2)->cd(1)->cd(1)->SetLogz(0);

  m_hbase_catch->Draw("colz");

  for (int layer=0; layer<NumOfPlaneCFT; layer++) {
    for (int seg=0; seg<NumOfSegCFT[layer]; seg++) {
      m_CFT_Arc_cont[layer][seg]->Draw("same");
    }
  }

  for (int seg=0; seg<NumOfSegBGO; seg++) {
    int size = m_BGO_Line_cont[seg].size();
    for (int i=0; i<size; i++)
      m_BGO_Line_cont[seg][i]->Draw("same");
  }

  for (int seg=0; seg<NumOfSegPiID; seg++) {
    int size = m_PiID_Line_cont[seg].size();
    for (int i=0; i<size; i++)
      m_PiID_Line_cont[seg][i]->Draw("same");
  }


  for (int i=0; i<m_CFTTrack_xy_cont.size(); i++)
    m_CFTTrack_xy_cont[i]->Draw();

  for (int i=0; i<m_SdcInTrack_Catch_xy_cont.size(); i++)
    m_SdcInTrack_Catch_xy_cont[i]->Draw();

  /*
  for (int i=0; i<m_BcOutTrack_Catch_xy_cont.size(); i++)
    m_BcOutTrack_Catch_xy_cont[i]->Draw();

  for (int i=0; i<m_SigmaTrack1_xy_cont.size(); i++)
    m_SigmaTrack1_xy_cont[i]->Draw("pl");

  for (int i=0; i<m_SigmaTrack2_xy_cont.size(); i++)
    m_SigmaTrack2_xy_cont[i]->Draw("pl");
  */
  //m_canvas_catch->cd(2);
  m_canvas->cd(2)->cd(1)->cd(2);

  for (int i=0; i<m_CFTTrack_zx_cont.size(); i++)
    m_CFTTrack_zx_cont[i]->Draw();

  for (int i=0; i<m_SdcInTrack_Catch_zx_cont.size(); i++)
    m_SdcInTrack_Catch_zx_cont[i]->Draw();

  /*
  for (int i=0; i<m_BcOutTrack_Catch_zx_cont.size(); i++)
    m_BcOutTrack_Catch_zx_cont[i]->Draw();

  for (int i=0; i<m_SigmaTrack1_zx_cont.size(); i++)
    m_SigmaTrack1_zx_cont[i]->Draw("pl");

  for (int i=0; i<m_SigmaTrack2_zx_cont.size(); i++)
    m_SigmaTrack2_zx_cont[i]->Draw("pl");

  for (int i=0; i<m_CFTTrackCand_zx_cont.size(); i++)
    m_CFTTrackCand_zx_cont[i]->Draw("pl");
  */
  //m_canvas_catch->cd(3);
  m_canvas->cd(2)->cd(1)->cd(3);

  for (int i=0; i<m_CFTTrack_zy_cont.size(); i++)
    m_CFTTrack_zy_cont[i]->Draw();

  for (int i=0; i<m_SdcInTrack_Catch_zy_cont.size(); i++)
    m_SdcInTrack_Catch_zy_cont[i]->Draw();

  /*
  for (int i=0; i<m_BcOutTrack_Catch_zy_cont.size(); i++)
    m_BcOutTrack_Catch_zy_cont[i]->Draw();

  for (int i=0; i<m_SigmaTrack1_zy_cont.size(); i++)
    m_SigmaTrack1_zy_cont[i]->Draw("pl");

  for (int i=0; i<m_SigmaTrack2_zy_cont.size(); i++)
    m_SigmaTrack2_zy_cont[i]->Draw("pl");

  for (int i=0; i<m_CFTTrackCand_zy_cont.size(); i++)
    m_CFTTrackCand_zy_cont[i]->Draw("pl");
  */

  //m_canvas_catch->cd();
  //m_canvas_catch->Update();
  //m_canvas_catch->Modified();

  m_canvas->cd(2)->cd(1)->cd(1);
  m_canvas->cd(2)->cd(1)->Update();
  m_canvas->cd(2)->cd(1)->Modified();
#endif

#if CATCH3d

  m_canvas_catch3d->cd();
  m_geometry_catch->Draw();
  
  for (int i=0; i<m_CFT_UV_cont.size(); i++)
    m_CFT_UV_cont[i]->Draw();

  //for (int i=0; i<m_BcOutTrack_Catch_cont.size(); i++)
  //m_BcOutTrack_Catch_cont[i]->Draw();

  for (int i=0; i<m_SdcInTrack_Catch_cont.size(); i++)
    m_SdcInTrack_Catch_cont[i]->Draw();
  
  for (int i=0; i<m_CFTTrack_cont.size(); i++)
    m_CFTTrack_cont[i]->Draw();
  
  //for (int i=0; i<m_vertex3d_cont.size(); i++)
  //m_vertex3d_cont[i]->Draw();

  gPad->GetView()->ZoomIn();
  gPad->GetView()->ZoomIn();
  gPad->GetView()->ZoomIn();
  
  m_canvas_catch3d->Update();
  m_canvas_catch3d->Modified();
#endif

#if CATCH
  m_canvas_dE_E->cd();
  for (int i=0; i<m_CATCH_dE_E_cont.size(); i++)
    m_CATCH_dE_E_cont[i]->Draw("p");
  m_canvas_dE_E->Update();
  m_canvas_dE_E->Modified();
#endif
  
  /*
  m_canvas_scat->cd();
  m_canvas_scat->cd(1);
  for (int i=0; i<m_scat_dp1_cont.size(); i++)
    m_scat_dp1_cont[i]->Draw("p");
  m_canvas_scat->cd(2);
  for (int i=0; i<m_scat_dp2_cont.size(); i++)
    m_scat_dp2_cont[i]->Draw("p");
  m_canvas_scat->cd(3);
  for (int i=0; i<m_scat_dE1_cont.size(); i++)
    m_scat_dE1_cont[i]->Draw("p");
  m_canvas_scat->cd(4);
  for (int i=0; i<m_scat_dE2_cont.size(); i++)
    m_scat_dE2_cont[i]->Draw("p");

  m_canvas_scat->Update();
  m_canvas_scat->Modified();

  */



}


//_____________________________________________________________________________
void
EventDisplay::CalcRotMatrix(Double_t TA, Double_t RA1, Double_t RA2, Double_t *rotMat)
{
  Double_t ct0 = TMath::Cos(TA*TMath::DegToRad() );
  Double_t st0 = TMath::Sin(TA*TMath::DegToRad() );
  Double_t ct1 = TMath::Cos(RA1*TMath::DegToRad());
  Double_t st1 = TMath::Sin(RA1*TMath::DegToRad());
  Double_t ct2 = TMath::Cos(RA2*TMath::DegToRad());
  Double_t st2 = TMath::Sin(RA2*TMath::DegToRad());

  Double_t rotMat1[3][3], rotMat2[3][3];

  /* rotation matrix which is same as DCGeomRecord.cc */

#if 1
  /* new definition of RA2 */
  rotMat1[0][0] =  ct0*ct2+st0*st1*st2;
  rotMat1[0][1] = -st0*ct2+ct0*st1*st2;
  rotMat1[0][2] =  ct1*st2;

  rotMat1[1][0] =  st0*ct1;
  rotMat1[1][1] =  ct0*ct1;
  rotMat1[1][2] = -st1;

  rotMat1[2][0] = -ct0*st2+st0*st1*ct2;
  rotMat1[2][1] =  st0*st2+ct0*st1*ct2;
  rotMat1[2][2] =  ct1*ct2;

#else
  /* old definition of RA2 */
  rotMat1[0][0] =  ct0*ct2-st0*ct1*st2;
  rotMat1[0][1] = -st0*ct2-ct0*ct1*st2;
  rotMat1[0][2] =  st1*st2;

  rotMat1[1][0] =  ct0*st2+st0*ct1*ct2;
  rotMat1[1][1] = -st0*st2+ct0*ct1*ct2;
  rotMat1[1][2] = -st1*ct2;

  rotMat1[2][0] =  st0*st1;
  rotMat1[2][1] =  ct0*st1;
  rotMat1[2][2] =  ct1;
#endif

  /* rotation matrix which rotate -90 deg at x axis*/
  rotMat2[0][0] =  1.0;
  rotMat2[0][1] =  0.0;
  rotMat2[0][2] =  0.0;

  rotMat2[1][0] =  0.0;
  rotMat2[1][1] =  0.0;
  rotMat2[1][2] =  1.0;

  rotMat2[2][0] =  0.0;
  rotMat2[2][1] = -1.0;
  rotMat2[2][2] =  0.0;

  for (Int_t i=0; i<9; i++)
    rotMat[i]=0.0;

  for (Int_t i=0; i<3; i++) {
    for (Int_t j=0; j<3; j++) {
      for (Int_t k=0; k<3; k++) {
        //rotMat[3*i+j] += rotMat1[i][k]*rotMat2[k][j];
        rotMat[i+3*j] += rotMat1[i][k]*rotMat2[k][j];
      }
    }
  }

}

//_____________________________________________________________________________
Int_t
EventDisplay::GetCommand()
{
  //Update();
  char ch;
  char data[100];
  static Int_t stat   = 0;
  static Int_t Nevent = 0;
  static Int_t ev     = 0;

  const Int_t kSkip = 1;
  const Int_t kNormal = 0;

  if (stat == 1 && Nevent > 0 && ev<Nevent) {
    ev++;
    return kSkip;
  }
  if (ev==Nevent) {
    stat=0;
    ev=0;
  }

  if (stat == 0) {
    hddaq::cout << "q|n|p>" << std::endl;

    scanf("%c",&ch);
    if (ch!='\n')
      while(getchar() != '\n');

    switch (ch) {
    case 'q':
      std::exit(EXIT_SUCCESS);
    case 'n':
      stat = 1;
      do {
        printf("event#>");
        scanf("%s",data);
      } while ((Nevent=atoi(data))<=0);
      //hddaq::cout << "Continue " << Nevent << "event" << std::endl;
      hddaq::cout << "Skip " << Nevent << "event" << std::endl;
      break;
    case 'p':
      m_theApp->Run(kTRUE);
      break;
    }
  }

  if (stat == 1)
    return kSkip;

  return kNormal;
}

//_____________________________________________________________________________
void
EventDisplay::Print(Int_t run_number, Int_t event_number)
{
  static Int_t prev_run_number = 0;
  static TString fig_dir;
  if(run_number != prev_run_number){
    fig_dir = Form("fig/evdisp/run%05d", run_number);
    gSystem->MakeDirectory(fig_dir);
  }
  m_canvas->Print(Form("%s/evdisp_run%05d_ev%d.png",
                       fig_dir.Data(), run_number, event_number));
  prev_run_number = run_number;
}

//_____________________________________________________________________________
void
EventDisplay::Run(Bool_t flag)
{
  hddaq::cout << FUNC_NAME << " TApplication is running" << std::endl;

  m_theApp->Run(flag);
}


//_____________________________________________________________________________
void
EventDisplay::DrawHSTrack(Int_t nStep, const std::vector<TVector3>& StepPoint,
                              Double_t q)
{
  del::DeleteObject(m_hs_step_mark);

  m_hs_step_mark = new TPolyMarker3D(nStep);
  for(Int_t i=0; i<nStep; ++i){
    m_hs_step_mark->SetPoint(i,
			     StepPoint[i].x(),
			     StepPoint[i].y(),
			     StepPoint[i].z());
  }

  Color_t color = (q > 0) ? kRed : kBlue;

  if(!m_is_save_mode){
    m_hs_step_mark->SetMarkerSize(1);
    m_hs_step_mark->SetMarkerStyle(6);
  }
  m_hs_step_mark->SetMarkerColor(color);

  m_canvas->cd(1)->cd(2);
  m_hs_step_mark->Draw();
  m_canvas->Update();

}
