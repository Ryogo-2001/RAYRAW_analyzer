#include "CFTPosParamMan.hh"

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>

#include <std_ostream.hh>

namespace
{
  const double offsetCATCH = 155;
}

CFTPosParamMan::CFTPosParamMan()
  : m_is_ready(false),
    m_file_name("")
{
  for (int i=0; i<NumOfPhiLayer; i++)
    for (int j=0; j<NumOfAngle; j++)
      for (int k=0; k<NumOfPhiPar; k++)
	PhiPosPar[i][j][k] = 0.;

  for (int i=0; i<NumOfULayer; i++)
    for (int j=0; j<NumOfSegCFT_UV4; j++)
      for (int k=0; k<NumOfUPar; k++)
	UPosPar[i][j][k] = 0.;

  for (int i=0; i<NumOfULayer; i++)
    for (int j=0; j<NumOfSegCFT_UV4; j++)
      for (int k=0; k<NumOfUDataPar; k++)
	UPosDataPar[i][j][k] = 0.;

  for (int i=0; i<NumOfULayer; i++)
    for (int j=0; j<NumOfSegCFT_UV4; j++)
      for (int k=0; k<NumOfTheta; k++)
	for (int l=0; l<NumOfUDataPar2; l++)
	  UPosDataPar2[i][j][k][l] = 0.;

}

CFTPosParamMan::~CFTPosParamMan()
{

}

bool CFTPosParamMan::Initialize( void )
{
  const std::string& class_name("CFTPosParamMan");

  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_ready ){
    std::cerr << "#W " << func_name
		<< " already initialied" << std::endl;
    return false;
  }


  std::ifstream ifs( m_file_name );
  if( !ifs.is_open() ){
    std::cerr << "#E " << func_name << " file open fail : "
		<< m_file_name << std::endl;
    return false;
  }

  int type_id;
  int plid, angle, seg, theta_id;
  double par[6];

  int invalid=0;
  std::string line;
  while( ifs.good() && std::getline(ifs, line) )
    {
      ++invalid;
      if(line[0] =='#' || line.empty())
	continue;

      std::istringstream input_line( line );
      input_line >> type_id ;

      if(type_id==0  || type_id==1){
	if( input_line >> plid >> angle
	    >> par[0] >> par[1] >> par[2] >> par[3] >> par[4] >> par[5]) {
	  if (type_id==0) {
	    for (int k=0; k<NumOfPhiPar; k++) {
	      PhiPosPar[plid][angle][k] = par[k];
	      //printf("%e ", par[k]);
	    }
	    //printf("\n");
	  } else if (type_id==1) {
	    for (int k=0; k<NumOfUPar; k++) {
	      UPosPar[plid][angle][k] = par[k];
	      //printf("%e ", par[k]);
	    }
	    //printf("\n");
	  }
	}
      } else if (type_id==2) {
	input_line >> plid >> seg;
	for (int k=0; k<NumOfUDataPar; k++) {
	  input_line >> UPosDataPar[plid][seg][k];
	}
      } else if (type_id==3) {
	input_line >> plid >> seg >> theta_id;
	for (int k=0; k<NumOfUDataPar2; k++) {
	  input_line >> UPosDataPar2[plid][seg][theta_id][k];
	}
      }

    }

  m_is_ready = true;
  return true;

}

//______________________________________________________________________________
bool
CFTPosParamMan::Initialize( const TString& file_name )
{
  m_file_name = file_name;
  return Initialize();
}


double  CFTPosParamMan::GetPhiShift(int layer, int angle, double z) const
{
  static const std::string funcname = "[CFTPosParamMan::GetPhiShift()]";
  if ( layer <0 || layer >= NumOfPhiLayer) {
    std::cerr << funcname << ": layer is invalid : " << layer << std::endl;
    return 0.;
  }

  if ( angle <0 || angle >= NumOfAngle) {
    std::cerr << funcname << ": angle is invalid : " << layer << std::endl;
    return 0.;
  }

  // z --> CFT local z pos
  z += offsetCATCH;


  if (z>= 0 && z <=400) {
    return PhiPosFunc(layer, angle, z);
  } else if (z < 0) {
    return PhiPosFunc(layer, angle, 0.);
  } else {
    return PhiPosFunc(layer, angle, 400.);
  }

  return 0.;
}

double  CFTPosParamMan::PhiPosFunc(int layer, int angle, double z) const
{
  return PhiPosPar[layer][angle][0]
    + PhiPosPar[layer][angle][1]*pow(z, 1.)
    + PhiPosPar[layer][angle][2]*pow(z, 2.)
    + PhiPosPar[layer][angle][3]*pow(z, 3.)
    + PhiPosPar[layer][angle][4]*pow(z, 4.)
    + PhiPosPar[layer][angle][5]*pow(z, 5.);
}

double  CFTPosParamMan::GetZposU(int layer, int seg, double phi, double theta) const
{
  double phi0 = UPosDataPar[layer][seg][0];

  double deltaShift2 = 0;
  if (phi-phi0 <= DivPhi/2)
    deltaShift2 = UPosDataPar[layer][seg][1];
  if (phi-phi0 >= 360 - DivPhi/2)
    deltaShift2 = UPosDataPar[layer][seg][NumOfUDataPar-1];
  else {
    int index0 = (int)(phi-phi0-DivPhi/2)/DivPhi+1;
    int index1 = index0+1;
    double y0 = UPosDataPar[layer][seg][index0];
    double delta_y = UPosDataPar[layer][seg][index1] - UPosDataPar[layer][seg][index0] ;
    double delta_x = phi-phi0 - ((index0-1)*(double)DivPhi + (double)DivPhi/2);

    deltaShift2 = y0 + delta_x/DivPhi*delta_y;
  }

  double deltaShift3 = 0;
  if (theta>=0 && theta<=180) {
    int thetaId = GetThetaId(theta);

    if (thetaId>=0 && thetaId<NumOfTheta) {
      if (phi-phi0 <= DivPhi2/2)
	deltaShift3 = UPosDataPar2[layer][seg][thetaId][1];
      if (phi-phi0 >= 360 - DivPhi2/2)
	deltaShift3 = UPosDataPar2[layer][seg][thetaId][NumOfUDataPar2-1];
      else {
	int index0 = (int)(phi-phi0-DivPhi2/2)/DivPhi2+1;
	int index1 = index0+1;
	double y0 = UPosDataPar2[layer][seg][thetaId][index0];
	double delta_y = UPosDataPar2[layer][seg][thetaId][index1] - UPosDataPar2[layer][seg][thetaId][index0] ;
	double delta_x = phi-phi0 - ((index0-1)*(double)DivPhi2 + (double)DivPhi2/2);

	deltaShift3 = y0 + delta_x/DivPhi2*delta_y;
	/*
	std::cout << "phi " << phi << ", index0 " << index0
		  << ", delta_y " << delta_y << ", delta_x " << delta_x
		  << ", deltaShift3 " << deltaShift3 << std::endl;
	*/
      }
    }
  }


  return UPosPar[layer][seg][0]
    + UPosPar[layer][seg][1]*pow(phi, 1.)
    + UPosPar[layer][seg][2]*pow(phi, 2.)
    + UPosPar[layer][seg][3]*pow(phi, 3.)
    + UPosPar[layer][seg][4]*pow(phi, 4.)
    + UPosPar[layer][seg][5]*pow(phi, 5.)
    + deltaShift2
    + deltaShift3;

}

int CFTPosParamMan::GetThetaId(double theta) const
{
  const std::string& class_name("CFTPosParamMan");
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  int index = -1;

  if (theta>=0 && theta<50)
    index = 0;
  else if (theta>=50 && theta<60)
    index = 1;
  else if (theta>=60 && theta<70)
    index = 2;
  else if (theta>=70 && theta<80)
    index = 3;
  else if (theta>=80 && theta<90)
    index = 4;
  else if (theta>=90 && theta<100)
    index = 5;
  else if (theta>=100 && theta<110)
    index = 6;
  else if (theta>=110 && theta<120)
    index = 7;
  else if (theta>=120 && theta<130)
    index = 8;
  else if (theta>=130 && theta<=180)
    index = 9;
  else {
    std::cerr << "#E " << func_name << " unexpected angle : "   << theta
	      << std::endl;
  }

  return index;

}
