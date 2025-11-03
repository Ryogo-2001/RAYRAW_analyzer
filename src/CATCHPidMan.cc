#include "CATCHPidMan.hh"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <std_ostream.hh>

//______________________________________________________________________________

CATCHPidMan::CATCHPidMan()
  : m_is_ready(false),
    m_file_name("")
{
}


CATCHPidMan::~CATCHPidMan()
{

}


//______________________________________________________________________________
Bool_t
CATCHPidMan::Initialize( void )
{
  const TString& class_name("CATCHPidMan");

  static const TString func_name("["+class_name+"::"+__func__+"()]");

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


  Int_t invalid=0;
  std::string line;
  while( ifs.good() && std::getline( ifs, line ) ){
    ++invalid;
    if( line[0]=='#' || line.empty() ) continue;
    std::istringstream input_line( line );
    Int_t parid=-1;
    Double_t p0, p1, p2, p3;

    if( input_line >> parid >> p0 >> p1 >> p2 >> p3 ){
      if (parid==0) {
	m_ParLow[0] = p0;
	m_ParLow[1] = p1;
	m_ParLow[2] = p2;
	m_ParLow[3] = p3;
      } else if (parid==1) {
	m_ParHigh[0] = p0;
	m_ParHigh[1] = p1;
	m_ParHigh[2] = p2;
	m_ParHigh[3] = p3;
      }
    }
  }

  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
Bool_t
CATCHPidMan::Initialize( const TString& file_name )
{
  m_file_name = file_name;
  return Initialize();
}

//______________________________________________________________________________
Bool_t
CATCHPidMan::CheckProton( Double_t BGO_E, Double_t CFT_dE, Double_t &delta ) const
{
  // lower
  Double_t cut1 = m_ParLow[0]*exp(BGO_E*m_ParLow[1]) + m_ParLow[2]*exp(BGO_E*m_ParLow[3]);

  // higher
  Double_t cut2 = m_ParHigh[0]*exp(BGO_E*m_ParHigh[1]) + m_ParHigh[2]*exp(BGO_E*m_ParHigh[3]);

  delta = CFT_dE - (cut1+cut2)/2.;

  if (CFT_dE >= cut1 && CFT_dE <= cut2)
    return true;
  else
    return false;


}

//______________________________________________________________________________
Bool_t
CATCHPidMan::CheckPi( Double_t BGO_E, Double_t CFT_dE) const
{
  // lower
  Double_t cut1 = m_ParLow[0]*exp(BGO_E*m_ParLow[1]) + m_ParLow[2]*exp(BGO_E*m_ParLow[3]);

  if (CFT_dE < cut1 && CFT_dE > 0)
    return true;
  else
    return false;


}

