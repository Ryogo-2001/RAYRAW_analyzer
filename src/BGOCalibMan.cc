#include "BGOCalibMan.hh"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <std_ostream.hh>

//______________________________________________________________________________

BGOCalibMan::BGOCalibMan()
  : m_is_ready(false),
    m_file_name("")
{
}


BGOCalibMan::~BGOCalibMan()
{

}


//______________________________________________________________________________
Bool_t
BGOCalibMan::Initialize( void )
{
  const TString& class_name("BGOCalibMan");

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

  for (Int_t i=0; i<NumOfSegBGO; i++) {
    m_Econt[i].clear();
    m_Hcont[i].clear();
  }

  Int_t invalid=0;
  std::string line;
  while( ifs.good() && std::getline( ifs, line ) ){
    ++invalid;
    if( line[0]=='#' || line.empty() ) continue;
    std::istringstream input_line( line );
    Int_t seg=-1, npoint=-1;
    Double_t E, PH;

    if( input_line >> seg >> npoint >> E >> PH ){
      m_Econt[seg].push_back(E);
      m_Hcont[seg].push_back(PH);
    }
  }

  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
Bool_t
BGOCalibMan::Initialize( const TString& file_name )
{
  m_file_name = file_name;
  return Initialize();
}

//______________________________________________________________________________
Bool_t
BGOCalibMan::GetEnergy( Int_t seg, Double_t pulse_height, Double_t &energy ) const
{
  if (!(seg>=0 && seg<NumOfSegBGO)) {
    return false;
  }

  Int_t icut = 0;

  Int_t Cut_Num = m_Hcont[seg].size();
  if(pulse_height<m_Hcont[seg][Cut_Num-1]){
    while((pulse_height>m_Hcont[seg][icut]) && (icut<Cut_Num-1)){
      icut+=1;
    }
  }else{
    icut = Cut_Num-1;
  }

  if(icut>0){
    energy  = m_Econt[seg][icut-1] + (m_Econt[seg][icut]-m_Econt[seg][icut-1])*(pulse_height-m_Hcont[seg][icut-1])/(m_Hcont[seg][icut]-m_Hcont[seg][icut-1]);
  }else{
    energy = 0;
  }

  return true;
}
