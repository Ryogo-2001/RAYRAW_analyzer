#include "EventSelectMan.hh"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <std_ostream.hh>
#include "FuncName.hh"

//______________________________________________________________________________

EventSelectMan::EventSelectMan()
  : m_is_ready(false),
    m_file_name(""),
    m_current_point(0)
{

}


EventSelectMan::~EventSelectMan()
{

}


//______________________________________________________________________________
Bool_t
EventSelectMan::Initialize( void )
{
  if( m_is_ready ){
    //std::cerr << "#W " << func_name
    hddaq::cerr << FUNC_NAME
		<< " already initialied" << std::endl;
    return false;
  }


  std::ifstream ifs(m_file_name);
  if( !ifs.is_open() ){
    hddaq::cerr << FUNC_NAME << " file open fail : "
		<< m_file_name << std::endl;
    return false;
  }


  TString line;
  while( ifs.good() && line.ReadLine(ifs) ){
    if( line[0]=='#' || line.IsNull() ) continue;
    std::istringstream input_line( line.Data() );
    Int_t eventid;

    if( input_line >> eventid ){
      m_EventNumCont.push_back(eventid);
    }
  }

  std::stable_sort(m_EventNumCont.begin(), m_EventNumCont.end());

  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
Bool_t
EventSelectMan::Initialize( const TString& file_name )
{
  m_file_name = file_name;
  return Initialize();
}

//______________________________________________________________________________
Bool_t
EventSelectMan::RunNumCheck( Int_t runNum ) const
{
  TString run = Form("%d", runNum);

  if (m_file_name.Contains(run))
    return true;

  return false;

}

//______________________________________________________________________________
Bool_t
EventSelectMan::IsGood( Int_t eventNum ) const
{
  for (Int_t i=m_current_point; i<m_EventNumCont.size(); i++) {
    m_current_point = i;
    if (eventNum == m_EventNumCont[i]) 
      return true;
    else if (eventNum < m_EventNumCont[i])
      return false;
  }

  return false;

}
