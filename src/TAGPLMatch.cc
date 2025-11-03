// -*- C++ -*-

#include "TAGPLMatch.hh"

#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>

#include <escape_sequence.hh>
#include <std_ostream.hh>
#include <UnpackerConfig.hh>
#include <UnpackerXMLReadDigit.hh>

#include "DetectorID.hh"
#include "FuncName.hh"
#include "Exception.hh"
#include "PrintHelper.hh"

//_____________________________________________________________________________
TAGPLMatch::Param::Param()
  : m_xmin(0.),
    m_xmax(0.)
{
}

//_____________________________________________________________________________
TAGPLMatch::Param::~Param()
{
}

//_____________________________________________________________________________
void
TAGPLMatch::Param::Print() const
{
  PrintHelper helper(1, std::ios::fixed);
  static const Int_t w = 5;
  hddaq::cout << " TAGPL seg " << std::right << std::setw(w) << m_seg << ":"
              << " (" << std::right
              << std::setw(w) << m_xmin << ", "
              << std::setw(w) << m_xmax << ")"
              << std::endl;
}

//_____________________________________________________________________________
TAGPLMatch::TAGPLMatch()
  : m_status(),
    m_param()
{
  m_status.reset();
}

//_____________________________________________________________________________
TAGPLMatch::~TAGPLMatch()
{
}

//_____________________________________________________________________________
Bool_t
TAGPLMatch::Initialize(const TString& file_name)
{
  std::ifstream ifs(file_name);
  if(!ifs.is_open()){
    std::cerr << FUNC_NAME << " No such parameter file("
              << file_name << ")" << std::endl;
    return false;
  }

  TString line;
  while(ifs.good() && line.ReadLine(ifs)){
    if(line.IsNull() || line[0] == '#') continue;
    std::istringstream iss(line.Data());
    std::istream_iterator<Double_t> itr_begin(iss);
    std::istream_iterator<Double_t> itr_end;
    std::vector<Double_t> cont(itr_begin, itr_end);
    if(cont.size() != kNParam){
      std::cerr << FUNC_NAME << std::endl
		<< " Number of parameters: " << cont.size() << std::endl
		<< " Required: " << kNParam << std::endl;
      throw Exception(FUNC_NAME + " invalid parameter.");
    }

    const Double_t tagplseg = cont[kTAGPLSegment];
    const Double_t xmin   = cont[kXMin];
    const Double_t xmax   = cont[kXMax];

    m_param[tagplseg].m_seg = tagplseg;
    m_param[tagplseg].m_xmin = xmin;
    m_param[tagplseg].m_xmax = xmax;
  }// read line

  if(m_status[kVerbose]) Print();
  m_status.set(kReady);

  return true;
}

//_____________________________________________________________________________
Bool_t
TAGPLMatch::Judge(Double_t bft_xpos, Double_t tagplseg)
{
  if(!m_status[kReady])
    throw Exception(FUNC_NAME + " is not initialized.");
  if(m_param.find(tagplseg) == m_param.end())
    return false;

  const Double_t xmin  = m_param.at(tagplseg).m_xmin;
  const Double_t xmax  = m_param.at(tagplseg).m_xmax;

  if(m_status[kVerbose]){
    hddaq::cout << FUNC_NAME << std::endl;
    m_param.at(tagplseg).Print();
    PrintHelper helper(2, std::ios::fixed);
    hddaq::cout << " SF pos " << std::setw(7) << bft_xpos << " -> ";
    if(xmin < bft_xpos && bft_xpos < xmax){
      hddaq::cout << hddaq::unpacker::esc::k_green
                  << "accept"
                  << hddaq::unpacker::esc::k_default_color << std::endl;
    } else {
      hddaq::cout << hddaq::unpacker::esc::k_red
                  << "reject"
                  << hddaq::unpacker::esc::k_default_color << std::endl;
    }
  }

  return (xmin < bft_xpos && bft_xpos < xmax);
}

//_____________________________________________________________________________
void
TAGPLMatch::Print() const
{
  hddaq::cout << FUNC_NAME << std::endl
              << "  seg: (xmin, xmax)" << std::endl;
  for(const auto& p : m_param) p.second.Print();
}
