#include "TemplateFitMan.hh"
#include "UnpackerManager.hh"
#include "Unpacker.hh"
#include "DeleteUtility.hh"

//namespace
//{
//  using namespace hddaq::unpacker::GUnpacker;
//}

TemplateFitFunction::TemplateFitFunction(TString filename)
  : m_area(0), m_sample_num(0), m_interval(0), m_center(0)
{
  std::string str;

  std::ifstream ifile(filename);

  Double_t xxx;
  Double_t yyy;

  if(!ifile)
    std::cout<<"Template File: "<< filename <<" is not exist"<<std::endl;

  Int_t nlines=0;

  m_area=0;

  Double_t max = 0;
  while(getline(ifile,str)){
    sscanf(str.data(),"%lf %lf",&xxx,&yyy);
    m_tempx.push_back(xxx);
    m_tempy.push_back(yyy);

    if(fabs(xxx) <0.0001)
      max = -yyy;

    m_area+=(Double_t) yyy;
    if(fabs(xxx)< 0.00001)
      m_center =nlines;

    nlines ++;

  }

  m_sample_num = nlines;

  for(int i=0;i<m_sample_num;i++)
    m_tempy[i] /= max;
  m_interval = m_tempx[1]-m_tempx[0];

}

TemplateFitFunction::~TemplateFitFunction()
{
  m_tempx.clear();
  m_tempy.clear();
  m_tempx.shrink_to_fit();
  m_tempy.shrink_to_fit();
}

Double_t TemplateFitFunction::myTemp(Double_t *x, Double_t *par)
{

  Double_t k=x[0];
  Double_t p=par[0];

  return par[1]*GetTemplateFunction(k-p);
}


Double_t TemplateFitFunction::GetTemplateFunction(Double_t x)
{

  Int_t xx = Int_t (x / m_interval) +m_center;
  if(x<m_tempx[0] || x >=m_tempx[m_sample_num-1])
    return 0;


  Int_t p;
  for(p= xx-5;p<xx+5;p++){
    if(x>= m_tempx[p] && x< m_tempx[p+1]){
      break;
    }
  }
  Double_t l = (x-m_tempx[p])/(m_tempx[p+1]-m_tempx[p]);

  return (m_tempy[p]+(m_tempy[p+1]-m_tempy[p])*l);

}


Double_t TemplateFitFunction::operator()(Double_t *x, Double_t *par)
{
  Double_t mix=0;
  for(int i=0; i<par[0] ; i++)
    mix += myTemp(&x[0],&par[2+2*i]);
  mix += par[1];

  return mix;
}

//______________________________________________________________________________

TemplateFitMan::TemplateFitMan()
  : m_is_ready(false),
    m_file_name(""),
    m_tempfunc_collection(),
    m_flag_ch14(false)
{
}


TemplateFitMan::~TemplateFitMan()
{
  for (int i=0; i<NumOfSegBGO; i++)
    if (m_fitFunction[i]) {
      delete m_fitFunction[i];
      m_fitFunction[i] = 0;
    }


  for (auto& elem: m_tempfunc_collection)
    del::ClearContainer(elem.second);

}


//______________________________________________________________________________
bool
TemplateFitMan::Initialize( void )
{
  const std::string& class_name("TemplateFitMan");

  static const TString func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_ready ){
    //hddaq::cerr << "#W " << func_name
    std::cerr << "#W " << func_name
		<< " already initialied" << std::endl;
    return false;
  }

  //std::vector <TString>filename_cont = split(m_file_name, '/');
  //int nCont = filename_cont.size();
  //TString filename = filename_cont[nCont-1];
  //TString runnum = filename.SubString(4, 4);
  //if (runnum == "7103")
  //m_flag_ch14 = true;

  for (int i=0; i<NumOfSegBGO; i++) {
    m_fitFunction[i] = new TemplateFitFunction(m_file_name+"."+std::to_string(i));
  }

  {
    TString name("BGO");
    auto& cont = m_tempfunc_collection[name];
    for (auto& hit: cont)
      delete hit;
    cont.clear();

    for (int i=0; i<NumOfSegBGO; i++) {
      TemplateFitFunction *tempfunc = new TemplateFitFunction(m_file_name+"."+std::to_string(i));
      cont.push_back(tempfunc);
    }
  }

  {
    TString name("TAG-PL");
    auto& cont = m_tempfunc_collection[name];
    for (auto& hit: cont)
      delete hit;
    cont.clear();

    for (int i=0; i<NumOfSegTagPL; i++) {
      // TemplateFitFunction *tempfunc = new TemplateFitFunction(m_file_name+"."+std::to_string(i));
      TemplateFitFunction *tempfunc = new TemplateFitFunction(m_file_name+"."+name);
      cont.push_back(tempfunc);
    }
  }

  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
bool
TemplateFitMan::Initialize( const TString& file_name )
{
  m_file_name = file_name;
  return Initialize();
}

//_____________________________________________________________________________
const TempFuncContainer&
TemplateFitMan::GetHitContainer(const TString& name) const
{
  auto itr = m_tempfunc_collection.find(name);
  if(itr == m_tempfunc_collection.end()){
    // throw Exception(FUNC_NAME + " No such detector: " + name);
    static TempFuncContainer null_container;
    return null_container;
  }else{
    return itr->second;
  }
}
