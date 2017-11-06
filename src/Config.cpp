#include "Config.hpp"

#include "tools/io.hpp"
#include "tools/util.hpp"

#include <TROOT.h>

//static
Config& Config::get(){
   static Config instance;
   return instance;
}

Config::Config()
{
   boost::property_tree::ptree pt;
   std::string cfgFile(CMAKE_SOURCE_DIR);
   cfgFile+="config.ini";
   boost::property_tree::read_ini(cfgFile,pt);

   treeVersion=pt.get<std::string>("input.version");
   treeName=pt.get<std::string>("input.treeName");
   dataBasePath=pt.get<std::string>("input.dataBasePath")+treeVersion+"/";
   gitHash=io::shellOutput("git log -1 --pretty=format:%h");
   if (TString(io::shellOutput("git status")).Contains("modified")) gitHash+="*";
   lumi=pt.get<float>("general.lumi");

   trigger_eff_Ph   =pt.get<float>("general.trigger_eff_Ph")   /100.0;
   trigger_eff_PhMET=pt.get<float>("general.trigger_eff_PhMET")/100.0;

   std::vector<float> vsf=util::to_vector<float>(pt.get<std::string>("sf.Vg"));
   assert(vsf.size()==2);
   sf.Vg=vsf[0];
   sf.e_Vg=vsf[1];
   vsf=util::to_vector<float>(pt.get<std::string>("sf.GJ"));
   assert(vsf.size()==2);
   sf.GJ=vsf[0];
   sf.e_GJ=vsf[1];
   sf.rho=pt.get<float>("sf.rho");
   sf.uncert_Vgamma=pt.get<float>("sf.uncert_Vgamma");   
   sf.uncert_gammaJ=pt.get<float>("sf.uncert_gammaJ");

   efake.f=pt.get<float>("efake.f")/100.0;
   efake.f_mc=pt.get<float>("efake.f_mc")/100.0;
   efake.syst_unc=pt.get<float>("efake.syst_unc");
   efake.label=pt.get<std::string>("efake.label");
   efake.color=gROOT->ProcessLine((pt.get<std::string>("efake.color")+";").c_str());

   lumiText=TString::Format("%.1f fb^{-1}",lumi*1e-3);
   sqrtsText=pt.get<std::string>("general.sqrtsText");
   extraText=pt.get<std::string>("general.extraText");

   outputDirectory=pt.get<std::string>("output.directory");
   datasets=DatasetCollection(pt,dataBasePath);
}

// static
Color& Color::instance_(){
   static Color instance;
   return instance;
}
// static
int Color::next(){
   return instance_().next_();
}
// static
void Color::reset(){
   instance_().pos_=-1;
}
int Color::next_(){
   pos_=(pos_+1)%size_;
   return cols_[pos_];
}

