#include "Dataset.hpp"

#include "Config.hpp"
#include "tools/util.hpp"
#include "tools/io.hpp"

#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TTree.h>

/*
 * Arguments:
 * - `name`: simple name to use/display
 * - `files`: the files belonging to this dataset
 * - `xsec`: list of cross sections in pb ({0} for data files)
 */
Dataset::Dataset(std::string name,std::string label,std::string color,std::vector<std::string> files,std::vector<float> xsecs,float syst_unc,TString dataBasePath,bool isData,bool isSignal)
   : name(name)
   , label(label)
   , color(gROOT->ProcessLine((color+";").c_str()))
   , syst_unc(syst_unc)
   , isData(isData)
   , isSignal(isSignal)
{
   subsets.clear();
   for (uint i=0;i<files.size();i++){
      subsets.push_back(Datasubset(files[i],xsecs[i],dataBasePath,name,isData,isSignal));
   }
}

std::vector<TString> Dataset::getSubsetNames() const
{
   std::vector<TString> v;
   for (auto dss: subsets){
      v.push_back(dss.name);
   }
   return v;
}

/*
 * Arguments:
 * - `filename`: the filename of this datasubset
 * - `xsec`: cross section in pb (0 for data)
 */
Datasubset::Datasubset(std::string filename,float xsec,TString dataBasePath,std::string datasetName,bool isData,bool isSignal)
   : filename(filename)
   , xsec(xsec)
   , isData(isData)
   , isSignal(isSignal)
   , datasetName(datasetName)
{
   std::vector<std::string> splitString;
   boost::split(splitString,filename,boost::is_any_of("."));
   assert(splitString.size()>=2);
   name=splitString[0];

   TFile f(dataBasePath+filename);
   if (!f.IsZombie()){
      TH1F*  h=(TH1F*)f.Get("TreeWriter/hCutFlow");
      TTree* t=(TTree*)f.Get("TreeWriter/eventTree");
      if (!h || !t) {
         debug<<TString::Format("%s is broken!",filename.c_str());
         Ngen=0;
         entries=0;
      } else {
         Ngen=h->GetBinContent(2);
         entries=t->GetEntries();
      }
      f.Close();
   } else {
      debug<<TString::Format("%s is broken!",filename.c_str());
   }
}

/* return the full path to the file */
TString Datasubset::getPath() const{
   return Config::get().dataBasePath+filename;
}


DatasetCollection::DatasetCollection(boost::property_tree::ptree const& pt,TString dataBasePath)
{
   // MC
   std::vector<std::string> filenames;
   std::vector<float> xsecs;
   std::vector<float> kfacts;
   std::vector<float> feffs;
   std::string label;
   float syst_unc;
   for (std::string sDs: util::to_vector<std::string>(pt.get<std::string>("input.mc_datasets"))){
      filenames = util::to_vector<std::string>(pt.get<std::string>(sDs+".files"));
      xsecs = util::to_vector<float>(pt.get<std::string>(sDs+".xsecs"));
      assert(filenames.size()==xsecs.size());
      kfacts = util::to_vector<float>(pt.get<std::string>(sDs+".xsecs"));
      boost::optional<std::string> s_kfacts = pt.get_optional<std::string>(sDs+".kfact");
      if (s_kfacts) { // apply k factors
         kfacts=util::to_vector<float>(*s_kfacts);
         assert(kfacts.size()==xsecs.size());
         for (unsigned i=0; i<kfacts.size();i++) xsecs[i]*=kfacts[i];
      }
      boost::optional<std::string> s_feffs = pt.get_optional<std::string>(sDs+".filter_effs");
      if (s_feffs){ // apply filter efficiencies
         feffs = util::to_vector<float>(*s_feffs);
         assert(feffs.size()==xsecs.size());
         for (unsigned i=0; i<feffs.size();i++) xsecs[i]*=feffs[i];
      }
      label=pt.get<std::string>(sDs+".label", sDs); // use Dataset name if no explicit label given
      syst_unc=pt.get<float>(sDs+".syst_unc", 1e6); // default is something huge, so set it if you want to use it!
      mc_datasets_.push_back(Dataset(sDs,label,pt.get<std::string>(sDs+".color"),filenames,xsecs,syst_unc,dataBasePath));
   }
   for (std::string sDs: util::to_vector<std::string>(pt.get<std::string>("input.mc_alternative_datasets"))){
      filenames = util::to_vector<std::string>(pt.get<std::string>(sDs+".files"));
      xsecs = util::to_vector<float>(pt.get<std::string>(sDs+".xsecs"));
      assert(filenames.size()==xsecs.size());
      boost::optional<std::string> s_kfacts = pt.get_optional<std::string>(sDs+".kfact");
      if (s_kfacts) { // apply k factors
         kfacts=util::to_vector<float>(*s_kfacts);
         assert(kfacts.size()==xsecs.size());
         for (unsigned i=0; i<kfacts.size();i++) xsecs[i]*=kfacts[i];
      }
      boost::optional<std::string> s_feffs = pt.get_optional<std::string>(sDs+".filter_effs");
      if (s_feffs){ // apply filter efficiencies
         feffs = util::to_vector<float>(*s_feffs);
         assert(feffs.size()==xsecs.size());
         for (unsigned i=0; i<feffs.size();i++) xsecs[i]*=feffs[i];
      }
      label=pt.get<std::string>(sDs+".label", sDs); // use Dataset name if no explicit label given
      syst_unc=pt.get<float>(sDs+".syst_unc", 1e6); // default is something huge, so set it if you want to use it!
      mc_alternative_datasets_.push_back(Dataset(sDs,label,pt.get<std::string>(sDs+".color"),filenames,xsecs,syst_unc,dataBasePath));
   }
   // Signals
   for (std::string sDs: util::to_vector<std::string>(pt.get<std::string>("input.signals"))){
      filenames = util::to_vector<std::string>(pt.get<std::string>(sDs+".files"));
      xsecs = util::to_vector<float>(pt.get<std::string>(sDs+".xsecs"));
      assert(filenames.size()==xsecs.size());
      boost::optional<std::string> s_kfacts = pt.get_optional<std::string>(sDs+".kfact");
      if (s_kfacts) { // apply k factors
         kfacts=util::to_vector<float>(*s_kfacts);
         assert(kfacts.size()==xsecs.size());
         for (unsigned i=0; i<kfacts.size();i++) xsecs[i]*=kfacts[i];
      }
      boost::optional<std::string> s_feffs = pt.get_optional<std::string>(sDs+".filter_effs");
      if (s_feffs){ // apply filter efficiencies
         feffs = util::to_vector<float>(*s_feffs);
         assert(feffs.size()==xsecs.size());
         for (unsigned i=0; i<feffs.size();i++) xsecs[i]*=feffs[i];
      }
      label=pt.get<std::string>(sDs+".label", sDs); // use Dataset name if no explicit label given
      syst_unc=pt.get<float>(sDs+".syst_unc", 1e6); // default is something huge, so set it if you want to use it!
      signal_datasets_.push_back(Dataset(sDs,label,pt.get<std::string>(sDs+".color"),filenames,xsecs,syst_unc,dataBasePath,false,true));
   }
   // Data
   for (std::string sDs: util::to_vector<std::string>(pt.get<std::string>("input.data_streams"))){
      filenames = util::to_vector<std::string>(pt.get<std::string>(sDs+".files"));
      xsecs=std::vector<float>(filenames.size(),-1);
      assert(filenames.size()==xsecs.size());
      label=pt.get<std::string>(sDs+".label", sDs); // use Dataset name if no explicit label given
      data_datasets_.push_back(Dataset(sDs,label,pt.get<std::string>(sDs+".color"),filenames,xsecs,0,dataBasePath,true));
   }
   fillReferenceMaps();
}
DatasetCollection::DatasetCollection(std::vector<Dataset> mc_datasets,std::vector<Dataset> mc_alternative_datasets,std::vector<Dataset> data_datasets,std::vector<Dataset> signal_datasets)
   : mc_datasets_(mc_datasets)
   , mc_alternative_datasets_(mc_alternative_datasets)
   , data_datasets_(data_datasets)
   , signal_datasets_(signal_datasets)
{
   fillReferenceMaps();
}

void DatasetCollection::fillReferenceMaps()
{
   mDatasets_.clear();
   mDatasubsets_.clear();

   for (Dataset &ds: mc_datasets_){
      mDatasets_.emplace(ds.name,std::ref(ds));
      for (Datasubset &dss:ds.subsets){
         mDatasubsets_.emplace(dss.name,std::ref(dss));
      }
   }
   for (Dataset &ds: mc_alternative_datasets_){
      mDatasets_.emplace(ds.name,std::ref(ds));
      for (Datasubset &dss:ds.subsets){
         mDatasubsets_.emplace(dss.name,std::ref(dss));
      }
   }
   for (Dataset &ds: signal_datasets_){
      mDatasets_.emplace(ds.name,std::ref(ds));
      for (Datasubset &dss:ds.subsets){
         mDatasubsets_.emplace(dss.name,std::ref(dss));
      }
   }
   for (Dataset &ds: data_datasets_){
      mDatasets_.emplace(ds.name,std::ref(ds));
      for (Datasubset &dss:ds.subsets){
         mDatasubsets_.emplace(dss.name,std::ref(dss));
      }
   }
}

std::vector<Datasubset> DatasetCollection::getDatasubsets(bool mc, bool signal, bool data) const
{
   std::vector<Datasubset> vDss;
   if (mc)for (Dataset const &ds: mc_datasets_) // only defaults, no alternatives
      for (Datasubset const& dss: ds.subsets) vDss.push_back(dss);
   if (signal) for (Dataset const &ds: signal_datasets_)
      for (Datasubset const& dss: ds.subsets) vDss.push_back(dss);
   if (data)for (Dataset const &ds: data_datasets_)
      for (Datasubset const& dss: ds.subsets) vDss.push_back(dss);

   return vDss;
}

std::vector<Datasubset> DatasetCollection::getDatasubsets(std::vector<TString> dsNames) const
{
   std::vector<Datasubset> vDss;
   for (TString const &dsName: dsNames){
      for (Datasubset const &dss: getDataset(dsName).subsets){
         vDss.push_back(dss);
      }
   }
   return vDss;
}


Dataset DatasetCollection::getDataset(TString const &dsName) const
{
   return mDatasets_.at(dsName.Data());
}

Datasubset DatasetCollection::getDatasubset(TString const &dssName) const
{
   return mDatasubsets_.at(dssName.Data());
}

std::vector<std::string> DatasetCollection::getDatasetNames() const
{
   std::vector<std::string> v;
   for (auto ds: mDatasubsets_){
      v.push_back(ds.first);
   }
   return v;
}
std::vector<std::string> DatasetCollection::getDatasubsetNames() const
{
   std::vector<std::string> v;
   for (auto const &dss: mDatasubsets_){
      v.push_back(dss.first);
   }
   return v;
}
std::vector<std::string> DatasetCollection::getDatasubsetNames(std::vector<TString> dsNames) const
{
   std::vector<std::string> v;
   for (auto const &dss: getDatasubsets(dsNames)){
      v.push_back(dss.name);
   }
   return v;
}

TString DatasetCollection::getLabel(TString const &dsName) const
{
   Config const &cfg=Config::get();
   if (dsName=="efake") return cfg.efake.label;
   return getDataset(dsName).label;
}

float DatasetCollection::getSystUncert(TString dsName) const
{
   Config const &cfg=Config::get();
   if (dsName=="efake") return cfg.efake.syst_unc;
   if (dsName=="TTcomb") dsName="TTGJets";
   return getDataset(dsName).syst_unc;
}
