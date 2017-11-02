#ifndef DATASET_HPP__
#define DATASET_HPP__

#include <string>
#include <vector>
#include <iostream>
#include <functional>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/algorithm/string.hpp>

#include <TString.h>

class Datasubset;

/* class representing the dataset of a physics process or data */
class Dataset
{
public:
   Dataset(std::string name,std::string label,std::string color,std::vector<std::string> files,std::vector<float> xsecs,float syst_unc,TString dataBasePath,bool isData=false,bool isSignal=false);
   std::vector<TString> getSubsetNames() const;
   std::string name;
   TString label;
   int color;
   int entries;
   float syst_unc;
   bool isData;
   bool isSignal;
   std::vector<Datasubset> subsets;
};

/* class representing a subset of a dataset as  processed, e.g. a HT bin of the whole process */
class Datasubset
{
public:
   Datasubset(std::string filename,float xsec,TString dataBasePath,std::string datasetName,bool isData,bool isSignal);
   TString getPath() const;
   std::string filename,name;
   float xsec;
   int entries,Ngen;
   bool isData;
   bool isSignal;
   std::string const datasetName;
};

/* Container to store datasets */
class DatasetCollection
{
public:
   DatasetCollection(){}
   DatasetCollection(boost::property_tree::ptree const& pt,TString dataBasePath);
   DatasetCollection(std::vector<Dataset> mc_datasets,std::vector<Dataset> mc_alternative_datasets,std::vector<Dataset> data_datasets,std::vector<Dataset> signal_datasets);
   std::vector<Dataset>    getDatasets() const {return mc_datasets_;}
   std::vector<Datasubset> getDatasubsets(bool mc=true, bool signal=false, bool data=false) const;
   std::vector<Datasubset> getDatasubsets(std::vector<TString> dsNames) const;
   Dataset    getDataset(TString const &dsName) const;
   Datasubset getDatasubset(TString const &dssName) const;
   std::vector<std::string> getDatasetNames() const;
   std::vector<std::string> getDatasubsetNames() const;
   std::vector<std::string> getDatasubsetNames(std::vector<TString> dsNames) const;
   TString getLabel(TString const &dsName) const;
   float getSystUncert(TString dsName) const;
private:
   std::vector<Dataset> mc_datasets_;
   std::vector<Dataset> mc_alternative_datasets_;
   std::vector<Dataset> data_datasets_;
   std::vector<Dataset> signal_datasets_;
   std::map<std::string,std::reference_wrapper<Dataset>> mDatasets_;
   std::map<std::string,std::reference_wrapper<Datasubset>> mDatasubsets_;

   void fillReferenceMaps();
};

#endif /* DATASET_HPP__ */
