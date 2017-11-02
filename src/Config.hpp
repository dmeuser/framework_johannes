#ifndef CONFIG_HPP__
#define CONFIG_HPP__

#include "Dataset.hpp"

#include <string>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <TString.h>

/* Singleton-like config class */
class Config
{
public:
   static Config& get();
   float lumi;
   float trigger_eff_PhMET;
   float trigger_eff_Ph;
   struct {
      float f,f_mc; // the fake rate
      float syst_unc;
      TString label;
      int color;
   } efake;
   struct {
      float Vg,e_Vg;
      float GJ,e_GJ;
      float rho;
      float uncert_Vgamma;     
      float uncert_gammaJ;
   } sf; // scale factors

   TString treeVersion;
   TString treeName;
   TString gitHash;
   TString dataBasePath;

   TString outputDirectory;

   float processFraction;
   bool  releaseMode;

   std::vector<std::string> modules;

   DatasetCollection datasets;

   // canvas decoration
   TString lumiText ;
   TString sqrtsText;
   TString extraText;
private:
   Config();
};

/* Singleton-like color class for automatic colors*/
class Color
{
public:
   static int next();
   static void reset();
private:
   std::vector<int> cols_={kGray,kGray+2,kAzure-8,kOrange-3,kCyan-3,kGreen-5,kRed-6,kGreen-6,kMagenta-6,kBlack};
   int size_=cols_.size();
   int pos_=-1;
   Color(){}
   static Color& instance_();
   int next_();
};
#endif /* CONFIG_HPP__ */
