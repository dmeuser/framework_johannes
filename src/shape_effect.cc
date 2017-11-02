#include <string>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <map>

#include <TH1.h>
#include <TStyle.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TFile.h>
#include <TLatex.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>

using namespace std;

std::map<int,float> getVg1rel()
{
   std::string path = "/user/jschulz/2016/photonmet/output/SR_Vg_shape_variance.tex";
   unsigned const nCol= 9;
   std::ifstream fstream (path,std::ifstream::in);
   std::string line;
   std::map<int,float> mVg_value;
   while (fstream.good()) {
      std::getline(fstream, line);
      if (line.find("#") == 0) continue;
      if (line.find("%") == 0) continue;
      std::vector<float> values;
      values.push_back(line);
      if (values.size() == nCol) {
         mVg_value[(int)values[0]] = values[3];
      }
   }
   return mVg_value;
}

/*
std::map<int,float> getVg2rel()
{
   std::string path = "/user/jschulz/2016/photonmet/output/SR_Vg_shape_variance.tex";
   unsigned const nCol= 9;
   std::ifstream fstream (path,std::ifstream::in);
   std::string line;
   std::map<int,float> mVg_value;
   while (fstream.good()) {
      std::getline(fstream, line);
      if (line.find("#") == 0) continue;
      if (line.find("%") == 0) continue;
      line=util::rm_duplicate_spaces(line);
      std::vector<float> values=util::to_vector<float>(line,' ');
      if (values.size() == nCol) {
         mVg_value[(int)values[0]] = values[5];
      }
   }
   return mVg_value;
}

std::map<int,float> getVg3rel()
{
   std::string path = "/user/jschulz/2016/photonmet/output/SR_Vg_shape_variance.tex";
   unsigned const nCol= 9;
   std::ifstream fstream (path,std::ifstream::in);
   std::string line;
   std::map<int,float> mVg_value;
   while (fstream.good()) {
      std::getline(fstream, line);
      if (line.find("#") == 0) continue;
      if (line.find("%") == 0) continue;
      line=util::rm_duplicate_spaces(line);
      std::vector<float> values=util::to_vector<float>(line,' ');
      if (values.size() == nCol) {
         mVg_value[(int)values[0]] = values[7];
      }
   }
   return mVg_value;
}

std::map<int,float> getVg4rel()
{
   std::string path = "/user/jschulz/2016/photonmet/output/SR_Vg_shape_variance.tex";
   unsigned const nCol= 9;
   std::ifstream fstream (path,std::ifstream::in);
   std::string line;
   std::map<int,float> mVg_value;
   while (fstream.good()) {
      std::getline(fstream, line);
      if (line.find("#") == 0) continue;
      if (line.find("%") == 0) continue;
      line=util::rm_duplicate_spaces(line);
      std::vector<float> values=util::to_vector<float>(line,' ');
      if (values.size() == nCol) {
         mVg_value[(int)values[0]] = values[9];
      }
   }
   return mVg_value;
}*/

int shape_effect(){ 
	 	TCanvas *c1 = new TCanvas("c1","Plots",790,600);
   c1->cd();
   const double SR_bin[] = {600,800,1000,1300,1600};
   TH1F *h_Vg = new TH1F("h_Vg","h_Vg",4,SR_bin);
   h_Vg->SetBinContent(1,5);
   h_Vg->Draw("hist");

   return 0;
}


