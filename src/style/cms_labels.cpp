#include "cms_labels.hpp"

#include "tools/io.hpp"

#include "Config.hpp"

#include <iostream>

#include "TLatex.h"

void style::draw_lumi(TPad &pad, bool simulation, bool drawLumiText)
{
   TString cmsText     = "CMS";
   TString extraText=Config::get().extraText;
   TString lumiText =Config::get().lumiText;
   TString sqrtsText=Config::get().sqrtsText;

   if (simulation){
      if (extraText!=""){
         extraText="simulation #bullet "+extraText;
      }else{
         extraText="simulation";
      }
   }

   float cmsTextFont   = 61;  // default is helvetic-bold
   float extraTextFont = 52;  // default is helvetica-italics
   // text sizes and text offsets with respect to the top frame
   // in unit of the top margin size
   float lumiTextSize     = 0.6;
   float lumiTextOffset   = 0.2;
   float cmsTextSize      = 0.75;
   // ratio of "CMS" and extra text size
   float extraOverCmsTextSize  = 0.76;

   ////////////////////////////////////

   // float H = pad.GetWh();
   // float W = pad.GetWw();
   float l = pad.GetLeftMargin();
   float t = pad.GetTopMargin();
   float r = pad.GetRightMargin();
   // float b = pad.GetBottomMargin();

   pad.cd();

   TLatex latex;
   latex.SetNDC();
   latex.SetTextAngle(0);
   latex.SetTextColor(kBlack);

   float extraTextSize = extraOverCmsTextSize*cmsTextSize;

   latex.SetTextFont(42);
   latex.SetTextAlign(31);
   latex.SetTextSize(lumiTextSize*t);
   if (drawLumiText) latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText+" "+sqrtsText);

   latex.SetTextAlign(11);
   latex.SetTextSize(cmsTextSize*t);
   TString cmsLabel;
   cmsLabel+=TString::Format("#font[%f]{%s}",cmsTextFont,cmsText.Data());
   cmsLabel+=TString::Format("#color[922]{#scale[%f]{#font[%f]{ %s}}}",
                             extraTextSize/float(cmsTextSize),extraTextFont,extraText.Data());

   latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsLabel);

   pad.Update();
}
