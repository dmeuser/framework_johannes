{
   // macro has to be loaded before:
   setTDRStyle();

   // own styles:

   // Calculating the same t+b and l+r margins does
   // not give quadratic result (height does not fit).
   // These settings are determined experimentally.
   gStyle->SetPadTopMargin   (.06);
   gStyle->SetPadBottomMargin(.13);
   gStyle->SetPadLeftMargin  (.16);
   gStyle->SetPadRightMargin (.06); // fits with 4 digits on y-axis

   // calling with kNAME does not work at PZ
   // kBird=57 does not look right at PZ...
   // kDeepSea=51, kRainBow=55,kBlueYellow= 54
 //  gStyle->SetPalette(54);
//   gStyle->SetPalette(112);
   gStyle->SetPalette(55);
   gStyle->SetHistLineWidth(2);
   gStyle->SetHistLineColor(kGray+2);
   gStyle->SetHistFillColor(kGray+2);
   gStyle->SetHistFillStyle(0);
   gStyle->SetNumberContours(255);

   gStyle->SetGridStyle(1);
   gStyle->SetGridColor(kGray+1);

   gStyle->SetNdivisions(505,"xy");

   TGaxis::SetExponentOffset(-.15,+.005,"y");
   TGaxis::SetExponentOffset(-.04,-.11 ,"x");
   TGaxis::SetMaxDigits(4);

   gStyle->SetTitleOffset(.92,"x");
   gStyle->SetTitleOffset(1.30,"y");
   gStyle->SetLabelOffset(-.006,"z");
   gStyle->SetTitleOffset(1.35,"z");
   gStyle->SetTitleSize(0.05, "z");

   gStyle->SetLegendFont(42);

   gStyle->SetHatchesSpacing(1);
   gStyle->SetHatchesLineWidth(2);

}
