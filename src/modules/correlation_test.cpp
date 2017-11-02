#include "Config.hpp"

#include "tools/io.hpp"

#include <TMath.h>
#include <TRandom3.h>

TRandom3 g_rand(12387983);

// http://stats.stackexchange.com/questions/38856/how-to-generate-correlated-random-numbers-given-means-variances-and-degree-of/38867#38867
float correlatedValue(float x, float rho) {
   float r2 = rho*rho;
   float ve = 1-r2;
   float SD = TMath::Sqrt(ve);
   float e = g_rand.Gaus(0,SD);
   return rho*x + e;
}

extern "C"
void run()
{
   float m1=130;
   float s1=30;
   float m2=40;
   float s2=3;

   float rho=-0.6;

   TGraph gr;
   float x,y;
   for (int i=0; i<10000; i++) {
      x=g_rand.Gaus(0,1);
      y=correlatedValue(x,rho);
      x=m1+x*s1;
      y=m2+y*s2;
      gr.SetPoint(i,x,y);
   }

   io::RootFileSaver saver("plots.root","correlation_test",false);
   TCanvas can;
   gr.Draw("AP");
   gr.SetMarkerSize(.2);
   gr.Fit("pol1");
   io::log<<TString::Format("Mean   %.2f %.2f",gr.GetMean(1),gr.GetMean(2));
   io::log<<TString::Format("Std    %.2f %.2f",gr.GetRMS(1),gr.GetRMS(2));
   io::log<<TString::Format("rho    %.2f",gr.GetCorrelationFactor());
   saver.save(can,"scatter");
}
