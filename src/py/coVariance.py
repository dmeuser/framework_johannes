import os, re, subprocess as sp, math
from distutils import spawn
from array import array
from operator import add,sub
import ROOT as rt

def covarinaceMatrix(datacard):
    newDatacard = "Newdatacard_TChiWG_1000.txt"
    with open(datacard) as f:
        lines = f.readlines()
        lines[4] = "shapes * * FAKE"
    with open(newDatacard, "w+") as f:
        f.write("\n".join(lines))
    ###subprocess.call(["text2workspace.py", "--X-allow-no-signal", "--X-allow-no-background", datacard])
    ###subprocess.call(["combine", "-M", "MaxLikelihoodFit", "--saveShapes", "--saveWithUnc", "--numToysForShape", "2000", "--saveOverall", "--preFitValue", "0", datacard.replace(".txt", ".root")])
#    sp.call(["text2workspace.py", "--channel-masks", "--X-allow-no-signal", "--X-allow-no-background", newDatacard])
#    sp.call(["combine", "-M", "MaxLikelihoodFit", "--saveShapes", "--saveWithUnc", "--numToysForShape", "2000", "--setPhysicsModelParameterRange", "mask_Signal1=1,mask_Signal2=1", "--saveOverall", newDatacard.replace(".txt", ".root")])
 #   sp.call(["text2workspace.py", "--X-allow-no-signal","--X-allow-no-background", newDatacard])
 #   sp.call(["combine", "-M", "MaxLikelihoodFit", "--saveShapes","--saveWithUnc", "--numToysForShape", "2000", "--saveOverallShapes",newDatacard.replace(".txt",".root"),  "--preFitValue", "0"])

    drawCovarinaceMatrix("mlfit.root")

def drawCovarinaceMatrix(filename):
    f = rt.TFile(filename)
    h2 = f.Get("shapes_prefit/overall_total_covar")    
    h2.SetTitle("    #bf{CMS} #it{Preliminary}                              35.9 fb^{-1} (13 TeV)")
    h2.SetTitleSize(100)
    h2.GetZaxis().SetTitle("Covariance between SR bins   ")
    h2.GetXaxis().SetBinLabel(1, "600-800")
    h2.GetYaxis().SetBinLabel(1, "600-800")
    h2.GetXaxis().SetBinLabel(2, "800-1000")
    h2.GetYaxis().SetBinLabel(2, "800-1000")
    h2.GetXaxis().SetBinLabel(3, "1000-1300")
    h2.GetYaxis().SetBinLabel(3, "1000-1300")
    h2.GetXaxis().SetBinLabel(4, "1300-#infty  ")
    h2.GetYaxis().SetBinLabel(4, "1300-#infty  ")
    h2.GetXaxis().SetTitle("S_{T}^{#gamma} (GeV)")
    h2.GetYaxis().SetTitle("S_{T}^{#gamma} (GeV)")
    h2.GetYaxis().LabelsOption("u")
    h2.GetXaxis().SetTitleSize(1.3*h2.GetZaxis().GetTitleSize())
    h2.GetYaxis().SetTitleSize(1.3*h2.GetZaxis().GetTitleSize())
    h2.GetZaxis().SetTitleSize(1.4*h2.GetZaxis().GetTitleSize())  
    h2.GetXaxis().SetLabelSize(1.7*h2.GetXaxis().GetLabelSize())
    h2.GetYaxis().SetLabelSize(1.2*h2.GetXaxis().GetLabelSize())
    h2.GetZaxis().SetLabelSize(0.7*h2.GetXaxis().GetLabelSize())   
    h2.GetXaxis().SetTitleOffset(1.2)
    h2.GetYaxis().SetTitleOffset(2.2)
    h2.GetZaxis().SetTitleOffset(1.1)
    h2.GetZaxis().SetRangeUser(0,1)
    h2.SetMarkerSize(1.8)
  #  h2.GetZaxis().SetRangeUser(1,750)
    h2.GetZaxis().SetRangeUser(0,1)
    rt.gStyle.SetOptStat(0000)
    rt.gStyle.SetPaintTextFormat("3.2f");
    rt.gStyle.SetPadLeftMargin(0.2)
    rt.gStyle.SetPadRightMargin(0.17)
    rt.gStyle.SetPadBottomMargin(0.15)
    rt.gStyle.SetPadTopMargin(0.07)     
    rt.gStyle.SetPalette(55)
    rt.gStyle.SetNdivisions(505,"xyz");
    c = rt.TCanvas()
    c.SetLogz(1)
    h2.Draw("colz text")
    rt.gPad.SaveAs("covariance.pdf")
   
    totalBkg = f.Get("shapes_prefit/total_background")
    bkgUnc = [totalBkg.GetBinError(bin) for bin in range(1,totalBkg.GetNbinsX()+1)]
    print bkgUnc
    for xbin in range(1,h2.GetNbinsX()+1):
        for ybin in range(1,h2.GetNbinsY()+1):
           xe = bkgUnc[xbin-1]
           ye = bkgUnc[ybin-1]
           if xe and ye:
               h2.SetBinContent(xbin,ybin,h2.GetBinContent(xbin,ybin)/xe/ye)
    c = rt.TCanvas()
  #  h2.GetZaxis().SetRangeUser(5,250)
    h2.GetZaxis().SetTitle("Correlation between SR bins   ")
    h2.Draw("colz text")
    c.RedrawAxis()
   # h2.GetZaxis().SetRangeUser(5,250)
    rt.gPad.SaveAs("correlation.pdf")

def main():
    covarinaceMatrix("datacard_TChiWG_1000.txt")

if __name__ == "__main__":
    main()
