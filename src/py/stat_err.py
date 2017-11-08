# example call: combine -M Asymptotic /user/lange/cmssw/CMSSW_7_6_5/src/HiggsAnalysis/CombinedLimit/data/tutorials/realistic-multi-channel.txt

import os, re, subprocess as sp, math
from distutils import spawn
from array import array
from operator import add,sub
import ROOT as rt

backgrounds=["GJ","TTcomb","Vg","diboson","efake"]

class Scan:
    T5gg,T5Wg,T6Wg,T6gg,GGM,TChiWg,TChiNg=range(7)

def prepareDatacard(obs,count,estat,esyst,eISR,correlated,mcUncertainties,pointName):
    dataCardPath="/tmp/datacard_%s.txt"%pointName
    nBins=len(obs)
    card="# "+pointName
    card+="""
imax %d  number of channels
jmax *  number of backgrounds ('*' = automatic)
kmax *  number of nuisance parameters (sources of systematical uncertainties)
------------
bin                        bin1       bin2       bin3       bin4
observation          """%nBins
    for o in obs:
        card+="%10d "%o

    card+="\n------------\n"
    allProc=["sig"]+backgrounds
    card+="bin                  "
    for i in range(nBins):
        card+=("      bin%d "%(i+1))*(len(backgrounds)+1)
    card+="\nprocess              "
    for i in range(nBins):
        for bkg in allProc:
            card+="%10s "%bkg
    card+="\nprocess              "
    for i in range(nBins):
        for j in range(len(backgrounds)+1):
            card+="%10d "%j
    card+="\nrate                 "
    for i in range(nBins):
        for b in allProc:
            card+="%10.2f "%count[b][i]

    card+="\n------------\n"
    # statistical errors
    col=0
    maxCol=nBins*len(allProc)
    fill="%10s "%"-"
    for i in range(nBins):
        for b in allProc:
            name="statB%d%s"%(i+1,b)
            card+="%-16s lnN "%name
            card+=fill*col
            card+="%10.2f "%estat[b][i]
            card+=fill*(maxCol-col-1)
            card+="\n"
            col+=1
    # systematic uncertainties
    for b,errors in esyst.iteritems():
        index=allProc.index(b)
        name="syst_%s"%b
        card+="%-16s lnN "%name
        # looping bins (all the same anyway...)
        for e in errors:
            card+=fill*index
            card+="%10.2f "%e
            card+=fill*(len(allProc)-index-1)
        card+="\n"
    for b,errors in eISR.iteritems():
        index=allProc.index(b)
        name="ISRsyst_%s"%b
        card+="%-16s lnN "%name
        # looping bins
        for e in errors:
            card+=fill*index
            card+="%10.2f "%e
            card+=fill*(len(allProc)-index-1)
        card+="\n"
    # correlated part of systematic uncertainties
    corr={}
    name="cor_"
    for b in correlated:
        corr[allProc.index(b)]=correlated[b]
        name+=b
    card+="%-16s lnN "%name
    for i in range(nBins):
        lastIndex=-1
        for ind in corr:
            card+=fill*(ind-lastIndex-1)
            card+="%10.2f "%corr[ind][i]
            lastIndex=ind
        card+=fill*(len(allProc)-lastIndex-1)
    card+="\n"
    # fully correlated uncertainties for all MC samples
    indices=[]
    for s in ["TTcomb","diboson","sig"]:
        indices.append(allProc.index(s))
    indices=sorted(indices)
    for name,val in mcUncertainties.iteritems():
        card+="%-16s lnN "%name
        for i in range(nBins):
            lastIndex=-1
            for ind in indices:
                card+=fill*(ind-lastIndex-1)
                card+="%10.2f "%val
                lastIndex=ind
            card+=fill*(len(allProc)-lastIndex-1)
        card+="\n"

    indices=[]
    for s in ["GJ"]:
        indices.append(allProc.index(s))
    indices=sorted(indices)
    card+="%-16s lnN "%"GJ_systematics"
    for i in range(nBins):
        lastIndex=-1
        for ind in indices:
            card+=fill*(ind-lastIndex-1)
            if i == 0:
               card+="%10.2f "%1.097
            elif i == 1:
               card+="%10.2f "%1.0797
            elif i == 2:
               card+="%10.2f "%1.0912
            elif i == 3:
               card+="%10.2f "%1.335
            else: print "error too many bins"
            lastIndex=ind
        card+=fill*(len(allProc)-lastIndex-1)
    card+="\n"
    indices=[]
    for s in ["Vg"]:
        indices.append(allProc.index(s))
    indices=sorted(indices)
    card+="%-16s lnN "%"Vg_systematics"
    for i in range(nBins):
        lastIndex=-1
        for ind in indices:
            card+=fill*(ind-lastIndex-1)
            if i == 0:
               card+="%10.2f "%1.0722
            elif i == 1:
               card+="%10.2f "%1.08045
            elif i == 2:
               card+="%10.2f "%1.08787
            elif i == 3:
               card+="%10.2f "%1.110
            else: print "error too many bins"
            lastIndex=ind
        card+=fill*(len(allProc)-lastIndex-1)
    card+="\n"
    with open(dataCardPath, "w") as f:
        f.write(card)
    return dataCardPath

def getSignalYield(point):
    f=rt.TFile(outdir+signal_scan,"read")
    hist=f.Get("pre_ph165/c_MET300/MT300/STg/"+point)
    hist_gen=f.Get("pre_ph165/c_MET300/MT300/STg/"+point+"_gen")
    hist_ISR=f.Get("pre_ph165/c_MET300/MT300/STg/"+point+"SRErrISR")  
    if math.isnan(hist.Integral()):
        # input tree for model was bad->Ngen=0->weight=inf
        f.Close()
        return None
    sigYield=[]
    statErr=[]
    metErr=[]
    ISRErr=[]
    for i in range(1,5):
        y=hist.GetBinContent(i)
        yg=hist_gen.GetBinContent(i)
        yISR=hist_ISR.GetBinContent(i)        
        sigYield.append(.5*(y+yg)) # use mean of reco and gen met yields
        if (y > 0):
           statErr.append(1+hist.GetBinError(i)/y)
           metErr.append(1+abs(y-yg)/(y+yg)) # uncertainty: half of the difference
           ISRErr.append(1+(abs(y-yISR)/y)) # uncertainty: full difference
        else:
           statErr.append(1+0.1)
           metErr.append(1+0.1) # uncertainty: half of the difference
           ISRErr.append(1+0.1) # uncertainty: full difference
    f.Close()
    return sigYield,statErr,metErr,ISRErr

def getMasses(point,scan):
    if scan==Scan.T5gg:   pattern="T5gg_(.*)_(.*)"
    elif scan==Scan.T5Wg:   pattern="T5Wg_(.*)_(.*)"
    elif scan==Scan.T6gg:   pattern="T6gg_(.*)_(.*)"   
    elif scan==Scan.T6Wg:   pattern="T6Wg_(.*)_(.*)"
    elif scan==Scan.GGM:    pattern=".*_(.*)_(.*)"
    elif scan==Scan.TChiWg: pattern="TChiWG_(.*)"
    elif scan==Scan.TChiNg: pattern="TChiNG_(.*)"
      
    m=re.search(pattern,point)
    masses=[]
    if m and scan==Scan.TChiWg and (len(m.groups())==1):
        masses=[int(m.groups()[0]),0]
    elif m and scan==Scan.TChiNg and (len(m.groups())==1):
        masses=[int(m.groups()[0]),0]        
    elif m and (len(m.groups())==2):
        for s in m.groups():
            masses.append(int(s))
    else:
        print "don't know what this is:",point
        exit(-1)
    return masses

def get_comb_xsec(mWino):
    # add C1C1 and N2C1 xs (in fb)
    xs=0.0
    with open(basedir+"xsec_C1C1_wino.csv") as f:
        for line in f:
            if line.startswith("#"): continue
            line=line.split()
            mW=int(line[0])
            if mW>=mWino:
                xs+=float(line[1])
                break
    with open(basedir+"xsec_N2C1_wino.csv") as f:
        for line in f:
            if line.startswith("#"): continue
            line=line.split()
            mW=int(line[0])
            if mW>=mWino:
                xs+=float(line[1])
                break
    return xs

def get_comb_xsec_unc(mWino):
    # add C1C1 and N2C1 uncertainties (in fb)
    xse=0.0
    with open(basedir+"xsec_C1C1_wino.csv") as f:
        for line in f:
            if line.startswith("#"): continue
            line=line.split()
            mW=int(line[0])
            if mW>=mWino:
                xse+=float(line[-1])
                break
    with open(basedir+"xsec_N2C1_wino.csv") as f:
        for line in f:
            if line.startswith("#"): continue
            line=line.split()
            mW=int(line[0])
            if mW>=mWino:
                xse+=float(line[-1])
                break
    return xse

def plot_uncert(scan):
    print "using",os.environ["CMSSW_VERSION"]
    print "combine tool",spawn.find_executable("combine")
    print

    scanFile=basedir
    if scan==Scan.T5gg:   scanFile+="T5gg_scan.txt"
    elif scan==Scan.T5Wg:   scanFile+="T5Wg_scan.txt"
    elif scan==Scan.T6gg:   scanFile+="T6gg_scan.txt"
    elif scan==Scan.T6Wg:   scanFile+="T6Wg_scan.txt"
    elif scan==Scan.GGM:    scanFile+="GGM_WinoBino_scan.txt"
    elif scan==Scan.TChiWg: scanFile+="TChiWg_scan.txt"
    elif scan==Scan.TChiNg: scanFile+="TChiNg_scan.txt"
    
    xsecFile=basedir
    if scan==Scan.T5gg:   xsecFile+="xsec_gluglu.csv"
    elif scan==Scan.T5Wg:   xsecFile+="xsec_gluglu.csv"
    elif scan==Scan.T6gg:   xsecFile+="xsec_sqsq.csv"
    elif scan==Scan.T6Wg:   xsecFile+="xsec_sqsq.csv"   
    elif scan==Scan.GGM:    xsecFile+="xsec_wino-bino.csv"
    elif scan==Scan.TChiWg: xsecFile+="xsec_N2C1_wino.csv"
    elif scan==Scan.TChiNg: xsecFile+="xsec_C1C1_wino.csv"
    
    sScan="unkown_scan"
    if scan==Scan.T5gg:   sScan="T5gg"
    elif scan==Scan.T5Wg:   sScan="T5Wg"
    elif scan==Scan.T6gg:   sScan="T6gg"
    elif scan==Scan.T6Wg:   sScan="T6Wg"
    elif scan==Scan.GGM:    sScan="GGM"
    elif scan==Scan.TChiWg: sScan="TChiWg"
    elif scan==Scan.TChiNg: sScan="TChiNg"
    
    count={}
    estat={}
    esyst={}
    eISR={}
    f=rt.TFile(outdir+"yields.root","read")
    for bkg in backgrounds:
        hist=f.Get("pre_ph165/c_MET300/MT300/STg/"+bkg)
        count[bkg]=[hist.GetBinContent(i) for i in range(1,5)]
        # subtract "additional" 1 from count index, because 1st bin is already left out
        estat[bkg]=[1+hist.GetBinError(i)/count[bkg][i-1] for i in range(1,5)]
        hist=f.Get("pre_ph165/c_MET300/MT300/STg/"+bkg+"_esyst")
        esyst[bkg]=[1+hist.GetBinContent(i)/count[bkg][i-1] for i in range(1,5)]
    hist=f.Get("pre_ph165/c_MET300/MT300/STg/data")
    obs=[int(hist.GetBinContent(i)) for i in range(1,5)]
    f.Close()

    # the same and 100% correlated for all MC
    mcUncertainties={
        "lumi"   : 1+ 2.6 /100,
        "phSF"   : 1+ 2. /100,
        "trigger": 1+ 0.43 /100,
    }

    points=[]
    with open(scanFile) as f:
        for p in f.read().split():
            p=p.split(".")[0]
            if scan==Scan.T5gg or scan==Scan.T5Wg or scan==Scan.T6Wg or scan==Scan.T6gg:
                m2,m1=getMasses(p,scan)
                if m2<1050: continue
            points.append(p)

    xsec={}
    xsec_err={}
    with open(xsecFile) as f:
        if scan==Scan.GGM:
            for line in f:
                if line.startswith("#"): continue
                line=line.split()
                key=int(line[0])*100000+int(line[1])
                xs=float(line[2])/1000 # to pb
                # xse=float(line[3]) # don't contain pdf
                xse=get_comb_xsec_unc(int(line[0]))/1000 # to pb
                xse/=xs
                xsec[key]=xs
                xsec_err[key]=xse
        elif scan==Scan.TChiWg:
            for line in f:
                if line.startswith("#"): continue
                xs=float(line.split()[1])/1000. # to pb
                xse=float(line.split()[2])/1000. # to pb
                xse/=xs
                xsec[int(line.split()[0])]=xs
                xsec_err[int(line.split()[0])]=xse
        elif scan==Scan.TChiNg:
            for line in f:
                if line.startswith("#"): continue
                xs=get_comb_xsec(int(line.split()[0]))/1000 # to pb
                xse=get_comb_xsec_unc(int(line.split()[0]))/1000 # to pb
                xse/=xs
                xsec[int(line.split()[0])]=xs
                xsec_err[int(line.split()[0])]=xse
        else:
            for line in f:
                if line.startswith("#"): continue
                xsec[int(line.split()[0])]=float(line.split()[1])
                xsec_err[int(line.split()[0])]=float(line.split()[2])/100.

    if scan==Scan.T5gg:   h_exp =rt.TH2F("","",35,800-25,2550+25,50,50,2550)
    if scan==Scan.T5Wg:   h_exp =rt.TH2F("","",35,800-25,2500+25,50,50,2550)
    if scan==Scan.T6gg:   h_exp =rt.TH2F("","",23,1000-25,2100+25,41,50,2100)    
    if scan==Scan.T6Wg:   h_exp =rt.TH2F("","",23,1000-25,2100+25,41,50,2100)
    if scan==Scan.GGM:    h_exp =rt.TH2F("","", 33, 205-12.5, 1005+12.5, 33, 215-12.5, 1015+12.5) # 32 points in each direction
    if scan==Scan.TChiWg: h_exp =rt.TH2F("","",41,300-12.5,1300+12.5,1,-1,1)
    if scan==Scan.TChiNg: h_exp =rt.TH2F("","",41,300-12.5,1300+12.5,1,-1,1)   
    h_exp_xs=rt.TH2F(h_exp)
    h_obs   =rt.TH2F(h_exp)
    h_obs_xs=rt.TH2F(h_exp)
    h_Sig=rt.TH2F(h_exp)
    h_StatB1=rt.TH2F(h_exp)
    h_StatB2=rt.TH2F(h_exp)
    h_StatB3=rt.TH2F(h_exp)
    h_StatB4=rt.TH2F(h_exp)
    
    g_StatB1=rt.TGraph2D()
    g_StatB2=rt.TGraph2D()
    g_StatB3=rt.TGraph2D()
    g_StatB4=rt.TGraph2D()

    g_MetB1=rt.TGraph2D()
    g_MetB2=rt.TGraph2D()
    g_MetB3=rt.TGraph2D()
    g_MetB4=rt.TGraph2D()

    g_ISRB1=rt.TGraph2D()
    g_ISRB2=rt.TGraph2D()
    g_ISRB3=rt.TGraph2D()
    g_ISRB4=rt.TGraph2D()
    
    print points            
    for i,point in enumerate(points):
        # print point,
        m2,m1=getMasses(point,scan)
        sigYield=getSignalYield(point)
        if not sigYield:
            print " broken!"
            continue # broken point
        c,st,sy,sISR=dict(count),dict(estat),dict(esyst),dict(eISR)
        # add signal
        c['sig']=sigYield[0]
        st['sig']=sigYield[1]
        sy['sig']=sigYield[2]
        sISR['sig']=sigYield[3]
       
        x,y=m2,m1
        if scan==Scan.GGM: x,y=m1,m2 # ggm display is transposed
        h_StatB1.SetBinContent(h_Sig.FindBin(x,y),100*(sigYield[1][0]-1))
        h_StatB2.SetBinContent(h_Sig.FindBin(x,y),100*(sigYield[1][1]-1))
        h_StatB3.SetBinContent(h_Sig.FindBin(x,y),100*(sigYield[1][2]-1))
        h_StatB4.SetBinContent(h_Sig.FindBin(x,y),100*(sigYield[1][3]-1))

        g_StatB1.SetPoint(i,x,y,100*(sigYield[1][0]-1))
        g_StatB2.SetPoint(i,x,y,100*(sigYield[1][1]-1))
        g_StatB3.SetPoint(i,x,y,100*(sigYield[1][2]-1))
        g_StatB4.SetPoint(i,x,y,100*(sigYield[1][3]-1))

        g_MetB1.SetPoint(i,x,y,100*(sigYield[2][0]-1))
        g_MetB2.SetPoint(i,x,y,100*(sigYield[2][1]-1))
        g_MetB3.SetPoint(i,x,y,100*(sigYield[2][2]-1))
        g_MetB4.SetPoint(i,x,y,100*(sigYield[2][3]-1))       

        if (scan == Scan.TChiWg or scan==Scan.TChiNg):
           g_ISRB1.SetPoint(i,x,y,100*(sigYield[3][0]-1))
           g_ISRB2.SetPoint(i,x,y,100*(sigYield[3][1]-1))
           g_ISRB3.SetPoint(i,x,y,100*(sigYield[3][2]-1))
           g_ISRB4.SetPoint(i,x,y,100*(sigYield[3][3]-1))

    f=rt.TFile(outdir+"uncert_%s.root"%sScan,"update")
    g_StatB1.Write("g_"+"statB1",rt.TObject.kOverwrite)
    g_StatB2.Write("g_"+"statB2",rt.TObject.kOverwrite)
    g_StatB3.Write("g_"+"statB3",rt.TObject.kOverwrite)
    g_StatB4.Write("g_"+"statB4",rt.TObject.kOverwrite)

    g_MetB1.Write("g_"+"MetB1",rt.TObject.kOverwrite)
    g_MetB2.Write("g_"+"MetB2",rt.TObject.kOverwrite)
    g_MetB3.Write("g_"+"MetB3",rt.TObject.kOverwrite)
    g_MetB4.Write("g_"+"MetB4",rt.TObject.kOverwrite)
    

    if (scan==Scan.TChiWg or scan==Scan.TChiNg):
       g_ISRB1.Write("g_"+"ISRB1",rt.TObject.kOverwrite)
       g_ISRB2.Write("g_"+"ISRB2",rt.TObject.kOverwrite)
       g_ISRB3.Write("g_"+"ISRB3",rt.TObject.kOverwrite)
       g_ISRB4.Write("g_"+"ISRB4",rt.TObject.kOverwrite)

    f.Close()

def redoHistogram(scan):
    sScan="unkown_scan"
    if scan==Scan.T5gg: sScan="T5gg"
    elif scan==Scan.T5Wg: sScan="T5Wg"
    elif scan==Scan.T6Wg: sScan="T6Wg"
    elif scan==Scan.T6gg: sScan="T6gg"
    elif scan==Scan.GGM: sScan="GGM"

    if (scan == Scan.T5gg):
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        label= "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma#tilde{G}"%(lsp_s,lsp_s)
        sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
    if (scan == Scan.T5Wg):
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        label= "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}"%(lsp_s,lsp_s)
        sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
    if (scan == Scan.T6gg):
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        label= "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma#tilde{G}"%(lsp_s,lsp_s)
        sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
    if (scan == Scan.T6Wg):
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        label= "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}"%(lsp_s,lsp_s)
        sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
    if (scan == Scan.GGM):
        lsp_s = ""
        label= "GGM"
        sParticle = "m_{#tilde{B}} (GeV)"
        LSP = "m_{#tilde{W}} (GeV)"
        

    f=rt.TFile(outdir+"limits_%s.root"%sScan,"update")
    gr=f.Get("gr_significance")
    h=gr.GetHistogram()
    h.GetZaxis().SetRangeUser(-3,3)
    h.GetXaxis().SetTitle(sParticle)
    h.GetYaxis().SetTitle(LSP) 
    for xbin in range(1,h.GetNbinsX()+1):
        for ybin in range(1,h.GetNbinsY()+1):
           if (h.GetBinContent(xbin,ybin) == 0):
              h.SetBinContent(xbin,ybin,-5)
    h.Write("h_significance_redone",rt.TObject.kOverwrite)
    f.Close()

def redoHistogramStat(scan):
    sScan="unkown_scan"
    if scan==Scan.T5gg: sScan="T5gg"
    elif scan==Scan.T5Wg: sScan="T5Wg"
    elif scan==Scan.T6Wg: sScan="T6Wg"
    elif scan==Scan.T6gg: sScan="T6gg"
    elif scan==Scan.GGM: sScan="GGM"

    if (scan == Scan.T5gg):
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        label= "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma#tilde{G}"%(lsp_s,lsp_s)
        sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
    if (scan == Scan.T5Wg):
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        label= "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}"%(lsp_s,lsp_s)
        sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
    if (scan == Scan.T6gg):
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        label= "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma#tilde{G}"%(lsp_s,lsp_s)
        sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
    if (scan == Scan.T6Wg):
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        label= "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}"%(lsp_s,lsp_s)
        sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
    if (scan == Scan.GGM):
        lsp_s = ""
        label= "GGM"
        sParticle = "m_{#tilde{B}} (GeV)"
        LSP = "m_{#tilde{W}} (GeV)"

    f=rt.TFile(outdir+"uncert_%s.root"%sScan,"update")
    gr=f.Get("g_statB1")
    h=gr.GetHistogram()
    h.SetTitle(sScan+" stat. uncertainty in bin1")
    rt.gStyle.SetPadRightMargin(0.15)
    rt.gStyle.SetPadBottomMargin(0.15)
    rt.gStyle.SetPadLeftMargin(0.15)
    h.GetZaxis().SetRangeUser(0,50)
    h.GetXaxis().SetTitleOffset(1.3)
    h.GetYaxis().SetTitleOffset(1.3)
    h.GetXaxis().SetTitleSize(1.5*h.GetXaxis().GetTitleSize())
    h.GetYaxis().SetTitleSize(1.5*h.GetYaxis().GetTitleSize())
    h.GetXaxis().SetTitle(sParticle)
    h.GetYaxis().SetTitle(LSP) 
    for xbin in range(1,h.GetNbinsX()+1):
        for ybin in range(1,h.GetNbinsY()+1):
           if (h.GetBinContent(xbin,ybin) == 0):
              h.SetBinContent(xbin,ybin,-5)
    c = rt.TCanvas()
    h.Draw("colz")
    rt.gPad.SaveAs(sScan+"_stat_B1.pdf")
    h.Write("g_StatB1_redone",rt.TObject.kOverwrite)
    gr=f.Get("g_statB2")
    h=gr.GetHistogram()
    h.GetZaxis().SetRangeUser(0,50)
    h.SetTitle(sScan+" stat. uncertainty in bin2")
    rt.gStyle.SetPadRightMargin(0.15)
    rt.gStyle.SetPadBottomMargin(0.15)
    rt.gStyle.SetPadLeftMargin(0.15)
    h.GetZaxis().SetRangeUser(0,50)
    h.GetXaxis().SetTitleOffset(1.3)
    h.GetYaxis().SetTitleOffset(1.3)
    h.GetXaxis().SetTitleSize(1.5*h.GetXaxis().GetTitleSize())
    h.GetYaxis().SetTitleSize(1.5*h.GetYaxis().GetTitleSize())
    h.GetXaxis().SetTitle(sParticle)
    h.GetYaxis().SetTitle(LSP)
    for xbin in range(1,h.GetNbinsX()+1):
        for ybin in range(1,h.GetNbinsY()+1):
           if (h.GetBinContent(xbin,ybin) == 0):
              h.SetBinContent(xbin,ybin,-5)
    c = rt.TCanvas()
    h.Draw("colz")
    rt.gPad.SaveAs(sScan+"_stat_B2.pdf")
    h.Write("g_StatB2_redone",rt.TObject.kOverwrite)
    gr=f.Get("g_statB3")
    h=gr.GetHistogram()
    h.GetZaxis().SetRangeUser(0,50)
    h.SetTitle(sScan+" stat. uncertainty in bin3")
    rt.gStyle.SetPadRightMargin(0.15)
    rt.gStyle.SetPadBottomMargin(0.15)
    rt.gStyle.SetPadLeftMargin(0.15)
    h.GetZaxis().SetRangeUser(0,50)
    h.GetXaxis().SetTitleOffset(1.3)
    h.GetYaxis().SetTitleOffset(1.3)
    h.GetXaxis().SetTitleSize(1.5*h.GetXaxis().GetTitleSize())
    h.GetYaxis().SetTitleSize(1.5*h.GetYaxis().GetTitleSize())
    h.GetXaxis().SetTitle(sParticle)
    h.GetYaxis().SetTitle(LSP)
    for xbin in range(1,h.GetNbinsX()+1):
        for ybin in range(1,h.GetNbinsY()+1):
           if (h.GetBinContent(xbin,ybin) == 0):
              h.SetBinContent(xbin,ybin,-5)
    c = rt.TCanvas()
    h.Draw("colz")
    rt.gPad.SaveAs(sScan+"_stat_B3.pdf")
    h.Write("g_StatB3_redone",rt.TObject.kOverwrite)
    gr=f.Get("g_statB4")
    h=gr.GetHistogram()
    h.GetZaxis().SetRangeUser(0,50)
    h.SetTitle(sScan+" stat. uncertainty in bin4")
    rt.gStyle.SetPadRightMargin(0.15)
    rt.gStyle.SetPadBottomMargin(0.15)
    rt.gStyle.SetPadLeftMargin(0.15)
    h.GetZaxis().SetRangeUser(0,50)
    h.GetXaxis().SetTitleOffset(1.3)
    h.GetYaxis().SetTitleOffset(1.3)
    h.GetXaxis().SetTitleSize(1.5*h.GetXaxis().GetTitleSize())
    h.GetYaxis().SetTitleSize(1.5*h.GetYaxis().GetTitleSize())
    h.GetXaxis().SetTitle(sParticle)
    h.GetYaxis().SetTitle(LSP)
    for xbin in range(1,h.GetNbinsX()+1):
        for ybin in range(1,h.GetNbinsY()+1):
           if (h.GetBinContent(xbin,ybin) == 0):
              h.SetBinContent(xbin,ybin,-5)
    c = rt.TCanvas()
    h.Draw("colz")
    rt.gPad.SaveAs(sScan+"_stat_B4.pdf")
    h.Write("g_StatB4_redone",rt.TObject.kOverwrite)
    f.Close()

def redoHistogramMet(scan):
    sScan="unkown_scan"
    if scan==Scan.T5gg: sScan="T5gg"
    elif scan==Scan.T5Wg: sScan="T5Wg"
    elif scan==Scan.T6Wg: sScan="T6Wg"
    elif scan==Scan.T6gg: sScan="T6gg"
    elif scan==Scan.GGM: sScan="GGM"

    if (scan == Scan.T5gg):
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        label= "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma#tilde{G}"%(lsp_s,lsp_s)
        sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
    if (scan == Scan.T5Wg):
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        label= "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}"%(lsp_s,lsp_s)
        sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
    if (scan == Scan.T6gg):
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        label= "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma#tilde{G}"%(lsp_s,lsp_s)
        sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
    if (scan == Scan.T6Wg):
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        label= "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}"%(lsp_s,lsp_s)
        sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
    if (scan == Scan.GGM):
        lsp_s = ""
        label= "GGM"
        sParticle = "m_{#tilde{B}} (GeV)"
        LSP = "m_{#tilde{W}} (GeV)"

    f=rt.TFile(outdir+"uncert_%s.root"%sScan,"update")
    gr=f.Get("g_MetB1")
    h=gr.GetHistogram()
    h.GetZaxis().SetRangeUser(0,20)
    h.SetTitle(sScan+" FastSim p_{T}^{miss} uncertainty in bin2")
    rt.gStyle.SetPadRightMargin(0.15)
    rt.gStyle.SetPadBottomMargin(0.15)
    rt.gStyle.SetPadLeftMargin(0.15)
    h.GetZaxis().SetRangeUser(0,50)
    h.GetXaxis().SetTitleOffset(1.3)
    h.GetYaxis().SetTitleOffset(1.3)
    h.GetXaxis().SetTitleSize(1.5*h.GetXaxis().GetTitleSize())
    h.GetYaxis().SetTitleSize(1.5*h.GetYaxis().GetTitleSize())
    h.GetXaxis().SetTitle(sParticle)
    h.GetYaxis().SetTitle(LSP)
    for xbin in range(1,h.GetNbinsX()+1):
        for ybin in range(1,h.GetNbinsY()+1):
           if (h.GetBinContent(xbin,ybin) == 0):
              h.SetBinContent(xbin,ybin,-5)
    c = rt.TCanvas()
    h.Draw("colz")
    rt.gPad.SaveAs(sScan+"_met_B1.pdf")
    h.Write("g_MetB1_redone",rt.TObject.kOverwrite)
    gr=f.Get("g_MetB2")
    h=gr.GetHistogram()
    h.GetZaxis().SetRangeUser(0,20)
    h.SetTitle(sScan+" FastSim p_{T}^{miss} uncertainty in bin1")
    rt.gStyle.SetPadRightMargin(0.15)
    rt.gStyle.SetPadBottomMargin(0.15)
    rt.gStyle.SetPadLeftMargin(0.15)
    h.GetZaxis().SetRangeUser(0,50)
    h.GetXaxis().SetTitleOffset(1.3)
    h.GetYaxis().SetTitleOffset(1.3)
    h.GetXaxis().SetTitleSize(1.5*h.GetXaxis().GetTitleSize())
    h.GetYaxis().SetTitleSize(1.5*h.GetYaxis().GetTitleSize())
    h.GetXaxis().SetTitle(sParticle)
    h.GetYaxis().SetTitle(LSP)
    for xbin in range(1,h.GetNbinsX()+1):
        for ybin in range(1,h.GetNbinsY()+1):
           if (h.GetBinContent(xbin,ybin) == 0):
              h.SetBinContent(xbin,ybin,-5)
    c = rt.TCanvas()
    h.Draw("colz")
    rt.gPad.SaveAs(sScan+"_met_B2.pdf")
    h.Write("g_MetB2_redone",rt.TObject.kOverwrite)
    gr=f.Get("g_MetB3")
    h=gr.GetHistogram()
    h.GetZaxis().SetRangeUser(0,20)
    h.SetTitle(sScan+" FastSim p_{T}^{miss} uncertainty in bin3")
    rt.gStyle.SetPadRightMargin(0.15)
    rt.gStyle.SetPadBottomMargin(0.15)
    rt.gStyle.SetPadLeftMargin(0.15)
    h.GetZaxis().SetRangeUser(0,50)
    h.GetXaxis().SetTitleOffset(1.3)
    h.GetYaxis().SetTitleOffset(1.3)
    h.GetXaxis().SetTitleSize(1.5*h.GetXaxis().GetTitleSize())
    h.GetYaxis().SetTitleSize(1.5*h.GetYaxis().GetTitleSize())
    h.GetXaxis().SetTitle(sParticle)
    h.GetYaxis().SetTitle(LSP)
    for xbin in range(1,h.GetNbinsX()+1):
        for ybin in range(1,h.GetNbinsY()+1):
           if (h.GetBinContent(xbin,ybin) == 0):
              h.SetBinContent(xbin,ybin,-5)
    c = rt.TCanvas()
    h.Draw("colz")
    rt.gPad.SaveAs(sScan+"_met_B3.pdf")
    h.Write("g_MetB3_redone",rt.TObject.kOverwrite)
    gr=f.Get("g_MetB4")
    h=gr.GetHistogram()
    h.GetZaxis().SetRangeUser(0,20)
    h.SetTitle(sScan+" FastSim p_{T}^{miss} uncertainty in bin4")
    rt.gStyle.SetPadRightMargin(0.15)
    rt.gStyle.SetPadBottomMargin(0.15)
    rt.gStyle.SetPadLeftMargin(0.15)
    h.GetZaxis().SetRangeUser(0,50)
    h.GetXaxis().SetTitleOffset(1.3)
    h.GetYaxis().SetTitleOffset(1.3)
    h.GetXaxis().SetTitleSize(1.5*h.GetXaxis().GetTitleSize())
    h.GetYaxis().SetTitleSize(1.5*h.GetYaxis().GetTitleSize())
    h.GetXaxis().SetTitle(sParticle)
    h.GetYaxis().SetTitle(LSP)
    for xbin in range(1,h.GetNbinsX()+1):
        for ybin in range(1,h.GetNbinsY()+1):
           if (h.GetBinContent(xbin,ybin) == 0):
              h.SetBinContent(xbin,ybin,-5)
    c = rt.TCanvas()
    h.Draw("colz")
    rt.gPad.SaveAs(sScan+"_met_B4.pdf")
    h.Write("g_MetB4_redone",rt.TObject.kOverwrite)
    f.Close()

if __name__ == '__main__':
    basedir="/user/jschulz/2016/photonmet/"
    outdir=basedir+"output/"
    signal_scan="signal_scan_v19.root"
    rho=-0.0
    
#    plot_uncert(Scan.TChiWg)
    
    plot_uncert(Scan.TChiNg)
    
 #   plot_uncert(Scan.T5gg)
 #   redoHistogram(Scan.T5gg)
   # redoHistogramStat(Scan.T5gg)
 #   redoHistogramMet(Scan.T5gg)

 #   plot_uncert(Scan.T5Wg)
 #   redoHistogram(Scan.T5Wg)
  #  redoHistogramStat(Scan.T5Wg)
 #   redoHistogramMet(Scan.T5Wg)

  #  plot_uncert(Scan.T6gg)
 #   redoHistogram(Scan.T6gg)
  #  redoHistogramStat(Scan.T6gg)
  #  redoHistogramMet(Scan.T6gg)

  #  plot_uncert(Scan.T6Wg)
 #   redoHistogram(Scan.T6Wg)
 #   redoHistogramStat(Scan.T6Wg)
  #  redoHistogramMet(Scan.T6Wg)

 #   plot_uncert(Scan.GGM)
 #   redoHistogram(Scan.GGM)
 #   redoHistogramStat(Scan.GGM)
 #   redoHistogramMet(Scan.GGM)
#


