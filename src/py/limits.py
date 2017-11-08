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

def parseCombineToolOutput(datacard):
    combineOutput=sp.check_output(["combine","-M","Asymptotic","--rMax","5","--rMin", "-5",datacard])
    def getLim(line):
        return float(line.split("<")[1])
    # print combineOutput
    limit={}
    for line in combineOutput.split("\n"):
        if "Observed Li" in line: limit["obs"] = getLim(line)
        if "Expected 50" in line: limit["exp"] = getLim(line)
        if "Expected 84" in line: limit["exp+1"] = getLim(line)
        if "Expected 16" in line: limit["exp-1"] = getLim(line)
        if "Expected 97" in line: limit["exp+2"] = getLim(line)
        if "Expected  2" in line: limit["exp-2"] = getLim(line)
    return limit

def parseCombineToolSignificanceOutput(datacard):
    combineOutput = sp.check_output(["combine", "-M", "ProfileLikelihood", "--significance", "--uncapped", "1", "--rMin", "-5", datacard], stderr=sp.STDOUT)
    significance=-5
    for line in combineOutput.split("\n"):
        if "Significance: " in line: significance = float(line.split(":")[1])
    return significance

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

def getSignalContamination(point):
    f=rt.TFile(outdir+signal_scan,"read")
    hist=f.Get("pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh/"+point)
    c=hist.Integral()
    if math.isnan(c):
        # input tree for model was bad->Ngen=0->weight=inf
        f.Close()
        return None
    f.Close()
    return c

def decomposeCorrelations(esyst,count):
    cor={"Vg":[],"GJ":[]}
    unc={"GJ":[]}
    for i in range(len(count["Vg"])):
        eX=(esyst["Vg"][i]-1)*count["Vg"][i]
        eY=(esyst["GJ"][i]-1)*count["GJ"][i]
        eXtoY=rho*eY
        eYunc=(1-rho**2)**.5 * eY
        # print eY,eXtoY,eYunc
        cor["Vg"].append((count["Vg"][i]+eX)/count["Vg"][i])
        cor["GJ"].append((count["GJ"][i]+eXtoY)/count["GJ"][i])
        unc["GJ"].append((count["GJ"][i]+eYunc)/count["GJ"][i])
    return cor,unc

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

def getContour(gr2):
    """ inspired by https://hypernews.cern.ch/HyperNews/CMS/get/susy/2122/2.html
    """
    c=rt.TCanvas()
    gr2.Draw("tri1")
    c.Update()
    cont=gr2.GetContourList(1.)
    if cont and len(cont)>0:
        N=0
        c=None
        for ci in cont:
            if ci.GetN() > N:
                N=ci.GetN()
                c=ci
        return c
    else:
        print "GetContourList failed"
        return None

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

def caclulateLimits(scan):
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

    gr={}
    #grSig={}
    grSig=rt.TGraph2D()
    for lvl in ["obs","obs+1","obs-1","exp","exp+1","exp-1","exp+2","exp-2","obs_xs"]:
        gr[lvl]=rt.TGraph2D()
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
    print points            
    for i,point in enumerate(points):
        print point,
        m2,m1=getMasses(point,scan)
        key=m2
        if scan==Scan.GGM: key=m2*100000+m1
        xs=xsec[key]
        sigYield=getSignalYield(point)
        contamin=getSignalContamination(point)
        if not sigYield:
            print " broken!"
            continue # broken point
        c,st,sy,sISR=dict(count),dict(estat),dict(esyst),dict(eISR)
        # decompose partially correlated parts:
        cor,unc=decomposeCorrelations(esyst,count)
        del sy["GJ"]
        del sy["Vg"]
        sy.update(unc) # re-add the uncorrelated part
        # add signal
        c['sig']=sigYield[0]
        st['sig']=sigYield[1]
        sy['sig']=sigYield[2]
        sISR['sig']=sigYield[3]
        # subtract bkg overestimation from signal contamination
        subtractGJ=[x*contamin for x in c["GJ"]]
        subtractVG=[x*contamin for x in c["Vg"]]
        subtract=map(add, subtractGJ, subtractVG)
        # actual subtraction done from S to not destroy the B-only hypothesis
        c["sig"]=map(sub,c["sig"],subtract)
        # avoid negative counts
        c["sig"]=[max(0,x) for x in c["sig"]]
        datacard=prepareDatacard(obs,c,st,sy,sISR,cor,mcUncertainties,point)
      #  continue
        rLimits=parseCombineToolOutput(datacard)
        Significance=parseCombineToolSignificanceOutput(datacard)
        e_xs_rel=xsec_err[key]
        rLimits["obs+1"]=rLimits["obs"]*(1+e_xs_rel)
        rLimits["obs-1"]=rLimits["obs"]*(1-e_xs_rel)
        x,y=m2,m1
        if scan==Scan.GGM: x,y=m1,m2 # ggm display is transposed
        for lvl in ["obs","obs+1","obs-1","exp","exp+1","exp-1","exp+2","exp-2"]:
            rLim=rLimits[lvl]
            gr[lvl].SetPoint(i,x,y,rLim)
        if scan==Scan.T5gg: h_Sig.SetBinContent(h_Sig.FindBin(x,y),Significance)
        if scan==Scan.T5Wg: h_Sig.SetBinContent(h_Sig.FindBin(x,y),Significance)
        if scan==Scan.T6gg: h_Sig.SetBinContent(h_Sig.FindBin(x,y),Significance)
        if scan==Scan.T6Wg: h_Sig.SetBinContent(h_Sig.FindBin(x,y),Significance)
        if scan==Scan.GGM: h_Sig.SetBinContent(h_Sig.FindBin(x,y),Significance)       
        grSig.SetPoint(i,x,y,Significance)
        rLim=rLimits['exp']
        h_exp.SetBinContent(h_exp.FindBin(x,y),rLim)
        h_exp_xs.SetBinContent(h_exp.FindBin(x,y),rLim*xs)
        print "exp. signal strength limit",rLim
        rLim=rLimits['obs']
        h_obs.SetBinContent(h_obs.FindBin(x,y),rLim)
        h_obs_xs.SetBinContent(h_obs.FindBin(x,y),rLim*xs)
        gr["obs_xs"].SetPoint(i,x,y,rLim*xs)
        print "obs. signal strength limit",rLim

    f=rt.TFile(outdir+"limits_%s.root"%sScan,"update")
    for lvl in ["obs","obs+1","obs-1","exp","exp+1","exp-1","exp+2","exp-2","obs_xs"]:
        gr[lvl].Write("gr_"+lvl,rt.TObject.kOverwrite)

    grSig.Write("gr_"+"significance",rt.TObject.kOverwrite)
    if scan==Scan.T5gg: h_Sig.Write("h_"+"significance",rt.TObject.kOverwrite)
    if scan==Scan.T5Wg: h_Sig.Write("h_"+"significance",rt.TObject.kOverwrite)
    if scan==Scan.T6gg: h_Sig.Write("h_"+"significance",rt.TObject.kOverwrite)
    if scan==Scan.T6Wg: h_Sig.Write("h_"+"significance",rt.TObject.kOverwrite)
    if scan==Scan.GGM: h_Sig.Write("h_"+"significance",rt.TObject.kOverwrite)    

    h_exp.Write("h_exp",rt.TObject.kOverwrite)
    h_exp_xs.Write("h_exp_xs",rt.TObject.kOverwrite)
    h_obs.Write("h_obs",rt.TObject.kOverwrite)
    h_obs_xs.Write("h_obs_xs",rt.TObject.kOverwrite)
    f.Close()

def getContours(scan):
    sScan="unkown_scan"
    if scan==Scan.T5gg: sScan="T5gg"
    elif scan==Scan.T5Wg: sScan="T5Wg"
    elif scan==Scan.T6Wg: sScan="T6Wg"
    elif scan==Scan.T6gg: sScan="T6gg"
    elif scan==Scan.GGM:  sScan="GGM"

    f=rt.TFile(outdir+"limits_%s.root"%sScan,"update")
    for lvl in ["obs","obs+1","obs-1","exp","exp+1","exp-1","exp+2","exp-2"]:
        gr= f.Get("gr_"+lvl)
        grC=getContour(gr)
        if grC:
            grC.Write("gr_"+lvl+"C",rt.TObject.kOverwrite)
        else:
            print "could not get contour for",lvl
    f.Close()

def redoHistogram(scan):
    sScan="unkown_scan"
    if scan==Scan.T5gg: sScan="T5gg"
    elif scan==Scan.T5Wg: sScan="T5Wg"
    elif scan==Scan.T6Wg: sScan="T6Wg"
    elif scan==Scan.T6gg: sScan="T6gg"
    elif scan==Scan.GGM:  sScan="GGM"

    f=rt.TFile(outdir+"limits_%s.root"%sScan,"update")
    gr=f.Get("gr_obs_xs")
    h=gr.GetHistogram()
    h.Write("h_obs_xs_redone",rt.TObject.kOverwrite)
    f.Close()

def smoothContours(scan):
    sScan="unkown_scan"
    if scan==Scan.T5gg: sScan="T5gg"
    elif scan==Scan.T5Wg: sScan="T5Wg"
    elif scan==Scan.T6Wg: sScan="T6Wg"
    elif scan==Scan.T6gg: sScan="T6gg" 
    elif scan==Scan.GGM:  sScan="GGM"
    f=rt.TFile(outdir+"limits_%s.root"%sScan,"update")
    gs=rt.TGraphSmooth()
    for lvl in ["obs","obs+1","obs-1","exp","exp+1","exp-1","exp+2","exp-2"]:
        gr=f.Get("gr_"+lvl+"C")
        # gr_sm=gs.SmoothLowess(gr)
        # gr_sm=gs.Approx(gr)
        gr_sm=gs.SmoothSuper(gr)
        gr_sm.Write("gr_"+lvl+"C_sm",rt.TObject.kOverwrite)
    f.Close()

if __name__ == '__main__':
    basedir="/user/jschulz/2016/photonmet/"
    outdir=basedir+"output/"
    signal_scan="signal_scan_v19.root"
    rho=-0.0
    
#    caclulateLimits(Scan.TChiWg)
    
#    caclulateLimits(Scan.TChiNg)
    
 #   caclulateLimits(Scan.T5gg)
 #   getContours(Scan.T5gg)
 #   redoHistogram(Scan.T5gg)
 #   smoothContours(Scan.T5gg)

    caclulateLimits(Scan.T5Wg)
 #   getContours(Scan.T5Wg)
 #   redoHistogram(Scan.T5Wg)
 #   smoothContours(Scan.T5Wg)

 #   caclulateLimits(Scan.T6Wg)
 #   getContours(Scan.T6Wg)
 #   redoHistogram(Scan.T6Wg)
 #   smoothContours(Scan.T6Wg)

#    caclulateLimits(Scan.T6gg)
#    getContours(Scan.T6gg)
#    redoHistogram(Scan.T6gg)
#    smoothContours(Scan.T6gg)

#    caclulateLimits(Scan.GGM)
#    getContours(Scan.GGM)
#    smoothContours(Scan.GGM)



