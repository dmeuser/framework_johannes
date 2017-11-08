#Tool to create datacards

import os, re, subprocess as sp, math
from distutils import spawn
from array import array
from operator import add,sub
import ROOT as rt

backgrounds=["GJ","TTcomb","Vg","diboson","efake"]

class Scan:
    T5gg,T5Wg,T6Wg,T6gg,GGM,TChiWg,TChiNg=range(7)

def prepareDatacard(obs,count,estat,esyst,eISR,correlated,mcUncertainties,pointName,sScan):
    dataCardPath=outdir+"datacards/"+sScan+"/datacard_%s.txt"%pointName
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
    elif scan==Scan.GGM:    pattern=".*_M2_(.*)_M1_(.*)"
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


def fillDatacards(scan):
    print "using",os.environ["CMSSW_VERSION"]
    print

    scanFile=basedir
    if scan==Scan.T5gg:   scanFile+="T5gg_scan.txt"
    elif scan==Scan.T5Wg:   scanFile+="T5Wg_scan.txt"
    elif scan==Scan.T6gg:   scanFile+="T6gg_scan.txt"
    elif scan==Scan.T6Wg:   scanFile+="T6Wg_scan.txt"
    elif scan==Scan.GGM:    scanFile+="GGM_WinoBino_scan.txt"
    elif scan==Scan.TChiWg: scanFile+="TChiWg_scan.txt"
    elif scan==Scan.TChiNg: scanFile+="TChiNg_scan.txt"
    
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
                #if m2<1100: continue
            points.append(p)

    for i,point in enumerate(points):
        print point,
        m2,m1=getMasses(point,scan)
        key=m2
        if scan==Scan.GGM: key=m2*100000+m1
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
        datacard=prepareDatacard(obs,c,st,sy,sISR,cor,mcUncertainties,point,sScan)


if __name__ == '__main__':
    basedir="../"
    outdir=basedir+"output/"
    signal_scan="signal_scan_v19.root"
    rho=-0.0
    #fillDatacards(Scan.T5gg)
    fillDatacards(Scan.T5Wg)
    """
    fillDatacards(Scan.T6Wg)
    fillDatacards(Scan.T6gg)
    fillDatacards(Scan.GGM)
    fillDatacards(Scan.TChiWg)
    fillDatacards(Scan.TChiNg)
    """

