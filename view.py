#!/usr/bin/env python2

import os,site,sys,time
import subprocess as sp
import ROOT as rt
thisFileDir=os.path.dirname(__file__)
site.addsitedir(thisFileDir+"/src/style")
rt.gROOT.LoadMacro(thisFileDir+"/src/style/tdrstyle.C")
rt.gROOT.Macro(thisFileDir+"/src/style/myStyle.C")

#site.addsitedir("/user/jschulz/2016/photonmet/src/style")
#rt.gROOT.LoadMacro("/user/jschulz/2016/photonmet/src/style/tdrstyle.C")
#rt.gROOT.Macro("/user/jschulz/2016/photonmet/src/style/myStyle.C")

# -r = remote:
# first download file
if len(sys.argv) > 1 and sys.argv[1]=="-r":
    sp.call(['scp','pz:ma/analysis/output/plots.root','output/plots.root'])
# open file in TBrowser
files=[]
for f in os.listdir(thisFileDir+"/output"):
    if f.endswith(".root"):files.append(rt.TFile("output/"+f,"read"))
br=rt.TBrowser()
if len(sys.argv) > 1 and "-p" in sys.argv:
    # persistent
    while True: time.sleep(100)
else:
    raw_input('...')
for f in files: f.Close()
