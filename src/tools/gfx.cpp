#include "gfx.hpp"

#include "Config.hpp"
#include "io.hpp"
#include "util.hpp"
#include "style/cms_labels.hpp"

#include <string>
#include <stdexcept>
#include <TGraph.h>
#include <THStack.h>
#include <TMultiGraph.h>
#include <TPaletteAxis.h>
#include <TStyle.h>
#include <TText.h>


void gfx::getDrawnAxes(TPad const &pad,TAxis*&out_xax,TAxis*&out_yax)
{
   for (TObject const *obj: *(pad.GetListOfPrimitives())){
      if (obj->InheritsFrom(TH1::Class())){
         out_xax=((TH1*)obj)->GetXaxis();
         out_yax=((TH1*)obj)->GetYaxis();
         break;
      }
      if (obj->InheritsFrom(TGraph::Class())){
         out_xax=((TGraph*)obj)->GetXaxis();
         out_yax=((TGraph*)obj)->GetYaxis();
         break;
      }
      if (obj->InheritsFrom(TMultiGraph::Class())){
         // TODO
         throw "Not implemented!";
         break;
      }
      if (obj->InheritsFrom(THStack::Class())){
         TH1* h=((THStack*)obj)->GetHistogram();
         out_xax=h->GetXaxis();
         out_yax=h->GetYaxis();
         break;
      }
   }
}

void gfx::setupAxes(TAxis &x,TAxis &y)
{
 //  x.CenterTitle();
 //  y.CenterTitle();
   TString unit=parseAxis(x);
   parseAxis(y);

   TString yt(y.GetTitle());
   if (yt.Contains("BINW")){
      if (unit.Length()==0) unit="unit";
      yt.ReplaceAll("BINW"," / "+unit);
      yt="#LT"+yt+"#GT";
      y.SetTitle(yt);
   }
   if (yt.Contains("BIN")){
      if (unit.Length()==0) unit="units";
      unit=TString::Format(" / %g ",x.GetBinWidth(1))+unit;
      yt.ReplaceAll("BIN",unit);
      y.SetTitle(yt);
   }
}

void gfx::setupDrawnAxes(SplitCan const &spcan)
{
   TAxis *x,*y;
   // upper pad
   float f=spcan.fraction_;
   gfx::getDrawnAxes(spcan.pU_, x, y);
   assert(x&&y);
   x->SetTitleSize(0);
   x->SetLabelSize(0);
   float size,off;

   // lower pad
   f=1-f;
   gfx::getDrawnAxes(spcan.pL_, x, y);
   assert(x&&y);
   size=gStyle->GetTitleSize("y")/f;
   off=gStyle->GetTitleOffset("y")*f;
   y->SetTitleSize(size);
 //  y->SetTitleSize(0.55);
   y->SetTitleOffset(off);
   size=gStyle->GetLabelSize("y")/f;
   off=gStyle->GetLabelOffset("y")*f;
   y->SetLabelSize(size);
   y->SetLabelOffset(0.008);
   x->SetLabelOffset(0.02);
   
   size=gStyle->GetTitleSize("x")/f;
   x->SetTitleSize(size);
   size=gStyle->GetLabelSize("x")/f;
   x->SetLabelSize(size);

   size=gStyle->GetTickLength("x")/f;
   x->SetTickLength(size);

   setupDrawnAxes(spcan.pL_);
}

void gfx::setupDrawnAxes(TPad const &pad)
{
   TAxis *x=0,*y=0;
   gfx::getDrawnAxes(pad,x,y);
   assert(x && y);
   setupAxes(*x,*y);
}

void gfx::setupZaxis(TH2 &h2, bool log){
   gPad->Update();
   TPaletteAxis* palAx=(TPaletteAxis*)h2.GetListOfFunctions()->FindObject("palette");
   if (!palAx){
      debug<<TString::Format("no z axis found for '%s'",h2.GetName());
      return;
   }
   // thinner:
   palAx->SetX2NDC(palAx->GetX2NDC()-.03);

   TAxis *ax=h2.GetZaxis();
   ax->CenterTitle();
   ax->SetTickLength(.01);
   if (!log) ax->SetLabelOffset(.005);
}

TString gfx::parseAxis(TAxis &ax)
{
   auto newTitle=parseTitle(ax.GetTitle());
   ax.SetTitle(newTitle.first);
   return newTitle.second;
}

std::pair<TString,TString> gfx::parseTitle(TString const &title)
{
   // (quantity name,displayed quantity, unit)
   std::vector<std::tuple<TString,TString,TString>> repl={
      // with units
      std::make_tuple("VPT","#vec{p}_{#scale[.8]{T}}","GeV")
      ,std::make_tuple("PT","p_{#scale[.8]{T}}","GeV")
      ,std::make_tuple("METSHT","%MET / %HT^{#scale[.8]{1/2}}","GeV^{#scale[.8]{1/2}}")
      ,std::make_tuple("METSIG","#font[52]{S}","")
      ,std::make_tuple("VMET","#vec{E}_{#scale[.8]{T}}^{#scale[.8]{miss}}","GeV")
      ,std::make_tuple("%MET","p_{#scale[.8]{T}}^{#scale[.8]{miss}}","GeV")
      //~ ,std::make_tuple("%HTG","H_{#scale[.8]{T}}^{#scale[.8]{gen}}","GeV")
      //~ ,std::make_tuple("HT","H_{#scale[.8]{T}}","GeV")
      ,std::make_tuple("%HTG","H_{#scale[.8]{T}}^{#scale[.8]{#gamma}}","GeV")
      ,std::make_tuple("MT","M_{#scale[.8]{T}}","GeV")
      ,std::make_tuple("STg","S_{#scale[.8]{T}}^{#scale[.8]{#gamma}}","GeV")
      ,std::make_tuple("SIEIE","#sigma_{i#etai#eta}","")
      ,std::make_tuple("SIPIP","#sigma_{i#phii#phi}","")
      ,std::make_tuple("TREFF","Trigger Efficiency","")
   };

   TString nt(title);
   TString unit;
   for (auto r:repl){
      if (nt.Contains("%"+std::get<0>(r))){
         nt.ReplaceAll("%"+std::get<0>(r),std::get<1>(r));
      }
      if (nt.Contains(std::get<0>(r))){
         nt.ReplaceAll(std::get<0>(r),std::get<1>(r));
         if (std::get<2>(r).Length()>0){
            nt+=TString::Format(" (%s)",std::get<2>(r).Data());
            unit=std::get<2>(r);
         }
      }
   }
   return {nt,unit};
}

void gfx::decorate(TPad &pad, bool simulation, bool lumiText)
{
   gfx::setupDrawnAxes(pad);
   if (!Config::get().releaseMode) drawVersionLabel();
   style::draw_lumi(pad,simulation,lumiText);
}

void gfx::drawVersionLabel(float x, float y)
{
   if (y<1) y=1-gPad->GetTopMargin();
   Config const&cfg=Config::get();
   TString version;
   if (cfg.processFraction<1.0)
      version+=TString::Format("%.1f%% ",cfg.processFraction*100);
   version+="("+cfg.treeVersion+" "+cfg.gitHash+")";
   TText l;
   l.SetTextColor(kGray+2);
   l.SetNDC();
   l.SetTextFont(52);
   l.SetTextSize(.025);
   l.SetTextAlign(31);
   l.SetTextAngle(90);
   l.DrawText(x,y,version);
}

TLatex gfx::cornerLabel(TString text, int pos, float textSize,float offset){
   float x,y;
   int textAlign;
   // 1=top left, 2=top right, 3=bottom left, 4=bottom right
   switch (pos){
   case 1:{x=  gPad->GetLeftMargin() +offset;y=1-gPad->GetTopMargin()   -offset,textAlign=13;break;};
   case 2:{x=1-gPad->GetRightMargin()-offset;y=1-gPad->GetTopMargin()   -offset,textAlign=33;break;};
   case 3:{x=  gPad->GetLeftMargin() +offset;y=  gPad->GetBottomMargin()+offset,textAlign=11;break;};
   case 4:{x=1-gPad->GetRightMargin()-offset;y=  gPad->GetBottomMargin()+offset,textAlign=31;break;};
   default: throw std::invalid_argument("pos");
   }
   TLatex l(x,y,text);
   l.SetNDC();
   l.SetTextFont(42);
   l.SetTextSize(textSize);
   l.SetTextAlign(textAlign);
   return l;
}

TPaveText gfx::paveText(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option)
{
   TPaveText pt(x1,y1,x2,y2,option);
   pt.SetBorderSize(0);
   pt.SetTextFont(42);
   pt.SetFillStyle(0);
   return pt;
}

TLegend gfx::mkLegend(float x1,float y1,float x2,float y2,int nCols,float offset)
{
   if (x2<0){
      x2=x1;
      if (x1 <.5) x1=gPad->GetLeftMargin()+offset;
      else        x2=1-gPad->GetRightMargin()-offset;
   }
   if (y2<0){
      y2=y1;
      if (y1 <.5) y1=gPad->GetBottomMargin()+offset;
      else        y2=1-gPad->GetTopMargin()-offset;
   }
   TLegend l(x1,y1,x2,y2);
   l.SetBorderSize(0);
   l.SetFillStyle(0);
   l.SetNColumns(nCols);
   return l;
}

void gfx::removeXerrors(TGraphErrors &gr)
{
   for (int i=0; i<gr.GetN(); i++)
      gr.SetPointError(i,0.,gr.GetErrorY(i));
}
void gfx::scaleXerrors(TGraphErrors &gr, float factor)
{
   for (int i=0; i<gr.GetN(); i++)
      gr.SetPointError(i,gr.GetErrorX(i)*factor,gr.GetErrorY(i));
}
void gfx::setYvalues(TGraphErrors &gr, float value)
{
   Double_t *x=gr.GetX();
   for (int i=0; i<gr.GetN(); i++)
      gr.SetPoint(i,x[i],value);
}


gfx::LegendEntries::~LegendEntries()
{
   for (TLine *l: lDummy_) delete l;
}

void gfx::LegendEntries::append(TObject const &obj, const char *label, Option_t *option)
{
   lEntr_.push_back(Entry{&obj,label,option});
}
void gfx::LegendEntries::append(const char *label)
{
   lEntr_.push_back(Entry{nullptr,label,""});
}

void gfx::LegendEntries::append(Width_t lwidth,Style_t lstyle,Color_t col,const char *label, Option_t *option)
{
   lDummy_.push_back(new TLine());
   TLine *l=lDummy_.back();
   l->SetLineWidth(lwidth);
   l->SetLineStyle(lstyle);
   l->SetLineColor(col);
   lEntr_.push_back(Entry{l,label,option});
}

void gfx::LegendEntries::prepend(TObject const &obj, const char *label, Option_t *option)
{
   lEntr_.push_front(Entry{&obj,label,option});
}

gfx::LegendEntries& gfx::LegendEntries::operator+=(LegendEntries const &other)
{
   lEntr_.insert(lEntr_.end(),other.lEntr_.begin(),other.lEntr_.end());
   return *this;
}

void gfx::LegendEntries::clear()
{
   lEntr_.clear();
}

void gfx::LegendEntries::pop_front()
{
   lEntr_.pop_front();
}

void gfx::LegendEntries::pop_back()
{
   lEntr_.pop_back();
}

TLegend gfx::LegendEntries::buildLegend(float x1,float y1,float x2,float y2, int nCols,Option_t *option)
{
   TLegend l=mkLegend(x1,y1,x2,y2,nCols);
   l.SetOption(option);
   for (Entry const &e: lEntr_)
      l.AddEntry(e.obj,e.label,e.opt);

   return l;
}


int gfx::SplitCan::s_count=0;
gfx::SplitCan::SplitCan(float fraction)
   : fraction_(fraction)
   , can_("splitCan"+TString::Itoa(s_count,10))
   , pU_("padU"+TString::Itoa(s_count,10),"padU",0,0,1,1)
   , pL_("padL"+TString::Itoa(s_count,10),"padL",0,0,1,1-fraction)
{
   SplitCan::s_count+=1;
   pU_.SetBottomMargin(1-fraction);
   pU_.Draw();
   pU_.Update();
   pL_.SetTopMargin(0);
   pL_.SetBottomMargin(gStyle->GetPadBottomMargin()/(1-fraction));
   pL_.SetFillColor(0);
   pL_.SetFillStyle(0);
   pL_.SetGridy();
   pL_.Draw();
   pL_.Update();
   cdUp();
}
