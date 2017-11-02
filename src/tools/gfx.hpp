#ifndef GFX_HPP__
#define GFX_HPP__

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPaveText.h>

#include <list>

namespace gfx
{
   class SplitCan
   {
      static int s_count;
   public:
      SplitCan(float fraction=0.65);
      void cdUp(){ pU_.cd(); }
      void cdLow(){ pL_.cd(); }
      void Update(){ pU_.Update(); pL_.Update(); }

      float fraction_;
      TCanvas can_;
      TPad pU_,pL_;
   };

   /* find the axes of the currently drawn object in the pad */
   void getDrawnAxes(TPad const &pad,TAxis*&out_xax,TAxis*&out_yax);
   /* format the axes of the */
   void setupAxes(TAxis &x,TAxis &y);
   /* format the drawn axes of the pad */
   void setupDrawnAxes(TPad const &pad);
   /* the same for a split canvas */
   void setupDrawnAxes(SplitCan const &spcan);

   /* Setup color palette z-axis for 2D histogram.
    * - Centering title
    * - Smaller width of the palette axis
    */
   void setupZaxis(TH2 &h2, bool log=true);
   /* parse the axis title (substring replacements), returns unit string */
   TString parseAxis(TAxis &ax);
   /* parse a title string (substring replacements), returns (title, unit)*/
   std::pair<TString,TString> parseTitle(TString const &title);

   /* decorate pad (axes parsing, labels) */
   void decorate(TPad &pad, bool simulation=true, bool lumiText=true);

   /* draw version label (tree version + git hash) */
   void drawVersionLabel(float x=0.995, float y=-1);

   /* label placed in the corner of the drawn axes' area
    * pos: 1=top left, 2=top right, 3=bottom left, 4=bottom right
    */
   TLatex cornerLabel(TString text, int pos=2, float textSize=0.04,float offset=.02);
   TPaveText paveText(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="ndc");

   /* Create TLegend. If only two coordinates are given, the others are set automatically */
   TLegend mkLegend(float x1,float y1,float x2=-1,float y2=-1,int nCols=1,float offset=.02);

   /* graph modifiers for displaying errors in ratios etc. */
   void removeXerrors(TGraphErrors &gr);
   void scaleXerrors(TGraphErrors &gr, float factor);
   void setYvalues(TGraphErrors &gr, float value);

   /* class for managing Legend entries */
   class LegendEntries
   {
   public:
      ~LegendEntries();
      void append(TObject const &obj, const char *label="", Option_t *option="lpf");
      void append(const char *label=""); // for label-only or empty slots
      void append(Width_t lwidth,Style_t lstyle,Color_t col,const char *label="", Option_t *option="l"); // append a line
      void prepend(TObject const &obj, const char *label="", Option_t *option="lpf");
      LegendEntries& operator+=(LegendEntries const &other);

      void clear();
      void pop_front();
      void pop_back();

      TLegend buildLegend(float x1,float y1,float x2=-1,float y2=-1,int nCols=1,Option_t *option="brNDC");
   private:
      struct Entry{
         TObject const* obj;
         TString label;
         TString opt;
      };
      std::list<Entry> lEntr_;

      std::vector<TLine*> lDummy_;
   };
} //namespace gfx

#endif /* GFX_HPP__ */
