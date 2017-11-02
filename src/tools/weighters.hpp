#ifndef WEIGHTERS_HPP
#define WEIGHTERS_HPP

// Abstract class for weight calculation (possibly) depending on different variables
class WeightCalculator
{
public:
   virtual float get() const=0;
protected:
   WeightCalculator(){}
};

class WeightCalculatorZgamma : public WeightCalculator
{
public:
   WeightCalculatorZgamma(float const &pt):pt_(pt){}
   float get() const override {
      // take from EXO-16-014
      // http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2016_078_v9.pdf
      if (pt_ < 175) return 1.39; // this was not calculated...
      if (pt_ < 190) return 1.39;
      if (pt_ < 250) return 1.35;
      if (pt_ < 400) return 1.30;
      if (pt_ < 700) return 1.23;
      return 1.23;
   }
private:
   float const &pt_;
};

class WeightCalculatorWgamma : public WeightCalculator
{
public:
   float get() const override {
      // take from EXO-16-014
      // http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2016_078_v9.pdf
      // was not calculated for pt<175...
      return 1.34;
   }
};


#endif /* WEIGHTERS_HPP */
