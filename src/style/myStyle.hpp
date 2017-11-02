#ifndef MYSTYLE_HPP__
#define MYSTYLE_HPP__

#include <TROOT.h>

void mystyle()
{
   TString stylePath(CMAKE_SOURCE_DIR);
   stylePath+="src/style/";
   gROOT->LoadMacro(stylePath+"tdrstyle.C");
   gROOT->Macro(stylePath+"myStyle.C");
}

#endif /* MYSTYLE_HPP__ */
