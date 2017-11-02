#include "util.hpp"

std::string util::rm_duplicate_spaces(std::string s)
{
   // remove leading/trailing space
   boost::trim(s);
   // replace duplicate spaces with single ones
   while (s.find( "  " ) != s.npos) boost::replace_all(s,"  "," ");
   return s;
}


