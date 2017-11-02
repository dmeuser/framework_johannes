#ifndef UTIL_HPP__
#define UTIL_HPP__

#include "io.hpp"

#include <cmath>
#include <string>
#include <sstream>
#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace util
{
   /* convert a comma-separated string to a vector of T */
   template<typename T>
   std::vector<T> to_vector(const std::string& s, char delim=',',bool eraseWhitespace=true)
   {
      std::vector<T> result;
      std::stringstream ss(s);
      std::string item;
      while(std::getline(ss, item, delim)){
         if (eraseWhitespace) boost::erase_all(item, " ");
         result.push_back(boost::lexical_cast<T>(item));
      }
      return result;
   }

   std::string rm_duplicate_spaces(std::string s);

   /* adds list of values quadratically */
   template<typename floatT>
   floatT quadSum(std::vector<floatT> values)
   {
      floatT qs=0;
      for (floatT v:values) qs+=v*v;
      return qs;
   }
   template<typename floatT>
   floatT sqrtQuadSum(std::vector<floatT> values)
   {
      return sqrt(quadSum(values));
   }
}

#endif /* UTIL_HPP__ */

