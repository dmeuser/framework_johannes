/* template implementation file for io */

template <class T>
T* io::RootFileReader::read(TString name) const
{
   if (intPath_!="") name=intPath_+"/"+name;
   T* h=(T*)file_->Get(name);
   if (!h){
      debug<<TString::Format("did not find '%s' in '%s'",(intPath_+"/"+name).Data(),fName_.Data());
      throw;
   }
   return h;
}

