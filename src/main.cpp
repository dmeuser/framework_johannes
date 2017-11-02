#include "Config.hpp"
#include "style/myStyle.hpp"

#include "tools/io.hpp"

#include <string>
#include <iostream>
#include <chrono>
#include <dlfcn.h>
#include <boost/program_options.hpp>

#include <TROOT.h>

static int parseCLIopt(int argc, char* argv[]);
static int runModule(std::string module);

int main(int argc, char* argv[]){
   gROOT->SetBatch();
   mystyle();

   parseCLIopt(argc,argv);

   Config const &cfg=Config::get();
   int exitValue;
   for (std::string m: cfg.modules){
      exitValue=runModule(m);
      if (exitValue!=0)return exitValue;
   }
   return 0;
}

static int parseCLIopt(int argc, char* argv[])
{
   Config &cfg=Config::get();

   // specifiy command line options
   namespace po = boost::program_options;
   po::options_description vis_desc("Allowed options");
   po::options_description hid_desc("Hidden options"); // for positional arguments
   vis_desc.add_options()
      ("help,h", "produce help message")
      ("fraction,f", po::value<float>(&cfg.processFraction)->default_value(1.0), "Fraction of events to process")
      ("release",  po::bool_switch(&cfg.releaseMode)->default_value(false), "Release mode (don't draw version labels)")
      ;
   hid_desc.add_options()
      ("modules", po::value<std::vector<std::string>>(), "Modules to run")
      ;
   po::positional_options_description pos_o;
   pos_o.add("modules", -1);

   // combining visible and hidden options
   po::options_description desc;
   desc.add(vis_desc).add(hid_desc);

   po::variables_map vm;
   po::store(po::command_line_parser(argc, argv).
             options(desc).positional(pos_o).run(), vm);
   po::notify(vm);

   if (vm.count("help")) {
      std::cout << "Usage: "<<argv[0]<<" module1[ module2[...]] [options]" << std::endl
                << vis_desc << std::endl;
      return 1;
   }

   cfg.modules.clear();
   if (vm.count("modules")){
      for (std::string m: vm["modules"].as<std::vector<std::string> >()) {
         cfg.modules.push_back(m);
      }
   }else{
      // run tests, if nothing specified
      cfg.modules.push_back("tests");
   }

   return 0;
}

static int runModule(std::string module)
{
   std::string libPath(LIBRARY_OUTPUT_PATH);
   libPath+="lib"+module+".so";
   void* handle = dlopen(libPath.c_str(), RTLD_LAZY);
   if (!handle) {
      std::cerr << "Cannot open library: " << dlerror() << std::endl
                << "Choose an available module" << std::endl;
      return 1;
   }
   typedef void (*run_t)();
   run_t run=(run_t)dlsym(handle,"run");
   const char *error=dlerror();
   if (error){
      std::cerr << "Cannot load symbol 'run': " << error << std::endl;
      dlclose(handle);
      return 1;
   }

   Config const &cfg=Config::get();
   io::log<<TString::Format("== running '%s' (%s %s) ==",module.c_str(),cfg.treeVersion.Data(),cfg.gitHash.Data());
   // measure time for running module
   namespace chr = std::chrono;
   chr::steady_clock::time_point tStart = chr::steady_clock::now();
   // execute the run function of the module
   (*run)();
   chr::steady_clock::time_point tEnd = chr::steady_clock::now();
   chr::duration<double> duration = chr::duration_cast<chr::duration<double>>(tEnd - tStart);
   // io::log*"Took"/duration.count()<<"s";
   io::log<<TString::Format("== %s (%.1f%%) took %.2fs ==",module.c_str(),cfg.processFraction*100,duration.count());
   dlclose(handle);
   return 0;
}
