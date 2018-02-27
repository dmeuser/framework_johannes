/*
 * Create Datacards
 */

#include "Config.hpp"
#include "tools/io.hpp"

#include <cstdlib>
#include <boost/python.hpp>

namespace py = boost::python;

static Config const &cfg=Config::get();

extern "C"
void run()
{
   Py_Initialize();

   try {
      py::object main_module = py::import("__main__");
      py::object main_namespace = main_module.attr("__dict__");
      main_namespace["basedir"]= CMAKE_SOURCE_DIR;
      main_namespace["outdir"] = (CMAKE_SOURCE_DIR+cfg.outputDirectory).Data();
      main_namespace["signal_scan"] = TString::Format("signal_scan_%s.root",cfg.treeVersion.Data()).Data();
      main_namespace["rho"] = cfg.sf.rho;
      py::exec_file((CMAKE_SOURCE_DIR+std::string("src/py/datacards.py")).c_str(),main_namespace);
   } catch (py::error_already_set const &) {
      assert(PyErr_Occurred());
      PyErr_Print();
   }
}
