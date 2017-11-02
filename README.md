## Requirements ##
- ROOT6
- C++11 compiler
- CMake >= 2.8 (lx-cluster: cmake28 executable)
- Boost

## Build ##
    mkdir build; cd build
    cmake ..
    make

On lx-cluster:

- use `CMSSW_7_4_X`
- use `cmake28` instead of `cmake`

## Run ##
Execute `run.x` (created in the build directory).
Without any arguments, a test module is run.
Meaningful modules can be chosen as arguments:

    $ run.x --help
    Usage: ./run.x module1[ module2[...]] [options]
    Allowed options:
      -h [ --help ]              produce help message
      -f [ --fraction ] arg (=1) Fraction of events to process
      --release                  Release mode (don't draw version labels)


## Output ##
Plots are stored as canvases in a root file.
The `view.py` script can be used for browsing the file,
it sets the right styles.

## Adding new modules ##
To add a new module, create a new `<module>.cpp` file in `src/modules`.
The entry function has to be `extern "C" void run() {...}`.
Add `<module>` to the list in `src/modules/CMakeLists.txt`


## PDF and Scale Weights ##
index description-from-LHEheader
0   <weight id="1"> mur=1 muf=1 </weight>
1   <weight id="2"> mur=1 muf=2 </weight>
2   <weight id="3"> mur=1 muf=0.5 </weight>
3   <weight id="4"> mur=2 muf=1 </weight>
4   <weight id="5"> mur=2 muf=2 </weight>
5   <weight id="6"> mur=2 muf=0.5 </weight>
6   <weight id="7"> mur=0.5 muf=1 </weight>
7   <weight id="8"> mur=0.5 muf=2 </weight>
8   <weight id="9"> mur=0.5 muf=0.5 </weight>

10  <weight id="11">Member 1</weight>
11  <weight id="12">Member 2</weight>
12  <weight id="13">Member 3</weight>
13  <weight id="14">Member 4</weight>
...
