# HMS_optics

Version 3.0

I made some modifications so that this code can be compiled and run on JLab machines. Probably some version of clang changed since the original writing. Here's how to compile:
```
 cd source
 mkdir build
 cd build
 cmake ..
 make
```
This makes the executables and puts everything into the build directory. A quick tool just to get yourself running:
```
 cd build
 ./hms_optics setup_optics_example.txt -o outputFile.root -a
```
This requires a config file which should be in the format as in the setup_optics_example.txt. Then you define the output root file. The new matrix elements (xTar independent and xTar dependent) are printed out at the end of the code. If there are memory or TBranch errors, you will need a consolidated root tree to work with. Including all of the variables will throw this error. 

See also https://github.com/brash99/HMS_optics to see the format for the configuration file.
