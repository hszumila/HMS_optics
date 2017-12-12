# HMS_optics -- from Jure's recode

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
 ./hms_optics ../data/setup_optics_z1cm_noxbeam_simple.txt -o outputFile
```
This requires a config file which should be in the format as in the setup_optics_z1cm_noxbeam_simple.txt. Then you define the output root file. The new matrix elements (xTar independent and xTar dependent) are printed out at the end of the code. For now, I've noticed that I have problems when using Root 6 and later. This should be fixed eventually.

See also https://github.com/brash99/HMS_optics
