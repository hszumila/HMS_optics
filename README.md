HMS_optics
==========

HMS Optics Optimization using ROOT compiled macro
-------------------------------------------------

The main program is HMS_ytarget_fit.C. It runs as a standard compiled ROOT macro: e.g.: 

.L HMS_ytarget_fit.C+

or

.x HMS_ytarget_fit.C+("setup_filename.txt")

There is one required argument and several more optional arguments: 

void HMS_ytarget_fit(const char* configfilename, const char *outrootfilename = "HMS_optics_fit_output.root", int niter_xtarcorr=1, double xsieve_offset_for_cut_centering=0.3, double ysieve_offset_for_cut_centering = 0.0)

The first argument is the configuration file name (required).

The second argument is the name of a ROOT file where a diagnostic ROOT tree is written. 

The third argument is the number of extra iterations of reconstruction to perform to get
the "xtarget" correction right.  You need at least one so that you know 
the vertex position so you can compute "xtar", the vertical position of the intersection 
of the spectrometer ray with the plane perpendicular to the HMS optical axis containing 
the origin. 

The other two optional arguments are offsets (in cm) to the x (vertical) and 
y (horizontal) positions of the sieve slit holes. These are not used to offset 
the "true" sieve hole positions, but are used to offset the center of the range of 
reconstructed track coordinates at the sieve slit used to assign events to each hole.

Configuration File Specfication
-------------------------------

The code parses the configuration file for a list of "runs", each of which is associated
with a given set of conditions such as beam position, HMS central angle, number of target 
foils, z position of the foils, cuts, etc. 

The code reads the file until the "endlist" keyword is encountered. The lines are expected 
to be in the form "keyword value". Possible keywords are: 

newrun XXXXX: This starts a new "run" definition. The "run number" XXXXX, which might or might not correspond with a CODA run number, must be unique or the properties will be overwritten. Everything read between this line and the next line starting with "newrun" will be assigned to this run. Whenever this keyword is read, all possible properties of the run are initialized with sensible defaults. 

filelist file1.root file2.root file*.root: This line contains a list of ROOT files containing hms ntuples associated with this run. Wildcards are accepted. This command can be used multiple times per run, subsequent invocations add to the file list. 

beampos x y dx/dz dy/dz: Expects four numbers after the "beampos" keyword. Defines the beam position on target and beam angles. Currently, the beam angles are not used and assumed to be zero.

thetaHMS thetadeg: "thetaHMS" keyword expects to be followed by a single parameter, the HMS central angle in degrees. HMS is always assumed to be on beam right. 

nfoil N: "nfoil" keyword defines number of target foils in the run. Expects one integer argument.

zfoil z1 z2 ... zN: "zfoil" keyword defines z position of target foils in the run. Expects nfoil arguments, assumed to be given in cm. If nfoil is not yet defined for the run, the command has no effect. 

sieveslit 1: "sieveslit" flag defines whether the sieve slit was in or out for this run. Expects an integer argument (0 or 1). If you have thin-foil data without sieve slit, the code will attempt to include the data in the ytarget fit, but generally you will only want to use data with sieve slit to fit ytarget, xptar and yptar, so usually you want to use "sieveslit 1". You can omit this command and it will default to 1. 

cut cut1 cut2 ... cutN: Defines event selection cuts, allows multiple cut definitions separated by spaces. Any expression that defines a valid TCut can be used here. 

gcut file gcut1 .... gcutN: Defines graphical cuts. The first argument is the name of a ROOT file containing TCutG objects and subsequent arguments are the names of TCutG objects to be loaded from the file. Note that the x and y variable names for the TCutG objects must match the names of variables in the ROOT tree for graphical cuts to work. 

In the case of keywords beampos, thetaHMS, nfoil, zfoil and sieveslit, if the keyword appears more than once, the last invocation supersedes any previous ones. In the case of filelist, cut, and gcut, subsequent invocations add files, TCut and TCutG objects to the list of files, cuts, gcuts for the run in question.  

After the "endlist" keyword is encountered, there are a few subsequent arguments expected. 

oldcoeff.dat: name of file containing "old" (starting) reconstruction matrix elements
newcoeff.dat: name of file containing "new" reconstruction matrix elements (file to which fit results will be written). 
htheta_offset_old: offset applied to yptar (in-plane) angle.
hphi_offset_old: offset applied to xptar (out-of-plane) angle.
fitorder: usually 5 or 6, this is the order of the fit to perform.
maxnperhole: integer, max number of events per sieve hole per target foil to include in the fit (generally a small number like 100-1000 events is good here, to ensure roughly equal weighting of sieve holes in the fit). 
maxnperfoil: max number of events per target foil to include in the fit (only applies to the "sieveslit 0" case, so generally has no effect). 
zoffset_foil: Global offset applied to the z position of all target foils. 
fit_xtar_coeffs_flag: 0 or 1, flag to choose whether to include xtarget-dependent coefficients in the fit (generally better to use 0 as these are computed from COSY, since in real data the system is underdetermined). 

In general, the code will then loop over all sieve holes for each foil for each run, and attempt to automatically determine a cut to select events going through a given sieve hole, but will require input from the user to validate each hole (you can tweak the range of the fit to get a better result, or you can reject a hole if a good fit cannot be achieved). Once cuts have been defined for all sieve holes for all foils for all runs, the program sets up and solves the equations for the new coefficients, and writes the file. The ROOT tree in the output root file contains diagnostic information, basically it has the focal plane track parameters, the "true" target track parameters determined from the sieve hole positions and the target foil/beam positions and the reconstructed target track parameters using the old (initial) reconstruction matrix elements. 
 
