# HMS_optics
HMS Optics Optimization using ROOT compiled macro

The main program is HMS_ytarget_fit.C. It runs as a standard compiled ROOT macro: e.g.: 

.L HMS_ytarget_fit.C+

or

.x HMS_ytarget_fit.C+("setup_filename.txt")

There is one required argument and several more optional arguments: 

void HMS_ytarget_fit(const char* configfilename, const char *outrootfilename = "HMS_optics_fit_output.root", int niter_xtarcorr=1, double xsieve_offset_for_cut_centering=-0.28, double ysieve_offset_for_cut_centering = 0.0){

The first argument is the configuration file name (required).

The second argument is the name of a ROOT file where a diagnostic ROOT tree is written. 

The third argument is the number of extra iterations of reconstruction to perform to get the "xtarget" correction right.  You need at least one so that you know the vertex position so you can compute "xtar", the vertical position of the intersection of the spectrometer ray with the plane perpendicular to the HMS optical axis containing the origin. 

The other two optional arguments are offsets (in cm) to the x (vertical) and y (horizontal) positions of the sieve slit holes. These are not used to offset the "true" sieve hole positions, but are used to offset the center of the range of reconstructed track coordinates at the sieve slit used to assign events to each hole. 
