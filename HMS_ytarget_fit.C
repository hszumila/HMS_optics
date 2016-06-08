#include "HMS_tree.h"
#include "HMS_tree.C"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TEventList.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TCut.h"
#include "TCutG.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TSystem.h"
#include "TDecompSVD.h"
#include "TPolyMarker.h"
#include "TList.h"
#include "TSpectrum.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLine.h"
#include "TPostScript.h"
#include "TMarker.h"
#include "TEllipse.h"
//#include <cstdlib>
//#include "unistd.h"
#include <stdio.h>

using namespace std;

const double yminhole = -6.096;
const double xminhole = -10.160;
const double yspacehole = 1.524;
const double xspacehole = 2.540;

//The following numbers are from HMS collimator survey dated Aug. 11, 2008.
const double y0sieve = 0.0; // y position of center sieve hole wrt HMS optical axis.
const double x0sieve = 0.0; // relative to HMS optical axis
const double xt0sieve = +0.311; // center sieve hole relative to ideal target center line. 
const double z0sieve = 166.032;

//The following numbers are from the HMS pointing survey: this is where the spectrometer points:
const double y_mispointing_HMS = 0.0;
const double x_mispointing_HMS = 0.14;

//octagon dimensions:
const double dx_coll = 11.646;
const double dy_coll = 4.575;
const double x0_coll = 0.106;
const double y0_coll = -0.013;


const int nrow=9;
const int ncol=9;

void HMS_ytarget_fit(const char* configfilename, const char *outrootfilename = "HMS_optics_fit_output.root", int niter_xtarcorr=1, double xsieve_offset_for_cut_centering=0.3, double ysieve_offset_for_cut_centering = 0.0){

  ofstream logfile("HMS_optics_fit_log.txt");
  
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadRightMargin(0.18);
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetTitleOffset(1.2,"Y");

  gStyle->SetNdivisions(505,"XYZ");
  
  gStyle->SetStatX(1);
  gStyle->SetStatY(1);
  gStyle->SetStatW(0.09);
  gStyle->SetStatH(0.3);
  
  ifstream setupfile(configfilename);
  
  TString currentline;

  //We want to define a list of runs, and with each file we associate a run number, a global cut and a list of foils
  vector<int> RunList; //run numbers
  map<int,vector<TString> > FileList; //list of files, can include wildcards.
  map<int,double> xbeam,xpbeam,ybeam,ypbeam; //Beam position at z = 0 and slope
  map<int,double> thetacentral; //central HMS angle
  map<int,int> nfoil; //number of target foils in this run
  map<int,vector<double> > zfoil; //z position of target foils in this run
  map<int,int> sieve_slit_flag; //Flag indicating presence or absence of sieve slit
  map<int,TCut> Cuts; //mapping of TCuts and run numbers (might add graphical cuts later)
  map<int,vector<TCutG*> > GCuts; //mapping of graphical cuts and run numbers
  //map<int,TEventList *> EventList; //mapping of run numbers and TEventLists 
  map<int,map<int,double> > xmean_hole, ymean_hole, sigx_hole, sigy_hole; //mean xy positions and sigx, sigy defining cuts to select events for each sieve hole:
  map<int,map<int,int> > use_hole; //flag to use hole or not use hole, by run number:
  map<int,map<int,int> > nevents_hole; //count number of events per foil/sieve hole.
  map<int,vector<double> > zmin_foil, zmax_foil, ymin_foil, ymax_foil; 
  map<int,vector<int> > use_foil; //flag by run, foil number to use data for that foil
  map<int,vector<int> > nevents_foil; //count number of events per foil, regardless of sieve hole
  
  int nruns = 0;
  int runnum = 0;
  
  while( currentline.ReadLine(setupfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      int ntokens = tokens->GetEntries();

      TString skeyword = ( (TObjString*) (*tokens)[0] )->GetString();
      
      if( skeyword == "newrun" && ntokens >= 2 ){
	TString srun = ( (TObjString*) (*tokens)[1] )->GetString();
	RunList.push_back( srun.Atoi() );
	runnum = srun.Atoi();
	nruns++;

	//Initialize properties of this run to sensible defaults:
	Cuts[runnum] = "";
	FileList[runnum].clear();
	xbeam[runnum] = ybeam[runnum] = xpbeam[runnum] = ypbeam[runnum] = 0.0;
	thetacentral[runnum] = 90.0*TMath::Pi()/180.0;
	nfoil[runnum] = 0;
	zfoil[runnum].clear();
	sieve_slit_flag[runnum] = 1; //assume sieve slit by default;
	TString elistname;
	elistname.Form("elist%d",runnum);
	//EventList[runnum] = new TEventList( elistname );
	GCuts[runnum].clear();
      }

      if( nruns > 0 ){
      
	if( skeyword == "filelist" && ntokens >= 2 ){
	  for( int i=1; i<ntokens; i++ ){
	    TString sfilename = ( (TObjString*) (*tokens)[i] )->GetString();
	    if( nruns > 0 ){
	      FileList[runnum].push_back( sfilename );
	    }
	  }
	}

	if( skeyword == "beampos" && ntokens >= 5 ){
	  TString sxbeam = ( (TObjString*) (*tokens)[1] )->GetString();
	  TString sybeam = ( (TObjString*) (*tokens)[2] )->GetString();
	  TString sxpbeam = ( (TObjString*) (*tokens)[3] )->GetString();
	  TString sypbeam = ( (TObjString*) (*tokens)[4] )->GetString();

	  xbeam[runnum] = sxbeam.Atof();
	  ybeam[runnum] = sybeam.Atof();
	  xpbeam[runnum] = sxpbeam.Atof();
	  ypbeam[runnum] = sypbeam.Atof();
	}

	if( skeyword == "thetaHMS" && ntokens >= 2 ){
	  TString stheta = ( (TObjString*) (*tokens)[1] )->GetString();
	  thetacentral[runnum] = stheta.Atof() * TMath::Pi()/180.0;
	}

	if( skeyword == "nfoil" && ntokens >= 2 ){
	  TString snfoil = ( (TObjString*) (*tokens)[1] )->GetString();
	  nfoil[runnum] = snfoil.Atoi();
	}

	if( skeyword == "zfoil" && nfoil.find( runnum ) != nfoil.end() ){
	  if( ntokens >= nfoil[runnum] + 1 ){
	    for( int i=1; i<=nfoil[runnum]; i++ ){
	      TString szfoil = ( (TObjString*) (*tokens)[i] )->GetString();
	      zfoil[runnum].push_back( szfoil.Atof() );
	    }
	  }
	}

	if( skeyword == "sieveslit" && ntokens >= 2 ){
	  TString ssieve = ( (TObjString*) (*tokens)[1] )->GetString();
	  sieve_slit_flag[runnum] = ssieve.Atoi();
	}

	if( skeyword == "cut" && ntokens >= 2 ){
	  //assume cuts separated by spaces:
	  for(int i=1; i<ntokens; i++){
	    TString scut = ( (TObjString*) (*tokens)[i] )->GetString();
	    Cuts[runnum] += scut;
	  }
	}

	if( skeyword == "gcut" && ntokens >= 3 ){
	  //assume second argument is a file name and subsequent arguments are gcuts to be loaded from the file:
	  TString sgcutsfilename = ( (TObjString*) (*tokens)[1] )->GetString();

	  TFile *ftemp = new TFile( sgcutsfilename, "READ" );
	  for( int i=2; i<ntokens; i++ ){
	    TString sgcutname = ( (TObjString*) (*tokens)[i] )->GetString();
	    TCutG *gcut_temp;

	    ftemp->GetObject( sgcutname, gcut_temp ); 
	    if( gcut_temp ){
	      Cuts[runnum] += gcut_temp->GetName();
	      GCuts[runnum].push_back( new TCutG( *gcut_temp ) );
	    }
	  }
	}
      }
    }
  }
  
  double htheta_offset_old,hphi_offset_old;

  string oldcoeffsfilename;
  string newcoeffsfilename;

  setupfile >> oldcoeffsfilename;
  setupfile >> newcoeffsfilename;
  setupfile >> htheta_offset_old;
  setupfile >> hphi_offset_old;

  ifstream oldcoeffsfile(oldcoeffsfilename.c_str());
  ofstream newcoeffsfile(newcoeffsfilename.c_str());
  
  int num_recon_terms_old;

  vector<double> xptarcoeffs_old;
  vector<double> yptarcoeffs_old;
  vector<double> ytarcoeffs_old;
  vector<double> deltacoeffs_old;
  vector<int> xfpexpon_old;
  vector<int> xpfpexpon_old;
  vector<int> yfpexpon_old;
  vector<int> ypfpexpon_old;
  vector<int> xtarexpon_old;

  vector<double> xptarcoeffs_new;
  vector<double> yptarcoeffs_new;
  vector<double> ytarcoeffs_new;
  vector<double> deltacoeffs_new;
  vector<int> xfpexpon_new;
  vector<int> xpfpexpon_new;
  vector<int> yfpexpon_new;
  vector<int> ypfpexpon_new;
  vector<int> xtarexpon_new;

  TString headerline;
  
  //line-by-line copy of optics coeffs file header:

  while( currentline.ReadLine(oldcoeffsfile,kFALSE) && !currentline.BeginsWith(" ----") ){
    newcoeffsfile << currentline.Data() << endl;
  }
  newcoeffsfile << currentline.Data() << endl;

  headerline = currentline;

  //in old version of the data 
  
  num_recon_terms_old = 0;

  while( currentline.ReadLine(oldcoeffsfile,kFALSE) && !currentline.BeginsWith(" ----") ){
    cout << currentline.Data() << endl;
    //extract the coeffs and exponents from the line:

    double coeff[4]; 
    //int expontemp[5];
    char cexpontemp[5];
    
    sscanf( currentline.Data(), "%lg %lg %lg %lg %5s", &(coeff[0]), &(coeff[1]), &(coeff[2]), &(coeff[3]),
	    cexpontemp );

    TString sexpon = cexpontemp;
    
    int expontemp[5];
    for(int i=0; i<5; i++){
      TString si(sexpon(i,1));
      expontemp[i] = si.Atoi();
    }
    
    //for(int i=0; i<4; i++){
    // cout << "Old coefficents (theta,y,phi,delta)=(" << coeff[0] << ", " << coeff[1] << ", " << coeff[2] << ", " << coeff[3] << ", " << coeff[4] << ")" << endl;
    // 	//}

    // cout << "Old exponents (x,xp,y,yp,xt)=(" << expontemp[0] << ", " << expontemp[1] << ", " << expontemp[2] << ", " << expontemp[3] << ", " << expontemp[4] << ")" << endl;
    
    // TString sc1(currentline(1,16));
    // TString sc2(currentline(17,16));
    // TString sc3(currentline(33,16));
    // TString sc4(currentline(49,16));
    
    xptarcoeffs_old.push_back(coeff[0]);
    ytarcoeffs_old.push_back(coeff[1]);
    yptarcoeffs_old.push_back(coeff[2]);
    deltacoeffs_old.push_back(coeff[3]);
    
    

    // for(int expon=0; expon<5; expon++){
    //   TString stemp(currentline(66+expon,1));
    //   expontemp[expon] = stemp.Atoi();
    // }

    xfpexpon_old.push_back(expontemp[0]);
    xpfpexpon_old.push_back(expontemp[1]);
    yfpexpon_old.push_back(expontemp[2]);
    ypfpexpon_old.push_back(expontemp[3]);
    xtarexpon_old.push_back(expontemp[4]);

    // cout << "(C_theta, C_y, C_phi, C_delta) = (" 
    // 	 << xptarcoeffs_old[num_recon_terms_old] << ", "
    // 	 << ytarcoeffs_old[num_recon_terms_old] << ", " 
    // 	 << yptarcoeffs_old[num_recon_terms_old] << ", "
    // 	 << deltacoeffs_old[num_recon_terms_old] << "), expon = "
    // 	 << xfpexpon_old[num_recon_terms_old] << ", " << xpfpexpon_old[num_recon_terms_old] << ", " 
    // 	 << yfpexpon_old[num_recon_terms_old] << ", " << ypfpexpon_old[num_recon_terms_old] << ", "
    // 	 << xtarexpon_old[num_recon_terms_old] << endl;
    
    num_recon_terms_old++;
  }

  TString endline = currentline;
  
  cout << "num recon terms in OLD matrix = " << num_recon_terms_old << endl;

  TFile* fout = new TFile(outrootfilename,"RECREATE");

  fout->cd();

  TTree* tout = new TTree("tout","");

  int trunnum;
  int tifoil;
  double tzfoil;
  double txfp,tyfp,txpfp,typfp,txtar,txbeam,tybeam;
  double txptarrecon,txptartrue,txptarfit;
  double typtarrecon,typtartrue,typtarfit;
  double tytarrecon,tytartrue,tytarfit;
  double txtartrue;

  //track quantities at fp:
  tout->Branch("xfp",&txfp,"xfp/D");
  tout->Branch("yfp",&tyfp,"yfp/D");
  tout->Branch("xpfp",&txpfp,"xpfp/D");
  tout->Branch("ypfp",&typfp,"ypfp/D");
  //xtar is reconstructed from ybeam and "old" recon quantities:
  tout->Branch("xtar",&txtar,"xtar/D");
  //xtartrue is reconstructed from beam position, foil position and sieve hole position
  
  //"OLD" reconstructed quantities:
  tout->Branch("xptar",&txptarrecon,"xptar/D");
  tout->Branch("yptar",&typtarrecon,"yptar/D");
  tout->Branch("ytar",&tytarrecon,"ytar/D");
  // tout->Branch("xptarfit",&txptarfit,"xptarfit/D");
  // tout->Branch("yptarfit",&typtarfit,"yptarfit/D");
  // tout->Branch("ytarfit",&tytarfit,"ytarfit/D");
  //"TRUE" angles and vertex coordinates:
  tout->Branch("xtartrue",&txtartrue,"xtartrue/D");
  tout->Branch("ytartrue",&tytartrue,"ytartrue/D");
  tout->Branch("xptartrue",&txptartrue,"xptartrue/D");
  tout->Branch("yptartrue",&typtartrue,"yptartrue/D");
  //"TRUE" beam position, foil position and foil number
  tout->Branch("xbeam", &txbeam, "xbeam/D");
  tout->Branch("ybeam", &tybeam, "ybeam/D");
  tout->Branch("zfoil", &tzfoil, "zfoil/D");
  tout->Branch("ifoil", &tifoil, "ifoil/I");
  tout->Branch("runnum", &trunnum, "runnum/I");
  
  vector<double> xtartrue,ytartrue,xptartrue,yptartrue,deltatrue;
  vector<double> xfptrue,yfptrue,xpfptrue,ypfptrue;
  //  vector<double> xtarrecon_old, ytarrecon_old, xptarrecon_old, yptarrecon_old,
    
  //gROOT->ProcessLine(".x ~/rootlogon.C");
  gStyle->SetOptFit();

  int fitorder;
  int maxnperhole, maxnperfoil; //max events per sieve hole and max events per target foil (for runs without sieve slit)
  //double nsigcut,nsigcut_y; //cut tolerance:
  
  setupfile >> fitorder;
  setupfile >> maxnperhole >> maxnperfoil;
  //setupfile >> nsigcut >> nsigcut_y;

  double zoffset_foil = 0.0;
  setupfile >> zoffset_foil; //Assume that foils are located at nominal positions + zoffset_foil;

  int fit_xtar_coeffs_flag = 0;
  setupfile >> fit_xtar_coeffs_flag; //fit xtar coeffs for xptar only!
  
  cout << fitorder << endl;

  double xhole[nrow][ncol], yhole[nrow][ncol];

  //initialize positions of sieve holes. hole positions are expressed in global coordinates
  for( int row=0; row<nrow; row++ ){
    for( int col=0; col<ncol; col++ ){
      yhole[row][col] = yminhole + yspacehole * col + y0sieve + y_mispointing_HMS;
      xhole[row][col] = xminhole + xspacehole * row + x0sieve + x_mispointing_HMS;
    }
  }
  
  //loop over all runs:
  //Start with ytarget/zbeam histograms by run:
  
  TClonesArray *ytarget_histos = new TClonesArray( "TH1D", RunList.size() );
  TClonesArray *zbeam_histos = new TClonesArray( "TH1D", RunList.size() );
  TClonesArray *xysieve_histos = new TClonesArray( "TH2D", RunList.size() );
  
  TCanvas *c1 = new TCanvas("c1","First Canvas",200,10,700,500);
  TCanvas *c2 = new TCanvas("c2","Second Canvas",200,10,700,500);

  c1->Update();
  c2->Update();
  
  int ihisto=0;
  int ihistosieve=0;
  
  int nterms_fit=0; //we only fit the xtarget-independent terms, we hold the xtar matrix elements fixed at their values calculated by COSY:
 
  // The point of this whole loop is to a) intitialize the NEW xptar, yptar, and ytar coefficients to zero, and b) to set the NEW delta
  // coefficients to be the same as the OLD delta coefficients.
  //
  for( int i=0; i<=fitorder; i++ ){
    for( int j=0; j<=fitorder-i; j++){
      for( int k=0; k<=fitorder-i-j; k++){
	for( int l=0; l<=fitorder-i-j-k; l++){
	  for( int m=0; m<=fitorder-i-j-k-l; m++){
	    xptarcoeffs_new.push_back( 0.0 );
	    yptarcoeffs_new.push_back( 0.0 );
	    ytarcoeffs_new.push_back( 0.0 );
	    //deltacoeffs_new.push_back( 
	    //determine delta coefficient for these exponents:
	    
	    xfpexpon_new.push_back( m );
	    yfpexpon_new.push_back( k );
	    xpfpexpon_new.push_back( l );
	    ypfpexpon_new.push_back( j );
	    xtarexpon_new.push_back( i );

	    int delta_term = -1;
	    for( int term=0; term<num_recon_terms_old; term++){
	      if( xfpexpon_old[term] == m &&
		  yfpexpon_old[term] == k &&
		  xpfpexpon_old[term] == l &&
		  ypfpexpon_old[term] == j &&
		  xtarexpon_old[term] == i ){
		delta_term = term;
		break;
	      }
	    }

	    if( delta_term >= 0 ){
	      deltacoeffs_new.push_back( deltacoeffs_old[delta_term] );
	    } else {
	      deltacoeffs_new.push_back( 0.0 );
	    }
	  
	    nterms_fit++;
	  }
	}
      }
    }
  }

  cout << "number of terms in new matrix (not counting xtar-dependent terms) = " << nterms_fit << endl;
  
  TMatrixD Fit_matrix_xptarget( nterms_fit, nterms_fit );
  TMatrixD Fit_matrix_yptarget( nterms_fit, nterms_fit );
  TMatrixD Fit_matrix_ytarget( nterms_fit, nterms_fit );
  TVectorD ytarget_vec( nterms_fit );
  TVectorD yptarget_vec( nterms_fit );
  TVectorD xptarget_vec( nterms_fit );
  
  for( int i=0; i<nterms_fit; i++){
    for( int j=0; j<nterms_fit; j++){
      Fit_matrix_xptarget(i,j) = 0.0;
      Fit_matrix_yptarget(i,j) = 0.0;
      Fit_matrix_ytarget(i,j) = 0.0;
    }
    ytarget_vec(i) = 0.0;
    yptarget_vec(i) = 0.0;
    xptarget_vec(i) = 0.0;
  }
  
  for(int irun=0; irun<RunList.size(); irun++ ){
    TChain *C = new TChain("h9010");
   
    int runnum = RunList[irun];

    for( int ifile=0; ifile<FileList[runnum].size(); ifile++ ){
      cout << "Adding file " << FileList[runnum][ifile].Data() << " to TChain" << endl;
      C->Add( FileList[runnum][ifile] );
    }

    cout << "Number of files added = " << C->GetNtrees() << endl;
    cout << "Number of events before cut = " << C->GetEntries() << endl;
    
    TString drawcmd;
    drawcmd.Form( ">>elist%d",runnum);

    //gDirectory->ls();
    //fout->ls();
    
    cout << "Draw selection command = " << drawcmd.Data() << endl;
    cout << "Cut expression = " << Cuts[runnum] << endl;
    
    C->Draw( drawcmd.Data(), Cuts[runnum] );

    //cout << "Number of events passing cut = " << EventList[runnum]->GetN() << endl;

    TString elist_name;
    elist_name.Form("elist%d",runnum);
    
    TEventList *elist_temp = (TEventList*) gDirectory->Get(elist_name);
    cout << "Number of events passing cut = " << elist_temp->GetN() << endl;
    
    HMS_tree* T = new HMS_tree(C);

    //Make histograms per foil? 

    TString histname;
    histname.Form( "ytarget_run%d", runnum );
    new( (*ytarget_histos)[ihisto] ) TH1D( histname, "", 250, -7.5,12.5 );

    ( (TH1D*) (*ytarget_histos)[ihisto] )->GetXaxis()->SetTitle("y_{target} (cm)");
    
    histname.Form( "zbeam_run%d", runnum );
    new( (*zbeam_histos)[ihisto] ) TH1D( histname, "", 250, -12.5, 19.0 );
    
    ( (TH1D*) (*zbeam_histos)[ihisto] )->GetXaxis()->SetTitle("z_{vertex} (cm)");

    for( int ifoil=0; ifoil<nfoil[runnum]; ifoil++ ){
      zfoil[runnum][ifoil] += zoffset_foil;
      
      histname.Form( "xysieve_run%d_foil%d", runnum, ifoil );
      new( (*xysieve_histos)[ihistosieve+ifoil] ) TH2D( histname, "", 250, -7.2, 7.2, 250, -12.5,12.5 );

      ( (TH2D*) (*xysieve_histos)[ihistosieve+ifoil] )->GetXaxis()->SetTitle("y_{sieve} (cm)");
      ( (TH2D*) (*xysieve_histos)[ihistosieve+ifoil] )->GetYaxis()->SetTitle("x_{sieve} (cm)");
    }
    
    //histname.Form( "xysieve_run%d", runnum );
    
    
    double thetaHMS = thetacentral[runnum];
    double x0beam = xbeam[runnum]; //+xbeam = beam right
    double y0beam = ybeam[runnum]; //+ybeam = vertical up

    //Correct equation for vertical beam position is ybeam = y0beam - fry
    // Correct equation for xbeam in spectrometer coordinates is xtarbeam = -ybeam = -y0beam + fry

    int ipass=0;

    for( ipass=0; ipass<2; ipass++ ){
      if( ipass == 1 && sieve_slit_flag[runnum] == 0 ) break;
      long nevents_run = 0;
      while( T->GetEntry( elist_temp->GetEntry( nevents_run++ ) ) ){ //This is the event loop:
	if( nevents_run % 10000 == 0 ) cout << nevents_run << endl;
	//reconstruct angles and vertex using the old coefficients:
	double ytarsum = 0.0;
	double yptarsum = 0.0;
	double xptarsum = 0.0;
	for( int term=0; term<num_recon_terms_old; term++ ){
	  //In reconstruction, we want xtarget to be with respect to the spectrometer axis, so to get xtarget in spectrometer coordinates from x in global coordinates,
	  //we SUBTRACT the x mispointing of the spectrometer:
	  double lambda =
	    pow( T->hsxfp / 100.0, xfpexpon_old[term] ) *
	    pow( T->hsyfp / 100.0, yfpexpon_old[term] ) *
	    pow( T->hsxpfp, xpfpexpon_old[term] ) *
	    pow( T->hsypfp, ypfpexpon_old[term] ) *
	    pow( (T->fry_cm - y0beam - x_mispointing_HMS)/100.0, xtarexpon_old[term] );
	  
	  ytarsum += ytarcoeffs_old[term] * lambda;
	  yptarsum += yptarcoeffs_old[term] * lambda;
	  xptarsum += xptarcoeffs_old[term] * lambda;
	}
	
	double ytarold = ytarsum * 100.0 + y_mispointing_HMS; //This ytarget is with respect to the origin; what is measured by the spectrometer is relative to HMS optical axis
	double yptarold = yptarsum + htheta_offset_old;
	double xptarold = xptarsum + hphi_offset_old;
	
	//The beam is assumed to have scattered from the point (xbeam + xpbeam * zfoil, ybeam + ypbeam*zfoil, zfoil)
	//This is the intersection of the horizontal component of the spectrometer ray with the nominal beam line.
	
	//ytarget = yvertex - yptar * zspec 
	// yvertex = zvertex * sin(theta) - xbeam * cos(theta)
	// zspec = zvertex * cos(theta) + xbeam * sin(theta)
	//ytarget = zvertex * sin(theta) - xbeam * cos(theta) - yptar*(zvertex*cos(theta) + xbeam*sin(theta))
	// zvertex*(sin(theta)-yptar*cos(theta)) = ytarget + xbeam*(cos(theta)+yptar*sin(theta))
	
	//A beam shift in the +x direction (beam right) shifts the z of the intersection point with the beam line by 
	double zvertexold = ( ytarold + (x0beam-T->frx_cm ) * ( cos(thetaHMS) + yptarold * sin(thetaHMS) ) )/ ( sin(thetaHMS) - yptarold * cos(thetaHMS) );
	double zspecold = zvertexold * cos(thetaHMS) + (x0beam-T->frx_cm) * sin(thetaHMS);
	
	double xtarold = T->fry_cm - y0beam - x_mispointing_HMS - zspecold * xptarold;
	
	// So, the basic idea seems to be as follows:
	// 1.  In the code just above, the old coefficients are used to 
	//     calculate a corrected estimate for xtar, based on calculations
	//     of xptar and zpec (which depends on the beam position).
	//
	// 2.  Now, armed with this new estimate of xtar, we will do at least one more iteration.
	//     The default is niter_xtarcorr = 1.  But, if do more than one iteration, presumably
	//     we will get increasingly convergent estimates for xtarold.
	
	for( int iter=0; iter<niter_xtarcorr; iter++ ){
	  ytarsum = 0.0;
	  yptarsum = 0.0;
	  xptarsum = 0.0;
	  for( int term=0; term<num_recon_terms_old; term++ ){
	    
	    double lambda =
	      pow( T->hsxfp / 100.0, xfpexpon_old[term] ) *
	      pow( T->hsyfp / 100.0, yfpexpon_old[term] ) *
	      pow( T->hsxpfp, xpfpexpon_old[term] ) *
	      pow( T->hsypfp, ypfpexpon_old[term] ) *
	      pow( xtarold/100.0, xtarexpon_old[term] );
	    
	    ytarsum += ytarcoeffs_old[term] * lambda;
	    yptarsum += yptarcoeffs_old[term] * lambda;
	    xptarsum += xptarcoeffs_old[term] * lambda;
	    
	  }
	  
	  ytarold = ytarsum*100.0 + y_mispointing_HMS;
	  yptarold = yptarsum + htheta_offset_old;
	  xptarold = xptarsum + hphi_offset_old;
	  
	  zvertexold = ( ytarold + (x0beam-T->frx_cm) * ( cos(thetaHMS) + yptarold * sin(thetaHMS) ) )/ ( sin(thetaHMS) - yptarold * cos(thetaHMS) );
	  zspecold = zvertexold * cos(thetaHMS) + (x0beam-T->frx_cm)*sin(thetaHMS);
	  
	  xtarold = T->fry_cm - y0beam - x_mispointing_HMS - zspecold * xptarold;
	}
	
	//at this point, xtarold is expressed relative to the HMS optical axis, ytarold is expressed in global coordinates, and sieve hole positions are expressed in
	//global coordinates.
	xtarold += x_mispointing_HMS;
	// and now, after the above statement, xtarold is expressed in global coordinates.

	if( ipass==0 ){ //fill ytarget/zbeam histograms:
	  ( (TH1D*) (*ytarget_histos)[ihisto] )->Fill( ytarold );
	  ( (TH1D*) (*zbeam_histos)[ihisto] )->Fill( zvertexold );
	} else {
	  for( int ifoil=0; ifoil<nfoil[runnum]; ifoil++){
	    if( use_foil[runnum][ifoil] == 1 &&
		zmin_foil[runnum][ifoil] <= zvertexold && zvertexold <= zmax_foil[runnum][ifoil] &&
		ymin_foil[runnum][ifoil] <= ytarold && ytarold <= ymax_foil[runnum][ifoil] ){
	      ( (TH2D*) (*xysieve_histos)[ihistosieve+ifoil] )->Fill( ytarold + yptarold * z0sieve, xtarold + xptarold * z0sieve );
	    }
	  }
	}
      }
      if( ipass == 0 ){ //determine target foil cuts:
	TLine L;
	L.SetLineWidth(2);
	L.SetLineColor(2);

	TH1D *htemp;
    
	for( int foil=0; foil<nfoil[runnum]; foil++ ){
	  //First, loop over xy sieve histograms for all runs with sieve slit:
	  c1->Clear();
	  //c1->Divide(2,1);
	  c1->SetGrid();
	  //c1->cd(2)->SetGrid();
	  
	  double zfoil_fit_spacing = 1.0; // Default fit spacing for single foil
	  
 	  if (foil == 0) {
		if (nfoil[runnum] >= 2) {
			zfoil_fit_spacing = fabs((zfoil[runnum][1]-zfoil[runnum][0])/2.0);
		}
	  }
	  
	  htemp = (TH1D*) (*zbeam_histos)[ihisto];
	  htemp->Draw();
	  L.DrawLine( zfoil[runnum][foil], 0.0, zfoil[runnum][foil], ( (TH1D*) (*zbeam_histos)[ihisto] )->GetMaximum() );
	  c1->Update();
	  c1->Modified();
	  cout << "Define z range to fit foil " << foil + 1 << ", nominal z = " << zfoil[runnum][foil] << ":" << endl;
	  double zmintemp, zmaxtemp;
	  zmintemp = zfoil[runnum][foil]-zfoil_fit_spacing;
	  zmaxtemp = zfoil[runnum][foil]+zfoil_fit_spacing;
	  //cout << "zmin = ";
	  //cin >> zmintemp;
	  //cout << "zmax = ";
	  //cin >> zmaxtemp;
	    
	  htemp->GetXaxis()->SetRangeUser(zmintemp,zmaxtemp);
	    
	  int binmax = htemp->GetMaximumBin(), binlow = binmax, binhigh = binmax;
	  //fit +/- HM of peak:
	  while( htemp->GetBinContent(binlow--) >= 0.5*htemp->GetMaximum() && binlow >= 1 ){};
	  while( htemp->GetBinContent(binhigh++) >= 0.5*htemp->GetMaximum() && binhigh <= htemp->GetNbinsX()){};
	  double zlowfit = htemp->GetBinCenter( binlow );
	  double zhighfit = htemp->GetBinCenter(binhigh );
	    
	  htemp->Fit("gaus","","",zlowfit,zhighfit);
	  c1->Update();
	  cout << "good fit (y/n)?";
	  TString reply="";
	  reply.ReadLine(cin,kTRUE);
	  if( reply.BeginsWith("n") ){
	    use_foil[runnum].push_back(0);
	  } else {
	    use_foil[runnum].push_back(1);
	    nevents_foil[runnum].push_back(0);
	  }
	  
	  htemp->GetXaxis()->SetRangeUser( -12.5, 19.0 );
	  
	  TF1 *functemp = (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus") );
	  
	  cout << "z peak - z foil = " << functemp->GetParameter(1) - zfoil[runnum][foil] << endl;
	  
	  zmin_foil[runnum].push_back( functemp->GetParameter(1) - 3.0*functemp->GetParameter(2) );
	  zmax_foil[runnum].push_back( functemp->GetParameter(1) + 3.0*functemp->GetParameter(2) );
	  
	  //ytarget:
	  
	  htemp = (TH1D*) (*ytarget_histos)[ihisto];
	  htemp->Draw();
	  c1->Update();

	  double ytarget_foil_position;
	  ytarget_foil_position = zfoil[runnum][foil]*sin(thetaHMS);

	  cout << "Define ytarget cuts for foil " << foil + 1 << endl;
	  
	  double ytarget_fit_spacing = 1.0; // Default fit spacing for single foil
	  
 	  if (foil == 0) {
		if (nfoil[runnum] >= 2) {
			ytarget_fit_spacing = fabs((zfoil[runnum][1]-zfoil[runnum][0])/2.0)*sin(thetaHMS);
		}
	  }

	  L.DrawLine( ytarget_foil_position, 0.0, ytarget_foil_position, ( (TH1D*) (*ytarget_histos)[ihisto] )->GetMaximum() );
	  c1->Update();
	  c1->Modified();
	  double ymintemp,ymaxtemp;
	  ymintemp=ytarget_foil_position-ytarget_fit_spacing;
	  ymaxtemp=ytarget_foil_position+ytarget_fit_spacing;
	  //cout << "ymin=";
	  //cin >> ymintemp;
	  //cout << "ymax=";
	  //cin >> ymaxtemp;

	  htemp->GetXaxis()->SetRangeUser( ymintemp,ymaxtemp );
	  binmax = htemp->GetMaximumBin(); binlow=binmax; binhigh=binmax;
	  while( htemp->GetBinContent(binlow--) >= 0.5*htemp->GetMaximum() && binlow >= 1 ){};
	  while( htemp->GetBinContent(binhigh++) >= 0.5*htemp->GetMaximum() && binhigh <= htemp->GetNbinsX() ){};
	  double ylowfit = htemp->GetBinCenter(binlow);
	  double yhighfit = htemp->GetBinCenter(binhigh);
	  htemp->Fit("gaus","","",ylowfit,yhighfit);

	  c1->Update();
	  cout << "good fit (y/n)";
	  reply.ReadLine(cin,kTRUE);
	  if( reply.BeginsWith("n") ){
	    use_foil[runnum][foil] = 0;
	  }
	  functemp = (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus") );
      
	  ymin_foil[runnum].push_back( functemp->GetParameter(1) - 3.0*functemp->GetParameter(2) );
	  ymax_foil[runnum].push_back( functemp->GetParameter(1) + 3.0*functemp->GetParameter(2) );
	  htemp->GetXaxis()->SetRangeUser( -7.5, 12.5 );

	  logfile << "Run " << runnum << ", foil " << foil << ", (zfoil,zmin,zmax,ymin,ymax)=(" << zfoil[runnum][foil] << ", " << zmin_foil[runnum][foil] << ", " << zmax_foil[runnum][foil] << ", "
		  << ymin_foil[runnum][foil] << ", " << ymax_foil[runnum][foil] << ")" << endl;
	}
    
	c1->Update();
	
      } else { //determine sieve hole cuts:

    //For y target fit, plot zbeam histograms and ytarget histograms; let user define a range of zbeam and ytarget to select for the fitting
    //For runs *with* sieve slit, compute ytrue and yptrue from hole position, foil position and beam position.
    //For runs without sieve slit, compute ytrue from foil position and assume that yptrue = yprecon. Perhaps for this reason we should fit the angles first for runs with sieve slit only. In either case, we need to select foils/holes:

    
	TH1D *htemp;
	TH2D *htemp2D;
    
	if( sieve_slit_flag[runnum] == 1 ){ //determine sieve hole cuts for this run:
	  bool auto_fitting = true;
 	  bool first = true;

	  for(int ifoil=0; ifoil<nfoil[runnum]; ifoil++){

	    bool continue_getting_sieve_holes = true;

	    while (continue_getting_sieve_holes) {

	    TEllipse etemp[nrow][ncol]; // array of ellipses representing sieve hole events to be used in fits.
	    TMarker M[nrow][ncol]; // array of markers representing sieve slit hole positions.
	    for( int col=0; col<ncol; col++ ){
	      for( int row=0; row<nrow; row++ ){
	    
		if( !(row==2&&col==3) && !(row==5&&col==5) ){

		  //Let's create a new function for each fit:

		  TString funcname;
		  funcname.Form("gauss2D_foil%d_row%d_col%d",ifoil,row,col);
		  
		  //We don't want ROOT to have any memory of the previous fit parameters of this function during the next fit:
		  TF2 *gauss2D = new TF2(funcname,"[0]*exp(-0.5*(pow((x-[1])/[2],2)+pow((y-[3])/[4],2)))",-7.5,7.5,-12.5,12.5);
		  
		  c1->Clear();
		  htemp2D = ( (TH2D*) (*xysieve_histos)[ihistosieve+ifoil] );

		  // htemp2D->Draw("col");

		  // c1->Update();
		  
		  double nevents_avg_per_hole = htemp2D->GetEntries() / double( nrow*ncol );
		  
		  //double xholetemp = xhole[row][col] + xsieve_offset_for_cut_centering;
		  //double yholetemp = yhole[row][col] + ysieve_offset_for_cut_centering;
		  double xholetemp = xhole[row][col] + xsieve_offset_for_cut_centering;
		  double yholetemp = yhole[row][col] + ysieve_offset_for_cut_centering;
	    
		  double startpar[5] = {1.0,yholetemp,0.2,xholetemp,0.3};
		  double starterr[5] = {TMath::Max(1.0,sqrt(nevents_avg_per_hole/3.0)),0.1,0.05,0.1,0.1}; 
		  // for(int ipar=0; ipar<5; ipar++){ //might be fixed later:
		  //   gauss2D.ReleaseParameter(ipar);
		  // }
		  //gauss2D.SetParErrors(starterr);
		  gauss2D->SetParameters(startpar);
		  gauss2D->SetParLimits( 2,0.01,yspacehole/5.0 );
		  gauss2D->SetParLimits( 4,0.01,xspacehole/5.0 );
		  gauss2D->SetParLimits( 1,yholetemp-1.0*yspacehole, yholetemp+1.0*yspacehole ); //won't work well if not already approximately calibrated
		  gauss2D->SetParLimits( 3,xholetemp-1.0*xspacehole, xholetemp+1.0*xspacehole ); 

		  gauss2D->Print();
		  
		  //gauss2D->SetRange( yholetemp-0.3*yspacehole, xholetemp-0.3*xspacehole,
		  //		       yholetemp+0.3*yspacehole, xholetemp+0.3*xspacehole );

		  int bin00 = htemp2D->FindBin( yholetemp-0.5*yspacehole, xholetemp-0.5*xspacehole );
		  int bin01 = htemp2D->FindBin( yholetemp-0.5*yspacehole, xholetemp+0.5*xspacehole );
		  int bin10 = htemp2D->FindBin( yholetemp+0.5*yspacehole, xholetemp-0.5*xspacehole );
		  int bin11 = htemp2D->FindBin( yholetemp+0.5*yspacehole, xholetemp+0.5*xspacehole );
	  
		  int binxlo,binylo,binxhi,binyhi,binz;
		  htemp2D->GetBinXYZ( bin00, binxlo, binylo, binz );
		  htemp2D->GetBinXYZ( bin11, binxhi, binyhi, binz );
	    
		  double xlo = yholetemp - 0.6*yspacehole;
		  double ylo = xholetemp - 0.6*xspacehole;
		  double xhi = yholetemp + 0.6*yspacehole;
		  double yhi = xholetemp + 0.6*xspacehole;

		  double nevents_hole_temp = htemp2D->Integral( binxlo, binxhi, binylo, binyhi );
		  
		  if( nevents_hole_temp >= 50 ){ //don't attempt a fit if number of events < 1000

		    bool fitok = false, aborthole=false;

		    TF1 *functemp;

		    TLine L;
		    L.SetLineColor(4);
		    L.SetLineWidth(2);

		    c1->cd();
		    htemp2D->Draw("col");
		    c1->Update();
	            
		    // Determine whether manual or automatic fitting is desired.
	
  		    if (first) { 
                       	first = false;
			cout << "Manual (m) or Automatic (a) Fitting?" << endl;
		        TString response;
        	    	response.ReadLine(cin);
		    	if((response.BeginsWith("m")||response.BeginsWith("M"))){
				auto_fitting = false;
		    	}
		    }
	 
		    while( !fitok ){ //fit 1D projections of x and y along y and x:
		      htemp = htemp2D->ProjectionX( "htemp_px", binylo, binyhi );
		      c2->cd();
		      htemp->GetXaxis()->SetRangeUser( yholetemp - 1.0*yspacehole, yholetemp + 1.0*yspacehole );
		      htemp->Draw();
		      htemp->Fit("gaus","","",xlo,xhi);
		      //draw vertical line at nominal hole position:
		      L.DrawLine( yhole[row][col], 0.0, yhole[row][col], htemp->GetMaximum() );
		      c2->Update();
		      TString reply;
		      if (auto_fitting){
			reply = "y";
		      }else{
		        cout << "fit ok (y/n/skip)?";
			reply.ReadLine(cin);
	              }

		      if( !(reply.BeginsWith("n")||reply.BeginsWith("s") )){
			fitok = true;
			functemp = (TF1*) htemp->GetListOfFunctions()->FindObject("gaus");
			gauss2D->SetParameter(1,functemp->GetParameter(1));
			gauss2D->SetParameter(2,TMath::Max(0.01,TMath::Min(fabs(functemp->GetParameter(2)),yspacehole/5.0)));
			gauss2D->SetParError(1,functemp->GetParError(1));
			gauss2D->SetParError(2,functemp->GetParError(2));

			if( functemp->GetParameter(0) < 100.0 ) gauss2D->FixParameter(2,TMath::Max(0.01,TMath::Min(fabs(functemp->GetParameter(2)),yspacehole/5.0)));
			
		      } else if (reply.BeginsWith("s") ){
			aborthole = true;
			fitok = true;
		      } else {
			cout << " xmin fit = " << endl;
			cin >> xlo;
			cout << "xmax fit = " << endl;
			cin >> xhi;
		      }
		    }

		    if( aborthole ) continue;
	    
		    fitok = false;
		    while( !fitok ){
		      htemp = htemp2D->ProjectionY( "htemp_py", binxlo, binxhi );
		      c2->cd();
		      htemp->GetXaxis()->SetRangeUser( xholetemp-1.5*xspacehole, xholetemp + 1.5*xspacehole );
		      htemp->Draw();
		      htemp->Fit("gaus","","",ylo,yhi);
		      L.DrawLine( xhole[row][col], 0.0, xhole[row][col], htemp->GetMaximum() );
		      c2->Update();
		      TString reply;
		      if (auto_fitting){
			reply = "y";
		      }else{
		        cout << "fit ok (y/n/skip)?";
			reply.ReadLine(cin);
	              }

		      if( !(reply.BeginsWith("n")||reply.BeginsWith("s") )){
			fitok = true;
			functemp = (TF1*) htemp->GetListOfFunctions()->FindObject("gaus");
			gauss2D->SetParameter(3,functemp->GetParameter(1));
			gauss2D->SetParameter(4,TMath::Max(0.01,TMath::Min(fabs(functemp->GetParameter(2)),xspacehole/5.0)));
			gauss2D->SetParError(3,functemp->GetParError(1));
			gauss2D->SetParError(4,functemp->GetParError(2));
			if( functemp->GetParameter(0) < 100.0 ) gauss2D->FixParameter(4,TMath::Max(0.01,TMath::Min(fabs(functemp->GetParameter(2)),xspacehole/5.0)));
		      } else if (reply.BeginsWith("s") ){
			aborthole = true;
			fitok = true;
		      } else {
			cout << " ymin fit = " << endl;
			cin >> ylo;
			cout << "ymax fit = " << endl;
			cin >> yhi;
		      }
		    }
	      
		    if( aborthole ) continue;

		    c1->cd();

		    if( !aborthole ){
		      gauss2D->SetRange( gauss2D->GetParameter(1)-2.5*gauss2D->GetParameter(2), gauss2D->GetParameter(3)-2.5*gauss2D->GetParameter(4),
					 gauss2D->GetParameter(1)+2.5*gauss2D->GetParameter(2), gauss2D->GetParameter(3)+2.5*gauss2D->GetParameter(4) );

		  
	      
		      htemp2D->Fit(gauss2D,"0R");
		      
		      etemp[row][col].SetFillStyle(0);
		      etemp[row][col].SetLineWidth(2);
		      etemp[row][col].SetLineColor(1);
		
		      M[row][col].SetMarkerStyle(21);
		      M[row][col].SetMarkerColor(kMagenta);
		      M[row][col].SetMarkerSize(1.5);
		      M[row][col].DrawMarker(yholetemp,xholetemp );
		
		      if( fabs( gauss2D->GetParameter(1) - yholetemp ) <= 1.0*yspacehole &&
			  fabs( gauss2D->GetParameter(3) - xholetemp ) <= 1.0*xspacehole &&
			  fabs( gauss2D->GetParameter(2) ) <= yspacehole/5.0 &&
			  fabs( gauss2D->GetParameter(4) ) <= xspacehole/5.0 ){
			  etemp[row][col].SetX1(gauss2D->GetParameter(1));
			  etemp[row][col].SetY1(gauss2D->GetParameter(3));
			  etemp[row][col].SetR1(2.5*gauss2D->GetParameter(2));
			  etemp[row][col].SetR2(2.5*gauss2D->GetParameter(4));
			  etemp[row][col].Draw("SAME");
			//etemp[row][col].DrawEllipse( gauss2D->GetParameter(1), gauss2D->GetParameter(3),
			//		   2.5*gauss2D->GetParameter(2),
			//		   2.5*gauss2D->GetParameter(4), 0, 360.0, 0.0, "SAME" );
			c1->Update();
		        TString reply;
		        if (auto_fitting){
			  reply = "y";
		        }else{
		          cout << "fit ok (y/n/skip)?";
			  reply.ReadLine(cin);
	                }
		
			int ihole = col + row*ncol + 1000*ifoil;
			if( !reply.BeginsWith("n") ){ //Then this is a good hole. Add cuts to the list:
			  
			  xmean_hole[runnum][ihole] = gauss2D->GetParameter(3);
			  ymean_hole[runnum][ihole] = gauss2D->GetParameter(1);
			  sigx_hole[runnum][ihole] = gauss2D->GetParameter(4);
			  sigy_hole[runnum][ihole] = gauss2D->GetParameter(2);
			  use_hole[runnum][ihole] = 1;
			  nevents_hole[runnum][ihole] = 0;
			  
			  logfile << "Run " << runnum << " foil " << ifoil << ", sieve hole at (xhole,yhole)=(" << xhole[row][col] << ", " << yhole[row][col] << "), (xmean,ymean,sigmax,sigmay)=("
				  << xmean_hole[runnum][ihole] << ", " << ymean_hole[runnum][ihole] << ", " << sigx_hole[runnum][ihole] << ", " << sigy_hole[runnum][ihole] << ")" 
				  << endl;
			  
			} else {
			  use_hole[runnum][ihole] = 0;
			}
		      }
		      //gauss2D->ReleaseParameter(2);
		      //gauss2D->ReleaseParameter(4);
		    }
		  }
		}

		if(row==nrow-1 && col==ncol-1) {
		  htemp2D->Draw("col");
		  c1->Modified();
		  c1->Update();
		  for (int trow=0; trow<nrow; trow++){
			for  (int tcol=0; tcol<ncol; tcol++){
				int iholetemp = tcol + trow*ncol + 1000*ifoil;
				if (use_hole[runnum][iholetemp] != 0){
		  		 etemp[trow][tcol].Draw("SAME");
				}
				c1->Modified();
				c1->Update();
			}
		  }
		  cout << "Finished finding all sieve slit holes ... Accept and move on? " << endl;
		  TString accept;
		  accept.ReadLine(cin);
		  if (accept.BeginsWith("y")) {
			continue_getting_sieve_holes = false;
		  } else {
			auto_fitting = false;
		  }
		}

	      } // end loop over columns
	     } // end loop over rows
	   } // end while continue_getting_sieve_holes
	 } // end loop over foils
        } // end if sieve_slit_flag == 1
      } // end else determine sieve hole cuts
    } // end if

    //Now we've determined sieve hole cuts, next, we loop over all events in the tree again, this time setting up the fit for scattering angles:

    // chi^2 = sum_i=1^N (xptrue_i - sum_jklmn C_xp^jklmn xfp^j yfp^k xpfp^l ypfp^m xtar^n)^2
    // dchi^2 / dC_jklmn = -2 sum_i=1^N (xptrue_i - sum_jklmn C_xp^jklmn xfp^j yfp^k xpfp^l ypfp^m xtar^n ) * term
    
    cout << "Looping over events with target foil and sieve hole cuts, setting up SVD fit matrices..." << endl;

    //In this second loop, fill histograms of OLD (recon - true) vs true
    
    long nevents_run = 0;
    while( T->GetEntry( elist_temp->GetEntry( nevents_run++ ) ) ){ //second event loop to set up the fits:
      if( nevents_run %10000 == 0 ) cout << nevents_run << endl;

      trunnum = runnum;
      
      //First, select events by foil:
      txfp = T->hsxfp;
      tyfp = T->hsyfp;
      txpfp = T->hsxpfp;
      typfp = T->hsypfp;
      txbeam = x0beam - T->frx_cm;
      tybeam = y0beam - T->fry_cm;

      //reconstruct (xptar, yptar, ytar, xtar, zvertex) again using old recon coeffs:
      double ytarsum = 0.0;
      double yptarsum = 0.0;
      double xptarsum = 0.0;
      for( int term=0; term<num_recon_terms_old; term++ ){
	//In reconstruction, we want xtarget to be with respect to the spectrometer axis, so to get xtarget in spectrometer coordinates from x in global coordinates,
	//we SUBTRACT the x mispointing of the spectrometer:
	double lambda =
	  pow( txfp / 100.0, xfpexpon_old[term] ) *
	  pow( tyfp / 100.0, yfpexpon_old[term] ) *
	  pow( txpfp, xpfpexpon_old[term] ) *
	  pow( typfp, ypfpexpon_old[term] ) *
	  pow( (-tybeam - x_mispointing_HMS)/100.0, xtarexpon_old[term] );
	
	ytarsum += ytarcoeffs_old[term] * lambda;
	yptarsum += yptarcoeffs_old[term] * lambda;  
	xptarsum += xptarcoeffs_old[term] * lambda;
      }

      double ytarold = ytarsum * 100.0 + y_mispointing_HMS; //This ytarget is with respect to the origin; what is measured by the spectrometer is relative to HMS optical axis
      double yptarold = yptarsum + htheta_offset_old;
      double xptarold = xptarsum + hphi_offset_old;

      //The beam is assumed to have scattered from the point (xbeam + xpbeam * zfoil, ybeam + ypbeam*zfoil, zfoil)
      //This is the intersection of the horizontal component of the spectrometer ray with the nominal beam line.
      
      //ytarget = yvertex - yptar * zspec 
      // yvertex = zvertex * sin(theta) - xbeam * cos(theta)
      // zspec = zvertex * cos(theta) + xbeam * sin(theta)
      //ytarget = zvertex * sin(theta) - xbeam * cos(theta) - yptar*(zvertex*cos(theta) + xbeam*sin(theta))
      // zvertex*(sin(theta)-yptar*cos(theta)) = ytarget + xbeam*(cos(theta)+yptar*sin(theta))

      //A beam shift in the +x direction (beam right) shifts the z of the intersection point with the beam line by 
      double zvertexold = ( ytarold + txbeam * ( cos(thetaHMS) + yptarold * sin(thetaHMS) ) )/ ( sin(thetaHMS) - yptarold * cos(thetaHMS) );
      double zspecold = zvertexold * cos(thetaHMS) + txbeam * sin(thetaHMS);

      double xtarold = -tybeam - x_mispointing_HMS - zspecold * xptarold;
      
      for( int iter=0; iter<niter_xtarcorr; iter++ ){
	ytarsum = 0.0;
	yptarsum = 0.0;
	xptarsum = 0.0;
	
	for( int term=0; term<num_recon_terms_old; term++ ){
	  double lambda =
	    pow( txfp / 100.0, xfpexpon_old[term] ) *
	    pow( tyfp / 100.0, yfpexpon_old[term] ) *
	    pow( txpfp, xpfpexpon_old[term] ) *
	    pow( typfp, ypfpexpon_old[term] ) *
	    pow( xtarold/100.0, xtarexpon_old[term] ); 
	  
	  ytarsum += ytarcoeffs_old[term] * lambda;
	  yptarsum += yptarcoeffs_old[term] * lambda;
	  xptarsum += xptarcoeffs_old[term] * lambda;
	  
	}
	
	ytarold = ytarsum*100.0 + y_mispointing_HMS;
	yptarold = yptarsum + htheta_offset_old;
	xptarold = xptarsum + hphi_offset_old;
	
	zvertexold = ( ytarold + txbeam * ( cos(thetaHMS) + yptarold * sin(thetaHMS) ) )/ ( sin(thetaHMS) - yptarold * cos(thetaHMS) );
	zspecold = zvertexold * cos(thetaHMS) + txbeam*sin(thetaHMS);
	
	xtarold = -tybeam - x_mispointing_HMS - zspecold * xptarold;
      }
      
      xtarold += x_mispointing_HMS; //back to global coordinates for xtar.

      txptarrecon = xptarold;
      typtarrecon = yptarold;
      tytarrecon = ytarold;
      txtar = xtarold - x_mispointing_HMS;
      
      tifoil = -1;
      
      for( int ifoil=0; ifoil<nfoil[runnum]; ifoil++){ //figure out which foil: must pass zvertex and ytarget cuts:
	if( use_foil[runnum][ifoil] == 1 &&
	    zmin_foil[runnum][ifoil] <= zvertexold && zvertexold <= zmax_foil[runnum][ifoil] &&
	    ymin_foil[runnum][ifoil] <= ytarold && ytarold <= ymax_foil[runnum][ifoil] ){

	  tzfoil = zfoil[runnum][ifoil];
	  tifoil = ifoil;
	  break;
	}
      }
      
      if( tifoil >= 0 ){ //then this event passed the ytarget and z vertex selection for at least one of the foils:
	//Check whether or not this event is a sieve slit event:
	if( sieve_slit_flag[runnum] == 1 ){ //sieve slit in, determine "true" angles and vertex from foil position and sieve hole position:
	  //compute xy sieve:
	  double xsieve = xtarold + xptarold * z0sieve;
	  double ysieve = ytarold + yptarold * z0sieve;
	  
	  int ihole = -1, irow=-1, icol=-1;
	  for(int row=0; row<nrow; row++){
	    for(int col=0; col<ncol; col++){
	      int hole = col+row*ncol + 1000*tifoil;
	      
	      //Use 2.5sigma
	      if( use_hole[runnum][hole] == 1 &&
		  sqrt( pow( (xsieve-xmean_hole[runnum][hole])/sigx_hole[runnum][hole], 2 ) + pow( (ysieve-ymean_hole[runnum][hole])/sigy_hole[runnum][hole], 2 ) ) <= 2.5 ){
		ihole = col+row*ncol;
		irow = row;
		icol = col;
		break;
	      }
	    }
	  }

	  if( ihole >= 0 ){ //event passes sieve hole cut: compute "true" angles and ytarget:
	    nevents_hole[runnum][ihole+1000*tifoil]++;
	    
	    //The vertex coordinate is (xbeam,ybeam,zfoil) (in global coordinates):
	    double zspec_vertex_true = tzfoil*cos(thetaHMS) + txbeam*sin(thetaHMS);
	    double xspec_vertex_true = -tybeam;
	    double yspec_vertex_true = tzfoil*sin(thetaHMS) - txbeam*cos(thetaHMS);
	    txptartrue = (xhole[irow][icol] - xspec_vertex_true)/(z0sieve - zspec_vertex_true);
	    typtartrue = (yhole[irow][icol] - yspec_vertex_true)/(z0sieve - zspec_vertex_true);
	    tytartrue = yspec_vertex_true - typtartrue * zspec_vertex_true - y_mispointing_HMS; //Why subtract the spectrometer mispointing for ytarget?
	    txtartrue = xspec_vertex_true - txptartrue * zspec_vertex_true - x_mispointing_HMS;

	    
	    //Fill the root tree:
	    tout->Fill();
	    //above is in global coordinates! In spec coordinates, HMS points at x = x_mispointing_HMS, y = y_mispointing_HMS
	    //Since xptar = (xhole - xvertex)/(zsieve - zvertex), a positive shift in xsieve corresponds to a positive shift in xptartrue, and therefore a negative shift in (xptarrecon - xptartrue);
	    //in the xy sieve plot, the center sieve hole shows up at a negative xsieve value using the "old" recon matrix elements.
	    // assuming that xtar is handled correctly, this implies that a positive correction to xptar is needed. In that case, why is the zero offset of xptar negative?

	    //don't include this event in the fit if the max. number of events per hole/foil has been reached:
	    //This helps to weight each hole equally in the fit:
	    if( nevents_hole[runnum][ihole+1000*tifoil] < maxnperhole ){
	      double ytarsum_xtar = 0.0;
	      double xptarsum_xtar = 0.0;
	      double yptarsum_xtar = 0.0;
	      //compute sum of all xtar-dependent terms in the OLD recon matrix:
	      for( int term=0; term<num_recon_terms_old; term++ ){
		if( xtarexpon_old[term] > 0 ){
		  double lambda =
		    pow( txfp / 100.0, xfpexpon_old[term] ) *
		    pow( tyfp / 100.0, yfpexpon_old[term] ) *
		    pow( txpfp, xpfpexpon_old[term] ) *
		    pow( typfp, ypfpexpon_old[term] ) *
		    pow( txtartrue/100.0, xtarexpon_old[term] ); 
		  
		  ytarsum_xtar += ytarcoeffs_old[term] * lambda;
		  yptarsum_xtar += yptarcoeffs_old[term] * lambda;
		  xptarsum_xtar += xptarcoeffs_old[term] * lambda;
		}
	      }
	      
	      int iterm=0;
	      
	      vector<double> lambda_fp(nterms_fit);
	      
	      for( iterm=0; iterm<nterms_fit; iterm++ ){
		lambda_fp[iterm] =
		  pow( txfp/100.0, xfpexpon_new[iterm] ) *
		  pow( tyfp/100.0, yfpexpon_new[iterm] ) *
		  pow( txpfp, xpfpexpon_new[iterm] ) *
		  pow( typfp, ypfpexpon_new[iterm] ) *
		  pow( txtartrue/100.0, xtarexpon_new[iterm] );
	      }
	      
	      // //compute lambda for each reconstruction term:
	      // for(int ixfp=0; ixfp<=fitorder; ixfp++){
	      //   for(int iyfp=0; iyfp<=fitorder-ixfp; iyfp++){
	      // 	for(int ixpfp=0; ixpfp<=fitorder-ixfp-iyfp; ixpfp++){
	      // 	  for(int iypfp=0; iypfp<=fitorder-ixfp-iyfp-ixpfp; iypfp++){
	      // 	    lambda_fp[iterm++] =
	      // 	      pow( txfp / 100.0, ixfp ) *
	      // 	      pow( tyfp / 100.0, iyfp ) *
	      // 	      pow( txpfp, ixpfp ) *
	      // 	      pow( typfp, iypfp );
	      // 	  }
	      // 	}
	      //   }
	      // }
	      
	      //Now, how do we handle the case where we want to hold the xtar coeffs fixed?
	      for( iterm=0; iterm<nterms_fit; iterm++ ){
		for( int jterm=0; jterm<nterms_fit; jterm++ ){
		  //If fit_xtar_coeffs_flag is zero, don't include xtar-dependent coeffs in the fit:
		  if( xtarexpon_new[iterm] == 0 && xtarexpon_new[jterm] == 0 ){ //xtar-independent terms:
		    Fit_matrix_xptarget(iterm,jterm) += lambda_fp[iterm]*lambda_fp[jterm];
		    Fit_matrix_yptarget(iterm,jterm) += lambda_fp[iterm]*lambda_fp[jterm];
		    Fit_matrix_ytarget(iterm,jterm) += lambda_fp[iterm]*lambda_fp[jterm];
		  } else if( fit_xtar_coeffs_flag != 0 ){ //xtar-dependent term: behavior depends on fit_xtar_coeffs_flag:
		    Fit_matrix_xptarget(iterm,jterm) += lambda_fp[iterm]*lambda_fp[jterm];
		  }
		}
		if( xtarexpon_new[iterm] == 0 ){//xtar-independent term:
		  ytarget_vec(iterm) += lambda_fp[iterm] * (tytartrue/100.0 - ytarsum_xtar);
		  yptarget_vec(iterm) += lambda_fp[iterm] * (typtartrue - yptarsum_xtar );
		  if( fit_xtar_coeffs_flag == 0 ){ //increment xptar sum while subtracting xtar correction:
		    xptarget_vec(iterm) += lambda_fp[iterm] * (txptartrue - xptarsum_xtar );
		  } else { //xtar coeffs included in the fit, don't subtract xtar correction:
		    xptarget_vec(iterm) += lambda_fp[iterm] * txptartrue;
		  }
		} else if( fit_xtar_coeffs_flag != 0 ){ //xtar-dependent term: if included in the fit, increment xptarget sum:
		  xptarget_vec(iterm) += lambda_fp[iterm] * txptartrue;
		}
	      }
	    }
	  }
	} else { //no sieve slit, use this event for ytarget fit only, assuming yptarget is already known:
	  double zspec_vertex_true = tzfoil*cos(thetaHMS) + txbeam*sin(thetaHMS);
	  double xspec_vertex_true = -tybeam - x_mispointing_HMS;
	  double yspec_vertex_true = tzfoil*sin(thetaHMS) - txbeam*cos(thetaHMS) - y_mispointing_HMS;

	  txtartrue = xspec_vertex_true - xptarold * zspec_vertex_true;
	  typtartrue = yptarold;
	  txptartrue = xptarold;
	  tytartrue = yspec_vertex_true - yptarold * zspec_vertex_true;

	  //To include a foil from a non-sieve-slit run in the tree/fit, the angle reconstruction has to be "good":
	  //This means the projection of the track to the collimator has to lie within the physical dimensions of the collimator:
	  //pass octagon cut:
	  double xsieve = txtartrue + txptartrue * z0sieve;
	  double ysieve = tytartrue + typtartrue * z0sieve;

	  bool passed_octagon = false;

	  //corners of octagon are at:
	  // (x,y) = (-dy, dx/2), (-dy/2, dx), (dy/2, dx), (dy, dx/2),
	  //         (dy, -dx/2), (dy/2, -dx), (-dy/2,-dx), (-dy, -dx/2)
	  
	  if( fabs(ysieve - y0_coll) <= dy_coll && fabs( xsieve - x0_coll ) <= dx_coll &&
	      ysieve >= -dy_coll + (dy_coll/dx_coll)*(xsieve-dx_coll/2.0) &&
	      ysieve <= dy_coll - (dy_coll/dx_coll)*(xsieve-dx_coll/2.0) &&
	      xsieve >= -dx_coll + (dx_coll/dy_coll)*(ysieve-dy_coll/2.0) &&
	      xsieve <= dx_coll - (dx_coll/dy_coll)*(ysieve-dy_coll/2.0 ) ){ //passed octagon cut:
	    
	    tout->Fill();

	    if( nevents_foil[runnum][tifoil] < maxnperfoil ){
	      double ytarsum_xtar = 0.0;
	      //double xptarsum_xtar = 0.0;
	      //double yptarsum_xtar = 0.0;
	      //compute sum of all xtar-dependent terms in the OLD recon matrix:
	      for( int term=0; term<num_recon_terms_old; term++ ){
		if( xtarexpon_old[term] > 0 ){
		  double lambda =
		    pow( txfp / 100.0, xfpexpon_old[term] ) *
		    pow( tyfp / 100.0, yfpexpon_old[term] ) *
		    pow( txpfp, xpfpexpon_old[term] ) *
		    pow( typfp, ypfpexpon_old[term] ) *
		    pow( txtartrue/100.0, xtarexpon_old[term] ); 
		  
		  ytarsum_xtar += ytarcoeffs_old[term] * lambda;
		  //yptarsum_xtar += yptarcoeffs_old[term] * lambda;
		  //xptarsum_xtar += xptarcoeffs_old[term] * lambda;
		}
	      }
	      
	      int iterm=0;
	      
	      vector<double> lambda_fp(nterms_fit);
	      
	      //compute lambda for each reconstruction term:
	      // for(int ixfp=0; ixfp<=fitorder; ixfp++){
	      //   for(int iyfp=0; iyfp<=fitorder-ixfp; iyfp++){
	      //     for(int ixpfp=0; ixpfp<=fitorder-ixfp-iyfp; ixpfp++){
	      // 	for(int iypfp=0; iypfp<=fitorder-ixfp-iyfp-ixpfp; iypfp++){
	      for( iterm=0; iterm<nterms_fit; iterm++){
		lambda_fp[iterm] =
		  pow( txfp / 100.0, xfpexpon_new[iterm] ) *
		  pow( tyfp / 100.0, yfpexpon_new[iterm] ) *
		  pow( txpfp, xpfpexpon_new[iterm] ) *
		  pow( typfp, ypfpexpon_new[iterm] );
	      }
	      // 	}
	      //     }
	      //   }
	      // }
	      
	      for( iterm=0; iterm<nterms_fit; iterm++ ){
		if( xtarexpon_new[iterm] == 0 ){
		  for( int jterm=0; jterm<nterms_fit; jterm++ ){
		    if( xtarexpon_new[jterm] == 0 ){
		      Fit_matrix_ytarget(iterm,jterm) += lambda_fp[iterm]*lambda_fp[jterm];
		    }
		  }
		  ytarget_vec(iterm) += lambda_fp[iterm]*(tytartrue/100.0 - ytarsum_xtar);
		//xptarget_vec(iterm) += lambda_fp[iterm]*(txptartrue - xptarsum_xtar);
		//yptarget_vec(iterm) += lambda_fp[iterm]*(typtartrue - yptarsum_xtar);
		}
	      }
	    }
	  }
	}
      }
    }

    //We don't fill the output tree during this loop over runs...
    
    elist_temp->Delete();
    
    //Histograms for each run: xy at collimator,
    //y target at each foil

    ihisto++;
    ihistosieve += nfoil[runnum];
    
    delete T;
    C->Delete();
  }

  //Now, do some post-processing of these matrices depending on xtar behavior:
  for( int iterm=0; iterm<nterms_fit; iterm++ ){
    if( xtarexpon_new[iterm] != 0 ){
      ytarget_vec[iterm] = 0.0;
      yptarget_vec[iterm] = 0.0;
      if( fit_xtar_coeffs_flag == 0 ){
	xptarget_vec[iterm] = 0.0;
      }
    }
    for( int jterm=0; jterm<nterms_fit; jterm++ ){
      if( xtarexpon_new[jterm] != 0 ){
	if( jterm != iterm ){
	  Fit_matrix_ytarget(iterm,jterm) = 0.0;
	  Fit_matrix_yptarget(iterm,jterm) = 0.0;
	  if( fit_xtar_coeffs_flag == 0 ){
	    Fit_matrix_xptarget(iterm,jterm) = 0.0;
	  }
	} else {
	  Fit_matrix_ytarget(iterm,jterm) = 1.0;
	  Fit_matrix_yptarget(iterm,jterm) = 1.0;
	  if( fit_xtar_coeffs_flag == 0 ){
	    Fit_matrix_xptarget(iterm,jterm) = 1.0;
	  }
	}
      }
    }
  }

  for( int iterm=0; iterm<nterms_fit; iterm++ ){
    for( int jterm=0; jterm<nterms_fit; jterm++ ){
	cout << "i = " << iterm << " j = " << jterm << " Fit_matrix_ytarget = " << Fit_matrix_ytarget(iterm,jterm) << endl;
	cout << "    " << iterm << "     " << jterm << " Fit_matrix_xptarget = " << Fit_matrix_xptarget(iterm,jterm) << endl;
	cout << "    " << iterm << "     " << jterm << " Fit_matrix_yptarget = " << Fit_matrix_yptarget(iterm,jterm) << endl;
    }
  }
  
  for( int iterm=0; iterm<nterms_fit; iterm++ ){
	cout << "i = " << iterm << "ytarget_vec = " << ytarget_vec[iterm] << endl;
	cout << "    " << iterm << "xptarget_vec = " << xptarget_vec[iterm] << endl;
	cout << "    " << iterm << "yptarget_vec = " << yptarget_vec[iterm] << endl;
  }
  
  TDecompSVD A_ytarget(Fit_matrix_ytarget);
  TDecompSVD A_yptarget(Fit_matrix_yptarget);
  TDecompSVD A_xptarget(Fit_matrix_xptarget);

  cout << "Fitting ytarget...";
  bool success_ytarget = A_ytarget.Solve( ytarget_vec );

  cout << " done, success = " << success_ytarget << endl;
  cout << "Fitting xptarget...";
  bool success_xptarget = A_xptarget.Solve( xptarget_vec );
  cout << " done, success = " << success_xptarget << endl;

  cout << "Fitting yptarget...";
  bool success_yptarget = A_yptarget.Solve( yptarget_vec );
  cout << " done, success = " << success_yptarget << endl;
  
  
  //At the end of the loop over runs, write out the new optics coefficients:
  for( int order=0; order<=fitorder; order++){
    for( int term=0; term<nterms_fit; term++){
      int degree = xtarexpon_new[term] + xfpexpon_new[term] + xpfpexpon_new[term] + yfpexpon_new[term] + ypfpexpon_new[term];

      if( degree == order ){ //output terms in order of increasing degree:
      
	if( xtarexpon_new[term] != 0 ){ //Modify solution vectors:
	  int old_term = -1;
	  for( int termold=0; termold<num_recon_terms_old; termold++ ){
	    if( xfpexpon_old[termold] == xfpexpon_new[term] &&
		yfpexpon_old[termold] == yfpexpon_new[term] &&
		xpfpexpon_old[termold] == xpfpexpon_new[term] &&
		ypfpexpon_old[termold] == ypfpexpon_new[term] &&
		xtarexpon_old[termold] == xtarexpon_new[term] ){
	      old_term = termold;
	      break;
	    }
	  }
	  if( old_term >= 0 ){
	    ytarget_vec(term) = ytarcoeffs_old[old_term];
	    yptarget_vec(term) = yptarcoeffs_old[old_term];
	    if( fit_xtar_coeffs_flag == 0 ){
	      xptarget_vec(term) = xptarcoeffs_old[old_term];
	    }
	  }
	}
    
	TString outputline;
	outputline.Form( " %16.9g%16.9g%16.9g%16.9g %d%d%d%d%d", xptarget_vec[term], ytarget_vec[term], yptarget_vec[term], deltacoeffs_new[term],
			     xfpexpon_new[term], xpfpexpon_new[term], yfpexpon_new[term], ypfpexpon_new[term], xtarexpon_new[term] );
	newcoeffsfile << outputline.Data() << endl;
      }
    }
  }
  // for( int term=0; term<num_recon_terms_old; term++ ){
  //   if(xtarexpon_old[term] > 0 ){
  //     TString outputline;
  //     outputline.Form( " %16.9g%16.9g%16.9g%16.9g %d%d%d%d%d", xptarcoeffs_old[term], ytarcoeffs_old[term], yptarcoeffs_old[term], deltacoeffs_old[term],
  // 		       xfpexpon_old[term], xpfpexpon_old[term], yfpexpon_old[term], ypfpexpon_old[term], xtarexpon_old[term] );
  //     newcoeffsfile << outputline.Data() << endl;
  //   }
  // }

  newcoeffsfile << endline.Data() << endl;    
  
  fout->cd();
  ytarget_histos->Write();
  zbeam_histos->Write();
  xysieve_histos->Write();

  fout->Write();
  
  fout->Close();
      
}
