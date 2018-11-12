
#include "myMath.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "TMath.h"
#include "TH1D.h"
#include "TSpectrum.h"
#include "TF1.h"

bool compare(Peak p1,Peak p2){return (p1.mean<p2.mean);}

bool compareHt(Peak p1,Peak p2){
  /*
  TF1 *f1 = new TF1("f1","[0]*exp(-0.5*((x-[1])/[2])**2) ",-15,15);
  TF1 *f2 = new TF1("f2","[0]*exp(-0.5*((x-[1])/[2])**2) ",-15,15);
  f1->FixParameter(0,p1.norm);
  f1->FixParameter(1,p1.mean);
  f1->FixParameter(2,p1.sigma);
  f2->FixParameter(0,p2.norm);
  f2->FixParameter(1,p2.mean);
  f2->FixParameter(2,p2.sigma);
  double max1 = f1->GetMaximum();
  double max2 = f2->GetMaximum();

  return (max1 > max2);
  */
  return (p1.norm>p2.norm);
}


// MultiPeakFunc.
MultiPeakFunc::MultiPeakFunc(int numPeaks) : nPeaks(numPeaks) {}


MultiPeakFunc::~MultiPeakFunc() {}


double MultiPeakFunc::Evaluate(double* x, double* pars) {
  double result = 0;
  for (int iPeak=0; iPeak<nPeaks; iPeak++) {
    double norm = pars[3*iPeak];
    double mean = pars[3*iPeak+1];
    double sigma = pars[3*iPeak+2];

    result += norm*TMath::Gaus(x[0], mean, sigma);
  }

  return result;
}


// Peak.
Peak::Peak() : norm(0.0), mean(0.0), sigma(0.0) {}


Peak::Peak(double norm, double mean, double sigma) :
  norm(norm), mean(mean), sigma(sigma)
{}


Peak::~Peak() {}


// Fitting.
Peak fitPeak(TH1D* histo, double normInit, double meanInit, double sigmaInit) {
  TF1* fitFunc = new TF1(
    "fitFunc",
    "gaus(0)",
    meanInit-2*TMath::Abs(sigmaInit),meanInit+2*TMath::Abs(sigmaInit)
    //histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax()
  );
  fitFunc->SetNpx(1000);

  fitFunc->SetParameter(0, normInit);
  fitFunc->SetParameter(1, meanInit);
  fitFunc->SetParameter(2, sigmaInit);
  histo->Fit("fitFunc", "QR");

  double aa = fitFunc->GetParameter(0);
  double mu = fitFunc->GetParameter(1);
  double ssig = TMath::Abs(fitFunc->GetParameter(2));
  double chi2 = fitFunc->GetChisquare();
  //std::cout<<"initial:\t"<<normInit<<"\t"<<meanInit<<"\t"<<sigmaInit<<std::endl;
  //std::cout<<"\tfitted:\t"<<aa<<"\t"<<mu<<"\t"<<ssig<<"\t"<<chi2<<std::endl;

  if ((mu>meanInit+0.2 || mu<meanInit-0.2) && (chi2>20 || ssig>1)){
    mu = meanInit;
  }
 
 
   if (aa>normInit*20 || chi2<0.1  ){//|| aa<normInit/10.0){
    ssig = 0;
  }
   if (ssig>sigmaInit){
    ssig = TMath::Abs(sigmaInit);
  }

   if (mu-meanInit>0.5){mu = meanInit;}
  //compare number of events in central peak to norm
   // std::cout<<"\tcorrected:\t"<<aa<<"\t"<<mu<<"\t"<<ssig<<std::endl;

  // if (ssig==0){std::cout<<"flagged event for removal!!!!"<<std::endl;}
  Peak peak(
    aa,
    mu,
    ssig
  );

  return peak;
}
/*
std::vector<Peak> sortByHeight(std::vector<Peak> peaksFound, int nFoil) {
  std::sort(peaksFound.begin(), peaksFound.end(),compareHt);

  std::vector<Peak> correctPeaks;
  for (int ii=0; ii<peaksFound.size(); ii++){
    //std::cout<<"peak after sort: "<<peaksFound.at(ii).norm<<","<<peaksFound.at(ii).mean<<","<<peaksFound.at(ii).sigma<<std::endl;
    if (peaksFound.at(ii).sigma>0.1){
      correctPeaks.push_back(peaksFound.at(ii));
    }
    if(correctPeaks.size()==nFoil) break;
  }

  std::sort(correctPeaks.begin(), correctPeaks.end(),compare);

  return correctPeaks;

}
*/

std::vector<Peak> findPeaks(TH1D* histo, int nfoil) {
  int nbins = histo->GetNbinsX();
  int range = nbins/nfoil;
  std::vector<Peak> peaksFound;

  double sigma = (histo->GetBinCenter(range) - histo->GetBinCenter(1))/6.0;
  // std::cout<<"nbins: "<<nbins<<" range: "<<range<<" sigma: "<<sigma<<std::endl;
  double prevMean = -100.0;
  double delta = 2.0*sigma;//0.3;

  int ll=1;
  for (uint ii=0; ii<nfoil; ii++){
    int maxrange = ll+range;
    if (ll+range>nbins-1){
      maxrange = nbins-2;
    }
    histo->GetXaxis()->SetRange(ll,maxrange);
    int maxBin = histo->GetMaximumBin();
    double mean = histo->GetBinCenter(maxBin);  
    double maxContent = histo->GetBinContent(maxBin);

    //double srange =  histo->GetBinCenter(ll);
    //double erange =  histo->GetBinCenter(maxrange);
    // std::cout<<"range: "<<srange<<" end range: "<<erange<<"peak at: "<<mean<<" start bin: "<<ll<<" max bin: "<<maxrange<<" maxbin: "<<maxBin<<std::endl;
    /*
    if (abs(prevMean-mean)<2.0*sigma){
      histo->GetXaxis()->SetRange(maxBin+1,maxrange);
      int maxBin2 = histo->GetMaximumBin();
      mean = histo->GetBinCenter(maxBin2);
      maxContent = histo->GetBinContent(maxBin);
    }
    */
    //histo->GetXaxis()->SetRange(maxBin,maxrange);
    //double sigma = (mean - histo->GetBinCenter(histo->GetMinimumBin()))/3.0;
    //  std::cout<<"here's what's happening: "<<maxBin<<std::endl;
    peaksFound.push_back(
			 Peak(
			      maxContent,
			      mean,
			      sigma
			      )
			 );
    prevMean = mean;
    ll+=range;

  }
  
  //1st check: if any two peaks are the same, we need to shift the window
  int flag = 0;
  for (int jj=1; jj<peaksFound.size(); jj++){
    if (abs(peaksFound.at(jj).mean - peaksFound.at(jj-1).mean)<delta){
      double mean0 = peaksFound.at(jj).mean;
      double mean1 = peaksFound.at(jj-1).mean;
         std::cout<<"too close together: "<<mean0<<" and "<<mean1<<std::endl;
      peaksFound.clear();
      flag ++;
      break;
    }
  }
  if (flag>0){
    //  std::cout<<"shifting the window...."<<std::endl;
    int kk= range/2;//1
    for (uint ii=0; ii<nfoil; ii++){
      int maxrange = kk+range;
      if (ii==0){
	maxrange = kk+range;///2.0;
      }
      if (maxrange>nbins-1){
	maxrange = nbins-2;
      }    
      histo->GetXaxis()->SetRange(kk,maxrange);
      double srange =  histo->GetBinCenter(kk);
      double erange =  histo->GetBinCenter(maxrange);
      //  std::cout<<"foil: "<<ii<<" kk: "<<kk<<std::endl;
      //   std::cout<<"range: "<<srange<<" end range: "<<erange<<std::endl;


      int maxBin = histo->GetMaximumBin();
      double binc = histo->GetBinCenter(maxBin);
      //  std::cout<<"\tmaxbin: "<<maxBin<<" binc: "<<binc<<std::endl;

      double mean = histo->GetBinCenter(maxBin);  
      double maxContent = histo->GetBinContent(maxBin);
      //histo->GetXaxis()->SetRange(maxBin,maxrange);
      //double sigma = (mean - histo->GetBinCenter(histo->GetMinimumBin()))/3.0;
      peaksFound.push_back(
			   Peak(
				maxContent,
				mean,
				sigma
				)
			   );
      //if (ii==0){
      //kk+=range/2.0;
      //}
      //else{
	kk+=range;
	//}
      //     std::cout<<"foil: "<<ii<<" kk: "<<kk<<std::endl;



    }
  }
    /*
    //2nd check: if any two peaks are the same, we need to shift the window
  int flag = 0;
  for (int jj=1; jj<peaksFound.size(); jj++){
    if (abs(peaksFound.at(jj).mean - peaksFound.at(jj-1).mean)<delta){
      double mean0 = peaksFound.at(jj).mean;
      double mean1 = peaksFound.at(jj-1).mean;
      std::cout<<"too close together: "<<mean0<<" and "<<mean1<<std::endl;
      peaksFound.clear();
      flag ++;
      break;
    }
  }
  if (flag>0){
    std::cout<<"shifting the window...."<<std::endl;
    int kk= 1;
    for (int ii=0; ii<nfoil; ii++){
      int maxrange = kk+range;
      if (ii==0){
	maxrange = kk+range/2.0;
      }
      if (maxrange>nbins-1){
	maxrange = nbins-2;
      }    
      histo->GetXaxis()->SetRange(kk,maxrange);
      double srange =  histo->GetBinCenter(kk);
      double erange =  histo->GetBinCenter(maxrange);
      std::cout<<"foil: "<<ii<<" kk: "<<kk<<std::endl;
      std::cout<<"range: "<<srange<<" end range: "<<erange<<std::endl;


      int maxBin = histo->GetMaximumBin();
      double binc = histo->GetBinCenter(maxBin);
      std::cout<<"\tmaxbin: "<<maxBin<<" binc: "<<binc<<std::endl;

      peaksFound.push_back(
			   Peak(
				histo->GetBinContent(maxBin),
				histo->GetBinCenter(maxBin),
				sigma
				)
			   );
      if (ii==0){
	kk+=range/2.0;
      }
      else{
	kk+=range;
      }




  }
    */
  
  return peaksFound;
}

std::vector<Peak> fitMultiPeak(TH1D* histo, double sigma) {
  TSpectrum* spec = new TSpectrum();
  int nPeaks = spec->Search(histo, sigma, "goff");

  double* peaksX = (double*)spec->GetPositionX();
  double* peaksY = (double*)spec->GetPositionY();

  MultiPeakFunc* peaksFunc = new MultiPeakFunc(nPeaks);
  TF1* fitFunc = new TF1(
    "fitFunc",
    peaksFunc, &MultiPeakFunc::Evaluate,
    histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax(),
    3*nPeaks
  );
  fitFunc->SetNpx(1000);

  double pars[100];
  for (int iPeak=0; iPeak<nPeaks; ++iPeak) {
    pars[3*iPeak] = peaksY[iPeak];
    pars[3*iPeak+1] = peaksX[iPeak];
    pars[3*iPeak+2] = sigma;
  }
  fitFunc->SetParameters(pars);
  histo->Fit("fitFunc", "Q");

  std::vector<Peak> peaks;
  for (int iPeak=0; iPeak<nPeaks; ++iPeak) { 
    peaks.push_back(
		    Peak(
			 fitFunc->GetParameter(3*iPeak),
			 fitFunc->GetParameter(3*iPeak+1),
			 TMath::Abs(fitFunc->GetParameter(3*iPeak+2))
			 )
		    );
  }
  

  //  std::sort(
  //	    peaks.begin(), peaks.end(),[](Peak p1, Peak p2) {
  //	      return p1.mean<p2.mean;
  //	    });

  std::sort(peaks.begin(), peaks.end(),compare);
  
  delete spec;

  return peaks;
}

std::vector<Peak> selectMultiPeakZ(TH1D* histo, int nfoil, double sinTheta) {
  std::vector<Peak> peaks;
  if (nfoil<3){//optics 2
    peaks.push_back(Peak(100,-5.0,1.5));
    peaks.push_back(Peak(100,5.0,1.5));
  }
  else {
    peaks.push_back(Peak(100,-12.5,1));
    peaks.push_back(Peak(100,-3.2,1));
    peaks.push_back(Peak(100,6.5,1));
  }
  return peaks;
}

std::vector<Peak> selectMultiPeakY(TH1D* histo, int nfoil, double sinTheta) {
  std::vector<Peak> peaks;
  if (nfoil<3){//optics 2
    peaks.push_back(Peak(100,-5.0*sinTheta,1.5));
    peaks.push_back(Peak(100,5.0*sinTheta,1.5));
  }
  else {
    peaks.push_back(Peak(100,-3.2,0.5));
    peaks.push_back(Peak(100,-1,0.5));
    peaks.push_back(Peak(100,1.5,0.5));
  }
  return peaks;
}

std::size_t getClosestIndex(double value, std::vector<double>& reference) {
  size_t index;
  double min_dist = INFINITY;
  double distance;

  for (size_t i=0; i<reference.size(); ++i) {
    distance = std::abs(value - reference.at(i));
    if (distance < min_dist) {
      index = i;
      min_dist = distance;
    }
  }

  return index;
}
