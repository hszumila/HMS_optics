#include "myMath.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "TMath.h"
#include "TH1D.h"
#include "TSpectrum.h"
#include "TF1.h"

bool compare(Peak p1,Peak p2){return (p1.mean<p2.mean);}

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
    histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax()
  );
  fitFunc->SetNpx(1000);

  fitFunc->SetParameter(0, normInit);
  fitFunc->SetParameter(1, meanInit);
  fitFunc->SetParameter(2, sigmaInit);
  histo->Fit("fitFunc", "Q");

  Peak peak(
    fitFunc->GetParameter(0),
    fitFunc->GetParameter(1),
    TMath::Abs(fitFunc->GetParameter(2))
  );

  return peak;
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
