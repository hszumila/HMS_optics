#ifndef myMath_h
#define myMath_h 1

#include <vector>


class TH1D;


class MultiPeakFunc {
  public:
    MultiPeakFunc(int numPeaks);
    ~MultiPeakFunc();

    double Evaluate(double* x, double* pars);

    double getNPeaks();

  private:
    int nPeaks;
};


class Peak {
  public:
    Peak();
    Peak(double norm, double mean, double sigma);
    ~Peak();

    double norm;
    double mean;
    double sigma;
};


Peak fitPeak(TH1D* histo, double normInit, double meanInit, double sigmaInit);

std::vector<Peak> fitMultiPeak(TH1D* histo, double sigma=0.5);

std::size_t getClosestIndex(double value, std::vector<double>& reference);


#endif  // myMath_h
