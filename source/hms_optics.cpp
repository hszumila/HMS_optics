// Standard includes.
#include <chrono>
#include <cmath>
#include <iostream>
  using std::cin;
  using std::cout;
  using std::endl;
#include <thread>
#include <utility>
#include <vector>

// ROOT includes.
#include "TCanvas.h"
#include "TDecompSVD.h"
#include "TDirectory.h"
#include "TEllipse.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMarker.h"
#include "TMatrixD.h"
#include "TRint.h"
#include "TString.h"
#include "TSystem.h"
#include "TVectorD.h"

// Project includes.
#include "cmdOptions.hpp"
#include "myConfig.hpp"
#include "myEvent.hpp"
#include "myMath.hpp"
#include "myOther.hpp"
#include "myRecMatrix.hpp"


int hms_optics(const cmdOptions::OptionParser_hmsOptics& cmdOpts);


int main(int argc, char* argv[]) {
  // Parse command line options for hms_optics.
  cmdOptions::OptionParser_hmsOptics cmdOpts;
  try {
    cmdOpts.init(argc, argv);
  }
  catch (const std::runtime_error& err) {
    cout << "hms_optics: " << err.what() << endl;
    cout << "hms_optics: Try `tmp -h` for more information." << endl;
    return 1;
  }
  if (cmdOpts.displayHelp) {
    cmdOpts.printHelp();
    return 0;
  }

  // Create command line options for ROOT.
  int argcRoot = 3;
  static char argvRoot[][100] = {"-q", "-l"};
  static char* argvRootList[] = {argv[0], argvRoot[0], argvRoot[1], NULL};

  // Run hms_optics as an application.
  TRint *theApp = new TRint("app", &argcRoot, argvRootList);
  int retCode = hms_optics(cmdOpts);
  theApp->Run(kTRUE);

  // Cleanup and exit.
  delete theApp;

  return retCode;
}


int hms_optics(const cmdOptions::OptionParser_hmsOptics& cmdOpts) {
  cout << "Reading config file: `" << cmdOpts.configFileName << "`." << endl;
  config::Config conf = config::loadConfigFile(cmdOpts.configFileName);

  cout << "Reading matrix file: `" << conf.recMatrixFileNameOld << "`." << endl;
  RecMatrix recMatrix = readMatrixFile(conf.recMatrixFileNameOld);


  cout << "Initializing new matrix." << endl;
  // Copy header and delta elements from old matrix.
  // Initialize other elements to 0.
  RecMatrix recMatrixNew;
  recMatrixNew.header = recMatrix.header;
  double C_D;
  std::vector<RecMatrixLine>::iterator lineIt;
  // Construct order by order.
  for (int order=0; order<=conf.fitOrder; ++order) {
    for (int m=0; m<=order; ++m) {
      for (int l=0; l<=order-m; ++l) {
        for (int k=0; k<=order-m-l; ++k) {
          for (int j=0; j<=order-m-l-k; ++j) {
            for (int i=0; i<=order-m-l-k-j; ++i) {
              if (i+j+k+l+m != order) continue;

              lineIt = recMatrix.findLine(i, j, k, l, m);
              if (lineIt != recMatrix.end()) C_D = lineIt->C_D;
              else C_D = 0.0;

              recMatrixNew.addLine(
                0.0, 0.0, 0.0, C_D,
                i, j, k, l, m
              );
            }
          }
        }
      }
    }
  }
  int recMatrixNewLen = static_cast<int>(recMatrixNew.size());
  cout << "  " << recMatrixNewLen << " terms" << endl;


  // Prepare for analysis.
  char tmp;

  TFile fo(cmdOpts.rootFileName.c_str(), "RECREATE");
  TDirectory* dir;

  std::vector<double> xSievePhys = conf.getSieveHolesX();
  std::vector<double> ySievePhys = conf.getSieveHolesY();

  TMatrixD xpTarFitMat(recMatrixNewLen, recMatrixNewLen);
  TMatrixD yTarFitMat(recMatrixNewLen, recMatrixNewLen);
  TMatrixD ypTarFitMat(recMatrixNewLen, recMatrixNewLen);
  TVectorD xpTarFitVec(recMatrixNewLen);
  TVectorD yTarFitVec(recMatrixNewLen);
  TVectorD ypTarFitVec(recMatrixNewLen);

  TCanvas* c1 = new TCanvas("c1", "c1", 100, 100, 600, 400);
  TCanvas* c2 = new TCanvas("c2", "c2", 100, 540, 600, 400);
  TCanvas* c3 = new TCanvas("c3", "c3", 702, 100, 600, 400);
  gPad->Update();


  cout << "Reading and analyzing root files:" << endl;
  for (const auto& runConf : conf.runConfigs) {  // run loop
    cout << "  " << runConf.runNumber << ":" << endl;

    const double& cosTheta = runConf.HMS.cosTheta;
    const double& sinTheta = runConf.HMS.sinTheta;
    const size_t nFoils = runConf.zFoils.size();

    // Create directory in output ROOT file for histograms.
    dir = fo.mkdir(
      TString::Format("run_%d", runConf.runNumber),
      TString::Format("histograms for run %d", runConf.runNumber)
    );
    dir->cd();


    // Reading events from input ROOT files.
    std::vector<Event> events = readEvents(runConf);
    size_t nEvents = events.size();
    cout << "    " << nEvents << " events survived cuts." << endl;


    cout << "    Reconstructing events: ";
    size_t iEvent = 0;
    double xpSum;
    double ySum;
    double ypSum;
    double lambda;

    reportProgressInit();
    for (auto& event : events) {  // reconstruction event loop
      if (iEvent%2000 == 0) reportProgress(iEvent, nEvents);

      xpSum = 0.0;
      ySum = 0.0;
      ypSum = 0.0;
      lambda = 0.0;

      event.xTar = -event.yVer - runConf.HMS.xMispointing;

      for (const auto& line : recMatrix.matrix) {  // line loop
        lambda =
          pow(event.xFp/100.0, line.E_x) *
          pow(event.xpFp, line.E_xp) *
          pow(event.yFp/100.0, line.E_y) *
          pow(event.ypFp, line.E_yp) *
          pow(event.xTar/100.0, line.E_xTar);

        xpSum += line.C_Xp * lambda;
        ySum += line.C_Y * lambda;
        ypSum += line.C_Yp * lambda;
      }   // line loop

      event.xpTar = xpSum + runConf.HMS.phiOffset;
      event.yTar = ySum*100.0 + runConf.HMS.yMispointing;
      event.ypTar = ypSum + runConf.HMS.thetaOffset;

      event.zVer =
        (event.yTar - event.xVer*(cosTheta + event.ypTar*sinTheta)) /
        (sinTheta - event.ypTar*cosTheta);

      event.xTarVer = -event.yVer;
      event.yTarVer = event.zVer*sinTheta + event.xVer*cosTheta;
      event.zTarVer = event.zVer*cosTheta - event.xVer*sinTheta;

      event.xTar = event.xTarVer - event.zTarVer*event.xpTar - runConf.HMS.xMispointing;

      // Iterate with updated value of xTar.
      for (int iIter=0; iIter<conf.xTarCorrIterNum; ++iIter) {  // iteration loop
        xpSum = 0.0;
        ySum = 0.0;
        ypSum = 0.0;
        lambda = 0.0;

        for (const auto& line : recMatrix.matrix) {  // line loop
          lambda =
            pow(event.xFp/100.0, line.E_x) *
            pow(event.xpFp, line.E_xp) *
            pow(event.yFp/100.0, line.E_y) *
            pow(event.ypFp, line.E_yp) *
            pow(event.xTar/100.0, line.E_xTar);

          xpSum += line.C_Xp * lambda;
          ySum += line.C_Y * lambda;
          ypSum += line.C_Yp * lambda;
        }  // line loop

        event.xpTar = xpSum + runConf.HMS.phiOffset;
        event.yTar = ySum*100.0 + runConf.HMS.yMispointing;
        event.ypTar = ypSum + runConf.HMS.thetaOffset;

        event.zVer =
          (event.yTar - event.xVer*(cosTheta + event.ypTar*sinTheta)) /
          (sinTheta - event.ypTar*cosTheta);

        event.xTarVer = -event.yVer;
        event.yTarVer = event.zVer*sinTheta + event.xVer*cosTheta;
        event.zTarVer = event.zVer*cosTheta - event.xVer*sinTheta;

        event.xTar = event.xTarVer - event.zTarVer*event.xpTar - runConf.HMS.xMispointing;
      }  // iteration loop

      event.xTar += runConf.HMS.xMispointing;

      event.xSieve = event.xTar + event.xpTar*conf.sieve.z0;
      event.ySieve = event.yTar + event.ypTar*conf.sieve.z0;

      ++iEvent;
    }  // reconstruction event loop
    reportProgressFinish();


    cout << "    Fitting target foils." << endl;

    // Setting historgams.
    double minx = runConf.zFoils.front() - 3.0;
    double maxx = runConf.zFoils.back() + 3.0;
    int binsx = 10 * static_cast<int>(maxx-minx);
    TH1D zVerHist(
      TString::Format("zVer"),
      TString::Format("zVer for run %d", runConf.runNumber),
      binsx, minx, maxx
    );
    zVerHist.GetXaxis()->SetTitle("z_{vertex}  [cm]");

    minx *= sinTheta;
    maxx *= sinTheta;
    TH1D yTarHist(
      TString::Format("yTar"),
      TString::Format("yTar for run %d", runConf.runNumber),
      binsx, minx, maxx
    );
    yTarHist.GetXaxis()->SetTitle("y_{target}  [cm]");

    // Filling the histograms.
    for (const auto& event : events) {
      zVerHist.Fill(event.zVer);
      yTarHist.Fill(event.yTar);
    }

    // Fitting the histograms.
    std::vector<Peak> zVerPeaks = fitMultiPeak(&zVerHist, 0.2);
    std::vector<Peak> yTarPeaks = fitMultiPeak(&yTarHist, 0.2);

    // Plotting the histograms.
    c1->cd();
    zVerHist.Draw();
    c1->Update();
    double miny = gPad->GetUymin();
    double maxy = gPad->GetUymax();
    // Add lines for physical positions of foils.
    std::vector<TLine> zFoilLines(nFoils);
    for (size_t iFoil=0; iFoil<nFoils; ++iFoil) {
      zFoilLines.at(iFoil) = TLine(
        runConf.zFoils.at(iFoil), miny,
        runConf.zFoils.at(iFoil), maxy
      );
      zFoilLines.at(iFoil).SetLineColor(6);
      zFoilLines.at(iFoil).SetLineWidth(2);
      zVerHist.GetListOfFunctions()->Add(&(zFoilLines.at(iFoil)));
    }
    c1->Update();
    gPad->Update();
    zVerHist.Write();

    c3->cd();
    yTarHist.Draw();
    c3->Update();
    miny = gPad->GetUymin();
    maxy = gPad->GetUymax();
    std::vector<TLine> yTarLines(nFoils);
    for (size_t iFoil=0; iFoil<nFoils; ++iFoil) {
      yTarLines.at(iFoil) = TLine(
        runConf.zFoils.at(iFoil)*sinTheta, miny,
        runConf.zFoils.at(iFoil)*sinTheta, maxy
      );
      yTarLines.at(iFoil).SetLineColor(6);
      yTarLines.at(iFoil).SetLineWidth(2);
      yTarHist.GetListOfFunctions()->Add(&(yTarLines.at(iFoil)));
    }
    c3->Update();
    gPad->Update();
    yTarHist.Write();

    if (cmdOpts.automatic) {
      std::this_thread::sleep_for(std::chrono::milliseconds(cmdOpts.delay));
    }
    else {
      cout << "    Continue? ";
      cin >> tmp;
    }

    c1->Clear();
    gPad->Update();
    c3->Clear();
    gPad->Update();


    cout << "    Fitting sieve holes." << endl;

    std::vector<TH2D> xySieveHists(nFoils);
    std::vector<TLine> xSieveLines(conf.sieve.nRow);
    std::vector<TLine> ySieveLines(conf.sieve.nCol);

    minx = xSievePhys.front() - 0.1*(xSievePhys.back()-xSievePhys.front());
    maxx = xSievePhys.back() + 0.1*(xSievePhys.back()-xSievePhys.front());
    binsx = 10 * static_cast<int>(maxx-minx);
    miny = ySievePhys.front() - 0.1*(ySievePhys.back()-ySievePhys.front());
    maxy = ySievePhys.back() + 0.1*(ySievePhys.back()-ySievePhys.front());
    int binsy = 10 * static_cast<int>(maxy-miny);

    // Construct lines for physical positions of sieve holes.
    for (size_t iRow=0; iRow<conf.sieve.nRow; ++iRow) {
      xSieveLines.at(iRow) = TLine(
        xSievePhys.at(iRow), miny,
        xSievePhys.at(iRow), maxy
      );
      xSieveLines.at(iRow).SetLineColor(6);
      xSieveLines.at(iRow).SetLineWidth(2);
    }
    for (size_t iCol=0; iCol<conf.sieve.nRow; ++iCol) {
      ySieveLines.at(iCol) = TLine(
        minx, ySievePhys.at(iCol),
        maxx, ySievePhys.at(iCol)
      );
      ySieveLines.at(iCol).SetLineColor(6);
      ySieveLines.at(iCol).SetLineWidth(2);
    }

    // Setting histograms.
    for (size_t iCol=0; iCol<conf.sieve.nCol; ++iCol) {
      ySieveLines.at(iCol) = TLine(
        minx, ySievePhys.at(iCol),
        maxx, ySievePhys.at(iCol)
      );
      ySieveLines.at(iCol).SetLineColor(6);
      ySieveLines.at(iCol).SetLineWidth(2);
    }

    for (size_t iFoil=0; iFoil<nFoils; ++iFoil) {
      xySieveHists.at(iFoil) = TH2D(
        TString::Format("xySieve_%d", static_cast<int>(iFoil)),
        TString::Format("xySieve for foil %d run %d", static_cast<int>(iFoil), runConf.runNumber),
        binsx, minx, maxx,
        binsy, miny, maxy
      );
      xySieveHists.at(iFoil).GetXaxis()->SetTitle("x_{sieve}  [cm]");
      xySieveHists.at(iFoil).GetYaxis()->SetTitle("y_{sieve}  [cm]");
      for (auto& line : xSieveLines) {
        xySieveHists.at(iFoil).GetListOfFunctions()->Add(&line);
      }
      for (auto& line : ySieveLines) {
        xySieveHists.at(iFoil).GetListOfFunctions()->Add(&line);
      }
    }

    // Filling the histograms.
    for (const auto& event : events) {
      for (size_t iFoil=0; iFoil<nFoils; ++iFoil) {
        if (
          zVerPeaks.at(iFoil).mean - 3.0*zVerPeaks.at(iFoil).sigma <= event.zVer &&
          event.zVer <= zVerPeaks.at(iFoil).mean + 3.0*zVerPeaks.at(iFoil).sigma &&
          yTarPeaks.at(iFoil).mean - 3.0*yTarPeaks.at(iFoil).sigma <= event.yTar &&
          event.yTar <= yTarPeaks.at(iFoil).mean + 3.0*yTarPeaks.at(iFoil).sigma
        ) {
          xySieveHists.at(iFoil).Fill(event.xSieve, event.ySieve);
          break;
        }
      }
    }

    // Setting things before starting.
    std::vector<std::vector<Peak>> xSievePeakss(nFoils);
    std::vector<std::vector<Peak>> ySievePeakss(nFoils);
    std::vector<std::vector<std::size_t>> xSieveIndexess(nFoils);
    std::vector<std::vector<std::size_t>> ySieveIndexess(nFoils);
    std::vector<std::vector<std::size_t>> nEventss(nFoils);
    std::vector<std::vector<TEllipse>> ellipsess(nFoils);

    TH1D* tmpHist;
    TMarker* tmpMark = new TMarker(0.0, 0.0, 22);
    tmpMark->SetMarkerColor(2);

    // Fit sieve holes for each foil.
    for (size_t iFoil=0; iFoil<nFoils; ++iFoil) {  // foil loop
      cout << "      Foil " << iFoil << "." << endl;

      TH2D& xySieveHist = xySieveHists.at(iFoil);

      c1->cd();
      xySieveHist.Draw("colz");
      c1->Update();
      gPad->Update();

      // Fit the projections to get position estimates.
      c2->cd();
      tmpHist = xySieveHist.ProjectionX();
      tmpHist->SetTitle("x_{fp} projection");
      tmpHist->Draw();
      std::vector<Peak> xSievePeaksFit = fitMultiPeak(tmpHist, 0.1);
      gPad->Update();

      c3->cd();
      tmpHist = xySieveHist.ProjectionY();
      tmpHist->SetTitle("y_{fp} projection");
      tmpHist->Draw();
      std::vector<Peak> ySievePeaksFit = fitMultiPeak(tmpHist, 0.1);
      gPad->Update();

      if (cmdOpts.automatic) {
        std::this_thread::sleep_for(std::chrono::milliseconds(cmdOpts.delay));
      }
      else {
        cout << "    Continue? ";
        cin >> tmp;
      }

      // Setup before fitting.
      std::vector<Peak>& xSievePeaks = xSievePeakss.at(iFoil);
      std::vector<Peak>& ySievePeaks = ySievePeakss.at(iFoil);
      std::vector<std::size_t>& xSieveIndexes = xSieveIndexess.at(iFoil);
      std::vector<std::size_t>& ySieveIndexes = ySieveIndexess.at(iFoil);
      std::vector<std::size_t>& nEvents = nEventss.at(iFoil);
      std::vector<TEllipse>& ellipses = ellipsess.at(iFoil);

      // Fit each individual hole.
      for (const auto& xSievePeak : xSievePeaksFit) {
        tmpMark->SetX(xSievePeak.mean);
        for (const auto& ySievePeak : ySievePeaksFit) {
          tmpMark->SetY(ySievePeak.mean);
          c1->cd();
          tmpMark->Draw();
          gPad->Update();

          // Find bounding box for current hole.
          int binXmin = xySieveHist.GetXaxis()->FindBin(xSievePeak.mean - 3*xSievePeak.sigma);
          int binXmax = xySieveHist.GetXaxis()->FindBin(xSievePeak.mean + 3*xSievePeak.sigma);
          int binYmin = xySieveHist.GetYaxis()->FindBin(ySievePeak.mean - 3*ySievePeak.sigma);
          int binYmax = xySieveHist.GetYaxis()->FindBin(ySievePeak.mean + 3*ySievePeak.sigma);

          // Want to have at least 50 events for fitting.
          double integral = xySieveHist.Integral(
            binXmin, binXmax,
            binYmin, binYmax
          );
          if (integral < 100) continue;

          // Fit x and y projection separately.
          c2->cd();
          tmpHist = xySieveHist.ProjectionX("_px", binYmin, binYmax);
          tmpHist->GetXaxis()->SetRange(binXmin, binXmax);
          tmpHist->Draw();
          Peak xSievePeakSingle = fitPeak(
            tmpHist,
            xSievePeak.norm,
            xSievePeak.mean,
            xSievePeak.sigma
          );
          gPad->Update();

          c3->cd();
          tmpHist = xySieveHist.ProjectionY("_py", binXmin, binXmax);
          tmpHist->GetXaxis()->SetRange(binYmin, binYmax);
          tmpHist->Draw();
          Peak ySievePeakSingle = fitPeak(
            tmpHist,
            ySievePeak.norm,
            ySievePeak.mean,
            ySievePeak.sigma
          );
          gPad->Update();

          // Construct bounding ellipse.
          TEllipse ellipse(
            xSievePeakSingle.mean, ySievePeakSingle.mean,
            2.5*xSievePeakSingle.sigma, 2.5*ySievePeakSingle.sigma
          );
          ellipse.SetLineColor(2);
          ellipse.SetLineWidth(2);
          ellipse.SetFillStyle(0);

          // Push everything to collection.
          xSievePeaks.push_back(xSievePeakSingle);
          ySievePeaks.push_back(ySievePeakSingle);
          xSieveIndexes.push_back(getClosestIndex(xSievePeakSingle.mean, xSievePhys));
          ySieveIndexes.push_back(getClosestIndex(ySievePeakSingle.mean, ySievePhys));
          nEvents.push_back(0);
          ellipses.push_back(ellipse);
        }
      }

      c1->cd();
      tmpMark->SetX(1000.0);
      tmpMark->Draw();
      for (auto& ellipse : ellipses) {
        xySieveHist.GetListOfFunctions()->Add(&ellipse);
      }
      gPad->Update();
      xySieveHist.Write();

      c2->Clear();
      gPad->Update();
      c3->Clear();
      gPad->Update();

      if (cmdOpts.automatic) {
        std::this_thread::sleep_for(std::chrono::milliseconds(cmdOpts.delay));
      }
      else {
        cout << "    Continue? ";
        cin >> tmp;
      }

      c1->Clear();
      gPad->Update();
    }  // foil loop

    // Cleanup of sieve fit.
    delete tmpHist;
    delete tmpMark;


    cout << "    Filling SVD matrices and vectors: ";
    iEvent = 0;

    reportProgressInit();
    for (const auto& event : events) {  // SVD filling loop
      if (iEvent%1000 == 0) reportProgress(iEvent, nEvents);
      ++iEvent;

      // Find which foil if any.
      uint iFoil = 0;
      for (iFoil=0; iFoil<nFoils; ++iFoil) {
        if (
          zVerPeaks.at(iFoil).mean - 3.0*zVerPeaks.at(iFoil).sigma <= event.zVer &&
          event.zVer <= zVerPeaks.at(iFoil).mean + 3.0*zVerPeaks.at(iFoil).sigma &&
          yTarPeaks.at(iFoil).mean - 3.0*yTarPeaks.at(iFoil).sigma <= event.yTar &&
          event.yTar <= yTarPeaks.at(iFoil).mean + 3.0*yTarPeaks.at(iFoil).sigma
        ) {
          break;
        }
      }

      // Skip event if it is too far from any foil.
      if (iFoil == nFoils) continue;

      // Find which sieve hole if any.
      uint iHole = 0;
      for (iHole=0; iHole<xSievePeakss.at(iFoil).size(); ++iHole) {
        Peak& xSieveP = xSievePeakss.at(iFoil).at(iHole);
        Peak& ySieveP = ySievePeakss.at(iFoil).at(iHole);
        if (
          xSieveP.mean - 2.5*xSieveP.sigma <= event.xSieve &&
          event.xSieve <= xSieveP.mean + 2.5*xSieveP.sigma &&
          ySieveP.mean - 2.5*ySieveP.sigma <= event.ySieve &&
          event.ySieve <= ySieveP.mean + 2.5*ySieveP.sigma
        ) {
          break;
        }
      }

      // Skip event if it is too far from any hole or if there is enough events
      // from this hole already.
      if (
        iHole == xSievePeakss.at(iFoil).size() ||
        nEventss.at(iFoil).at(iHole) > 50
      ) {
        continue;
      }
      ++nEventss.at(iFoil).at(iHole);


      // Calculate the real or "physical" event quantities.
      double zFoil = runConf.zFoils.at(iFoil);

      double xTarVerPhy = -event.yVer;
      double yTarVerPhy = zFoil*sinTheta + event.xVer*cosTheta;
      double zTarVerPhy = zFoil*cosTheta - event.xVer*sinTheta;

      double xpTarPhy =
        (xSievePhys.at(xSieveIndexess.at(iFoil).at(iHole)) - xTarVerPhy) /
        (conf.sieve.z0 - zTarVerPhy);
      double ypTarPhy =
        (ySievePhys.at(ySieveIndexess.at(iFoil).at(iHole)) - yTarVerPhy) /
        (conf.sieve.z0 - zTarVerPhy);
      double xTarPhy = xTarVerPhy - xpTarPhy*zTarVerPhy - runConf.HMS.xMispointing;
      double yTarPhy = yTarVerPhy - ypTarPhy*zTarVerPhy - runConf.HMS.yMispointing;


      // Calculate contributions of xTar dependent terms.
      // Use old reconstruction matrix and xTarPhy.
      xpSum = 0.0;
      ySum = 0.0;
      ypSum = 0.0;

      for (const auto& line : recMatrix.matrix) {
        if (line.E_xTar == 0) continue;

        lambda =
          pow(event.xFp/100.0, line.E_x) *
          pow(event.xpFp, line.E_xp) *
          pow(event.yFp/100.0, line.E_y) *
          pow(event.ypFp, line.E_yp) *
          pow(xTarPhy/100.0, line.E_xTar);

        xpSum += line.C_Xp * lambda;
        ySum += line.C_Y * lambda;
        ypSum += line.C_Yp * lambda;
      }

      // Calculate lambdas for xTar independent terms.
      // Use new matrix and xTarPhy.
      std::vector<double> lambdas;
      for (const auto& line : recMatrixNew.matrix) {
        if (line.E_xTar == 0) {
          lambdas.push_back(
            pow(event.xFp/100.0, line.E_x) *
            pow(event.xpFp, line.E_xp) *
            pow(event.yFp/100.0, line.E_y) *
            pow(event.ypFp, line.E_yp) *
            pow(xTarPhy/100.0, line.E_xTar)
          );
        }
        else {
          lambdas.push_back(0.0);
        }
      }

      // Add lambda_i * lambda_j to (i,j)-th element of SVD matrices.
      // Add (_TarPhy - _Sum) to SVD vectors.
      // For xTar dependent terms lambdas are 0.0, so we need not care.
      Int_t i = 0;
      Int_t j = 0;
      for (const auto& lambda_i : lambdas) {
        j = 0;
        for (const auto& lambda_j : lambdas) {
          xpTarFitMat(i, j) += lambda_i * lambda_j;
          yTarFitMat(i, j) += lambda_i * lambda_j;
          ypTarFitMat(i, j) += lambda_i * lambda_j;

          ++j;
        }

        xpTarFitVec(i) += lambda_i * (xpTarPhy - xpSum);
        yTarFitVec(i) += lambda_i * (yTarPhy/100.0 - ySum);
        ypTarFitVec(i) += lambda_i * (ypTarPhy - ypSum);

        ++i;
      }
    }  // SVD filling loop
    reportProgressFinish();

    //break;  // TMP
  }  // run loop


  cout << "Post-processing SVD matrices and vectors." << endl;
  // For xTar dependent coefficients:
  //   set matrix diagonal values to 1.0 and offdiagonal values to 0.0;
  //   set vector values to 0.0.
  // Due to calculation of lambdas, only diagonal values need to be set now.
  for (Int_t i=0; i<recMatrixNewLen; ++i) {
    std::size_t i_ = static_cast<std::size_t>(i);

    if (recMatrixNew.matrix.at(i_).E_xTar == 0) continue;

    xpTarFitMat(i, i) = 1.0;
    yTarFitMat(i, i) = 1.0;
    ypTarFitMat(i, i) = 1.0;
  }

  cout << "Solving SVD problems:" << endl;
  TDecompSVD xpTarSVD(xpTarFitMat);
  TDecompSVD yTarSVD(yTarFitMat);
  TDecompSVD ypTarSVD(ypTarFitMat);

  bool xpTarSuccess = xpTarSVD.Solve(xpTarFitVec);
  cout << "  xpTar: " << (xpTarSuccess ? "success" : "failure") << endl;
  bool yTarSuccess = yTarSVD.Solve(yTarFitVec);
  cout << "  yTar: " << (yTarSuccess ? "success" : "failure") << endl;
  bool ypTarSuccess = ypTarSVD.Solve(ypTarFitVec);
  cout << "  ypTar: " << (ypTarSuccess ? "success" : "failure") << endl;


  cout << "Constructing new optics matrix." << endl;
  lineIt = recMatrixNew.begin();
  std::vector<RecMatrixLine>::iterator lineOldIt;
  for (Int_t iTerm=0; iTerm<recMatrixNewLen; ++iTerm) {
    if (lineIt->E_xTar == 0) {
      lineIt->C_Xp = xpTarFitVec(iTerm);
      lineIt->C_Y = yTarFitVec(iTerm);
      lineIt->C_Yp = ypTarFitVec(iTerm);
    }
    else {
      lineOldIt = recMatrix.findLine(*lineIt);
      if (lineOldIt != recMatrix.end()) {
        lineIt->C_Xp = lineOldIt->C_Xp;
        lineIt->C_Y = lineOldIt->C_Y;
        lineIt->C_Yp = lineOldIt->C_Yp;
      }
    }

    ++lineIt;
  }

  cout << "  matrix saved to: " << conf.recMatrixFileNameNew << endl;
  writeMatrixFile(conf.recMatrixFileNameNew, recMatrixNew);

  // Cleanup and exit.
  delete c3;
  delete c2;
  delete c1;

  return 0;
}
