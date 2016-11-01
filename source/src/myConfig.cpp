#include "myConfig.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

// BeamConfig implementation.

config::BeamConfig::BeamConfig() : x0(0.0), y0(0.0), xp0(0.0), yp0(0.0) {}


config::BeamConfig::~BeamConfig() {}


// HMSconfig implementation.

config::HMSconfig::HMSconfig() :
  thetaCentral(0.0), cosTheta(1.0), sinTheta(0.0),
  thetaOffset(0.0), phiOffset(0.0),
  xMispointing(0.14), yMispointing(0.0)
{}


config::HMSconfig::~HMSconfig() {}


// RunConfig implementation.

config::RunConfig::RunConfig() :
  runNumber(0), fileList(), cuts(""), zFoils(),
  beam(), HMS()
{}


config::RunConfig::~RunConfig() {}


// SieveConfig implementation.

config::SieveConfig::SieveConfig() :
  nRow(9), nCol(9),
  xHoleMin(-10.160), yHoleMin(-6.096), xHoleSpace(2.540), yHoleSpace(1.524),
  x0(0.0), y0(0.0), z0(166.032)
{}


config::SieveConfig::~SieveConfig() {}


// Config implementation.

config::Config::Config() :
  recMatrixFileNameOld(""), recMatrixFileNameNew(""),
  fitOrder(0), maxEventsPerHole(0), zFoilOffset(0.0),
  xTarCorrIterNum(1), runConfigs(), sieve()
{}


config::Config::~Config() {}


std::vector<double> config::Config::getSieveHolesX() const {
  std::vector<double> xSieveHoles(sieve.nRow);

  for (size_t i=0; i<sieve.nRow; ++i) {
    xSieveHoles.at(i) =
      sieve.xHoleMin + static_cast<double>(i)*sieve.xHoleSpace +
      sieve.x0 + runConfigs.at(0).HMS.xMispointing;
  }

  return xSieveHoles;
}


std::vector<double> config::Config::getSieveHolesY() const {
  std::vector<double> ySieveHoles(sieve.nCol);

  for (size_t i=0; i<sieve.nCol; ++i) {
    ySieveHoles.at(i) =
      sieve.yHoleMin + static_cast<double>(i)*sieve.yHoleSpace +
      sieve.y0 + runConfigs.at(0).HMS.yMispointing;
  }

  return ySieveHoles;
}


// Implementation of functions.

std::vector<std::string> config::tokenize(const std::string& line) {
  std::vector<std::string> tokens;
  std::string item;
  std::stringstream ss;
  ss.str(line);

  while(getline(ss, item, ' ')) {
    if (!item.empty()) tokens.push_back(item);
  }

  return tokens;
}


config::Config config::loadConfigFile(const std::string& fname) {
  Config conf;

  std::ifstream ifs(fname);
  std::string line;
  std::vector<std::string> tokens;

  if (!ifs.is_open()) {
    throw std::runtime_error("Could not open file: `"+fname+"`!");
  }

  while (getline(ifs, line)) {
    if (line.empty()) continue;

    tokens = tokenize(line);
    if (tokens[0][0] == '#') continue;
    if (tokens[0] == "endlist") break;

    if (tokens[0] == "newrun") {
      conf.runConfigs.push_back(RunConfig());
      conf.runConfigs.back().runNumber = stoi(tokens[1]);
    }
    else if (tokens[0] == "filelist") {
      conf.runConfigs.back().fileList = std::vector<std::string> (tokens.begin()+1, tokens.end());
    }
    else if (tokens[0] == "beampos") {
      conf.runConfigs.back().beam.x0 = stod(tokens[1]);
      conf.runConfigs.back().beam.y0 = stod(tokens[2]);
      conf.runConfigs.back().beam.xp0 = stod(tokens[3]);
      conf.runConfigs.back().beam.yp0 = stod(tokens[4]);
    }
    else if (tokens[0] == "thetaHMS") {
      conf.runConfigs.back().HMS.thetaCentral = stod(tokens[1]);
      conf.runConfigs.back().HMS.cosTheta = cos(conf.runConfigs.back().HMS.thetaCentral*TMath::DegToRad());
      conf.runConfigs.back().HMS.sinTheta = sin(conf.runConfigs.back().HMS.thetaCentral*TMath::DegToRad());
    }
    else if (tokens[0] == "zfoil") {
      for (size_t i=1; i<tokens.size(); ++i) {
        conf.runConfigs.back().zFoils.push_back(stod(tokens.at(i)));
      }
    }
    else if (tokens[0] == "cut") {
      conf.runConfigs.back().cuts = tokens[1];
    }
  }

  double thetaOffset;
  double phiOffset;
  double tmp;

  ifs >> conf.recMatrixFileNameOld;
  ifs >> conf.recMatrixFileNameNew;
  ifs >> thetaOffset >> phiOffset;
  ifs >> conf.fitOrder;
  ifs >> conf.maxEventsPerHole >> tmp;
  ifs >> conf.zFoilOffset;

  // Post-process.
  for (auto& runConf : conf.runConfigs) {
    runConf.HMS.thetaOffset = thetaOffset;
    runConf.HMS.phiOffset = phiOffset;

    for (auto& zFoil : runConf.zFoils) {
      zFoil += conf.zFoilOffset;
    }
  }

  return conf;
}
