#include "myRecMatrix.hpp"

#include "myConfig.hpp"

#include <iomanip>
#include <stdexcept>
#include <string>

#include <fstream>
#include <iostream>
  using std::cout;
  using std::endl;
#include <sstream>


// RecMatrixLine implementation.

RecMatrixLine::RecMatrixLine() :
  C_Xp(), C_Y(), C_Yp(), C_D(),
  E_x(), E_xp(), E_y(), E_yp(), E_xTar()
{}


RecMatrixLine::RecMatrixLine(
  double C_Xp, double C_Y, double C_Yp, double C_D,
  int E_x, int E_xp, int E_y, int E_yp, int E_xTar
) :
  C_Xp(C_Xp), C_Y(C_Y), C_Yp(C_Yp), C_D(C_D),
  E_x(E_x), E_xp(E_xp), E_y(E_y), E_yp(E_yp), E_xTar(E_xTar)
{}


RecMatrixLine::~RecMatrixLine() {}


// RecMatrix implementation.

RecMatrix::RecMatrix() : header(), matrix() {}


RecMatrix::~RecMatrix() {}


size_t RecMatrix::size() const {
  return matrix.size();
}


void RecMatrix::addLine(const RecMatrixLine& line) {
  matrix.push_back(line);
}


void RecMatrix::addLine(
  double C_Xp, double C_Y, double C_Yp, double C_D,
  int E_x, int E_xp, int E_y, int E_yp, int E_xTar
) {
  matrix.push_back(
    RecMatrixLine(
      C_Xp, C_Y, C_Yp, C_D,
      E_x, E_xp, E_y, E_yp, E_xTar
    )
  );
}


// Iterators.
std::vector<RecMatrixLine>::iterator RecMatrix::begin() {
  return matrix.begin();
}


std::vector<RecMatrixLine>::iterator RecMatrix::end() {
  return matrix.end();
}


std::vector<RecMatrixLine>::iterator RecMatrix::findLine(
  int E_x, int E_xp, int E_y, int E_yp, int E_xTar
) {
  for (std::vector<RecMatrixLine>::iterator it=begin(); it!=end(); ++it) {
    if (
      it->E_x == E_x && it->E_xp == E_xp &&
      it->E_y == E_y && it->E_yp == E_yp &&
      it->E_xTar == E_xTar
    ) return it;
  }

  return end();
}


std::vector<RecMatrixLine>::iterator RecMatrix::findLine(
  const RecMatrixLine& line
) {
  for (std::vector<RecMatrixLine>::iterator it=begin(); it!=end(); ++it) {
    if (
      it->E_x == line.E_x && it->E_xp == line.E_xp &&
      it->E_y == line.E_y && it->E_yp == line.E_yp &&
      it->E_xTar == line.E_xTar
    ) return it;
  }

  return end();
}


// Other functions.

std::ostream& operator<<(std::ostream& os, const RecMatrixLine& line) {
  std::ios::fmtflags f(os.flags());
  std::streamsize prevPrec = os.precision(9);

  os
    << std::scientific
    << std::setw(17) << line.C_Xp
    << std::setw(17) << line.C_Y
    << std::setw(17) << line.C_Yp
    << std::setw(17) << line.C_D << "  "
    << line.E_x
    << line.E_xp
    << line.E_y
    << line.E_yp
    << line.E_xTar;

  os.precision(prevPrec);

  os.flags(f);

  return os;
}


RecMatrix readMatrixFile(const std::string& fileName) {
  std::ifstream ifs(fileName);
  if (!ifs.is_open()) {
    throw std::runtime_error("Could not open file: `" + fileName + "`!");
  }

  RecMatrix recMatrix = RecMatrix();

  std::string line;

  while (getline(ifs, line)) {
    if (line.find("!") == std::string::npos) break;
    recMatrix.header += line + '\n';
  }

  double Xp, Y, Yp, D;
  std::string exponents;

  while (getline(ifs, line)) {
    if (line.find("---") != std::string::npos) break;
    std::stringstream ss;
    ss.str(line);
    ss >> Xp >> Y >> Yp >> D >> exponents;

    recMatrix.addLine(
      Xp, Y, Yp, D,
      exponents[0] - '0', exponents[1] - '0', exponents[2] - '0',
      exponents[3] - '0', exponents[4] - '0'
    );
  }

  ifs.close();

  return recMatrix;
}


void writeMatrixFile(const std::string& fileName, const RecMatrix& recMatrix) {
  std::ofstream ofs(fileName);

  ofs << recMatrix.header;
  ofs << " ---------------------------------------------------------------------" << endl;
  for (const auto& line : recMatrix.matrix) {
    ofs << line << endl;
  }
  ofs << " ---------------------------------------------------------------------" << endl;

  ofs.close();
}
