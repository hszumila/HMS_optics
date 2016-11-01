#ifndef myRecMatrixIO_h
#define myRecMatrixIO_h 1

#include <iostream>
#include <string>
#include <vector>


class RecMatrixLine {
  public:
    RecMatrixLine();
    RecMatrixLine(
      double C_Xp, double C_Y, double C_Yp, double C_D,
      int E_x, int E_xp, int E_y, int E_yp, int E_xTar
    );
    ~RecMatrixLine();

    friend std::ostream& operator<<(std::ostream& os, const RecMatrixLine& line);

    double C_Xp;
    double C_Y;
    double C_Yp;
    double C_D;

    int E_x;
    int E_xp;
    int E_y;
    int E_yp;
    int E_xTar;
};


class RecMatrix {
  public:
    RecMatrix();
    ~RecMatrix();

    size_t size() const;

    void addLine(const RecMatrixLine& line);
    void addLine(
      double C_Xp, double C_Y, double C_Yp, double C_D,
      int E_x, int E_xp, int E_y, int E_yp, int E_xTar
    );

    std::vector<RecMatrixLine>::iterator begin();
    std::vector<RecMatrixLine>::iterator end();
    std::vector<RecMatrixLine>::iterator findLine(
      int E_x, int E_xp, int E_y, int E_yp, int E_xTar
    );
    std::vector<RecMatrixLine>::iterator findLine(const RecMatrixLine& line);

    std::string header;
    std::vector<RecMatrixLine> matrix;
};


std::ostream& operator<<(std::ostream& os, const RecMatrixLine& line);

RecMatrix readMatrixFile(const std::string& fileName);
void writeMatrixFile(const std::string& fileName, const RecMatrix& recMatrix);

#endif  // myRecMatrixIO_h
