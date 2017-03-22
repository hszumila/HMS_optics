#ifndef myEvent_h
#define myEvent_h 1

#include <vector>

#include "myConfig.hpp"


class Event {
  public:
    Event();
    ~Event();

    void reset();

    // Focal plane variables.
    double xFp;  // cm
    double yFp;  // cm
    double xpFp;
    double ypFp;

    // Vertex variables in laboratory system.
    double xVer;  // cm
    double yVer;  // cm
    double zVer;  // cm

    // Reconstructed variables in spectrometer target system.
    double xTar;  // cm
    double yTar;  // cm
    double xpTar;
    double ypTar;

    // Vertex variables in spectrometer target system.
    double xTarVer;  // cm
    double yTarVer;  // cm
    double zTarVer;  // cm

    // Sieve variables.
    double xSieve;  // cm
    double ySieve;  // cm
};


std::vector<Event> readEvents(const config::RunConfig& runConf);


#endif  // myEvent_h
