#ifndef cmdOptions_h
#define cmdOptions_h 1

#include <string>


//! Interface for dealing with command line arguments.
namespace cmdOptions {

  class OptionParser_reconstruct {
    public:
      OptionParser_reconstruct();
      ~OptionParser_reconstruct();

      void init(const int& argc, const char* const* argv);
      void printHelp();

      bool displayHelp;
      bool automatic;

      std::string rootFileName;
      unsigned long delay;

      std::string configFileName;
      std::string matrixFileName;
  };

  class OptionParser_hmsOptics {
    public:
      OptionParser_hmsOptics();
      ~OptionParser_hmsOptics();

      void init(const int& argc, const char* const* argv);
      void printHelp();

      bool displayHelp;
      bool automatic;

      std::string rootFileName;
      unsigned long delay;

      std::string configFileName;
  };

}

#endif  // cmdOptions_h
