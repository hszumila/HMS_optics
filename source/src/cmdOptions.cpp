#include "cmdOptions.hpp"

#include <cstring>
#include <iostream>
#include <stdexcept>


// Implementation of OptionParser_reconstruct.

cmdOptions::OptionParser_reconstruct::OptionParser_reconstruct() :
  displayHelp(false), automatic(false),
  rootFileName("out.root"), delay(2000),
  configFileName(), matrixFileName()
{}


cmdOptions::OptionParser_reconstruct::~OptionParser_reconstruct() {}


void cmdOptions::OptionParser_reconstruct::init(
  const int& argc, const char* const* argv
) {
  int operands = 0;

  // First check for -h flag and ignore others.
  for (int i=1; i<argc; ++i) {
    if (strcmp(argv[i], "-h") == 0) {
      displayHelp = true;
      return;
    }
  }

  for (int i=1; i<argc; ++i) {
    // Check for flags without arguments.
    if (strcmp(argv[i], "-a") == 0) {
      automatic = true;
    }
    // Check for flags with arguments.
    else if (strcmp(argv[i], "-o") == 0) {
      if (i == argc-1 || argv[i+1][0] == '-') {
        std::string errorMsg = "Missing operand after `" + std::string(argv[i]) + "`.";
        throw std::runtime_error(errorMsg.c_str());
      }
      else {
        rootFileName = std::string(argv[i+1]);
        ++i;
      }
    }
    else if (strcmp(argv[i], "-d") == 0) {
      try {
        delay = std::stoul(std::string(argv[i+1]));
        ++i;
      }
      catch (const std::invalid_argument& err) {
        std::string errorMsg = "Wrong type of operand after `" + std::string(argv[i]) + "` : `" + std::string(argv[i+1]) + "`.";
        throw std::runtime_error(errorMsg.c_str());
      }
    }
    // Check for invalid flags.
    else if (argv[i][0] == '-') {
      std::string errorMsg = "Invaid option `" + std::string(argv[i]) + "`.";
      throw std::runtime_error(errorMsg.c_str());
    }
    // Here are our two filenames.
    else if (operands == 0) {
      configFileName = std::string(argv[i]);
      ++operands;
    }
    else if (operands == 1) {
      matrixFileName = std::string(argv[i]);
      ++operands;
    }
    // Only two filenames :)
    else if (operands > 1) {
      std::string errorMsg = "Extra operand `" + std::string(argv[i]) + "`.";
      throw std::runtime_error(errorMsg.c_str());
    }
    // If anything else goes wrong...
    else {
      std::string errorMsg = "Something wrong here `" + std::string(argv[i]) + "`.";
      throw std::runtime_error(errorMsg.c_str());
    }
  }

  // Check if we got two filenames.
  if (operands != 2) {
    std::string errorMsg = "Missing operand after `" + std::string(argv[argc-1]) + "`.";
    throw std::runtime_error(errorMsg.c_str());
  }
}


void cmdOptions::OptionParser_reconstruct::printHelp() {
  std::cout << "Usage: reconstruct [OPTION]... CONFIG_F MATRIX_F" << std::endl << std::endl;
  std::cout << "CONFIG_F : configuration file name" << std::endl;
  std::cout << "MATRIX_F : reconstruction matrix file name" << std::endl << std::endl;
  std::cout << "[OPTION] :" << std::endl;
  std::cout << "  -h : display this help" << std::endl;
  std::cout << "  -a : proceed automatically, do not wait for input" << std::endl;
  std::cout << "  -o ROOTout : save output ROOT file to `ROOTout`" << std::endl;
  std::cout << "  -d DELAY : delay when showing key plots (in miliseconds)" << std::endl;
  std::cout << "             default is `2000`" << std::endl;
}


// Implementation of OptionParser_hmsOptics.

cmdOptions::OptionParser_hmsOptics::OptionParser_hmsOptics() :
  displayHelp(false), automatic(false),
  rootFileName("out.root"), delay(2000),
  configFileName()
{}


cmdOptions::OptionParser_hmsOptics::~OptionParser_hmsOptics() {}


void cmdOptions::OptionParser_hmsOptics::init(
  const int& argc, const char* const* argv
) {
  int operands = 0;

  // First check for -h flag and ignore others.
  for (int i=1; i<argc; ++i) {
    if (strcmp(argv[i], "-h") == 0) {
      displayHelp = true;
      return;
    }
  }

  for (int i=1; i<argc; ++i) {
    // Check for flags without arguments.
    if (strcmp(argv[i], "-a") == 0) {
      automatic = true;
    }
    // Check for flags with arguments.
    else if (strcmp(argv[i], "-o") == 0) {
      if (i == argc-1 || argv[i+1][0] == '-') {
        std::string errorMsg = "Missing operand after `" + std::string(argv[i]) + "`.";
        throw std::runtime_error(errorMsg.c_str());
      }
      else {
        rootFileName = std::string(argv[i+1]);
        ++i;
      }
    }
    else if (strcmp(argv[i], "-d") == 0) {
      try {
        delay = std::stoul(std::string(argv[i+1]));
        ++i;
      }
      catch (const std::invalid_argument& err) {
        std::string errorMsg = "Wrong type of operand after `" + std::string(argv[i]) + "` : `" + std::string(argv[i+1]) + "`.";
        throw std::runtime_error(errorMsg.c_str());
      }
    }
    // Check for invalid flags.
    else if (argv[i][0] == '-') {
      std::string errorMsg = "Invaid option `" + std::string(argv[i]) + "`.";
      throw std::runtime_error(errorMsg.c_str());
    }
    // Here is our one filename.
    else if (operands == 0) {
      configFileName = std::string(argv[i]);
      ++operands;
    }
    // Only one filename :)
    else if (operands > 0) {
      std::string errorMsg = "Extra operand `" + std::string(argv[i]) + "`.";
      throw std::runtime_error(errorMsg.c_str());
    }
    // If anything else goes wrong...
    else {
      std::string errorMsg = "Something wrong here `" + std::string(argv[i]) + "`.";
      throw std::runtime_error(errorMsg.c_str());
    }
  }

  // Check if we got one filename.
  if (operands != 1) {
    std::string errorMsg = "Missing operand after `" + std::string(argv[argc-1]) + "`.";
    throw std::runtime_error(errorMsg.c_str());
  }
}


void cmdOptions::OptionParser_hmsOptics::printHelp() {
  std::cout << "Usage: hms_optics [OPTION]... CONFIG_F" << std::endl << std::endl;
  std::cout << "CONFIG_F : configuration file name" << std::endl;
  std::cout << "[OPTION] :" << std::endl;
  std::cout << "  -h : display this help" << std::endl;
  std::cout << "  -a : proceed automatically, do not wait for input" << std::endl;
  std::cout << "  -o ROOTout : save output ROOT file to `ROOTout`" << std::endl;
  std::cout << "  -d DELAY : delay when showing key plots (in miliseconds)" << std::endl;
  std::cout << "             default is `2000`" << std::endl;
}
