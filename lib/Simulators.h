#ifndef SIMULATORS_H
#define SIMULATORS_H

#include <filesystem>
#include <iostream>
#include <string>



namespace simulators_support
{
  const std::filesystem::path result_dirpath = "result";

  void CCheck()
  {
    const std::filesystem::path current_path = std::filesystem::current_path();
    if (!std::filesystem::exists(current_path / "lib" / "Simulators.h")) {
      std::cout
        << "Current working path is not expected one."
        << "Do you continue? [y/N]";
      std::string input;
      std::cin >> input;
      if (input != "y") {
        std::exit(EXIT_FAILURE);
      }
    }

    std::filesystem::create_directory(result_dirpath);
  }
}



#endif // SIMULATORS_H
