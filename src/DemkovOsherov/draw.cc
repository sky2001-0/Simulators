#include <cmath>
#include <cstddef>
#include <iostream>

#include "../../lib/include/src/globals.h"
#include "../../lib/include/src/ProgressBar.h"
#include "../../lib/Simulators.h"

#include "Model.h"

using namespace simulators_support;
using namespace rs::utils;
using namespace DemkovOsherovModel;



constexpr std::size_t size_of_model = 2;
const std::size_t first_state = 0;
const int magnetic_qnumber = 1;
const bool with_decay = true;

const double adiabatic = 1.2;
const double rabi_rate = 2.;



void draw()
{
  CCheck();
  SetStyle();
  IgnoreWarning();

  const double rabi = rabi_rate * ps::split_02;
  const double chirp_rate = std::pow(rabi / adiabatic, 2) / 2.;

  std::cout
    << "Chirp rate = "
    << chirp_rate / u::tau / u::GHz * u::ns
    << " GHz/ns"
    << std::endl;

  auto sim = Simulator<size_of_model>::PsOrtho1S2P(
    chirp_rate, rabi, magnetic_qnumber, with_decay
  );

  ResultData result = {sim, chirp_rate, result_dirpath};
  {
    auto dmat = sim.InitialState(first_state);
    double time_n = sim.time_start();
    csp::ProgressBar<std::size_t> pb(sim.size());

    for (int i = 0; i < sim.size(); ++i, ++pb) {
      result.Set(dmat, time_n);
      sim.RK4(dmat, time_n);
    }
  }

  result.Save();
}
