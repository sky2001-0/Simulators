#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "TF1.h"
#include "TTree.h"

#include "../../lib/include/src/globals.h"
#include "../../lib/include/src/ProgressBar.h"
#include "../../lib/include/src/Range.h"
#include "../../lib/RootSupport/RootSupport.h"
#include "../../lib/Simulators.h"

#include "Model.h"

using namespace simulators_support;
using namespace rs::utils;
using namespace DemkovOsherovModel;



constexpr std::size_t size_of_model = 2;
const std::size_t first_state = 0;
const int magnetic_qnumber = 1;
const bool with_decay = true;

const std::vector<std::string> branch_names = {
  "adiabatic/D", "rabi_rate/D", "max_prob/D", "max_time_ns/D"
};

const auto adiabatic_range = csp::Range<double>(0.5, 2.9, 36, false);
const std::vector<double> rabi_rate_range = (
  []()
  {
    std::vector<double> result;
    for (auto d : csp::Range<double>(-2., 2., 21, true)) {
      result.push_back(std::pow(10., d));
    }
    return result;
  }()
);



void search()
{
  CCheck();
  SetStyle();
  IgnoreWarning();

  rs::TreeHelper th = {branch_names};

  csp::ProgressBar<std::size_t> pb(
    adiabatic_range.size() * rabi_rate_range.size()
  );

  for (auto adiabatic : adiabatic_range) {
    th.get("adiabatic/D") = adiabatic;
    for (auto rabi_rate : rabi_rate_range) {
      th.get("rabi_rate/D") = rabi_rate;

      const double rabi = rabi_rate * ps::split_02;
      const double chirp_rate = std::pow(rabi / adiabatic, 2) / 2.;

      {
        auto sim = Simulator<size_of_model>::PsOrtho1S2P(
          chirp_rate, rabi, magnetic_qnumber, with_decay
        );
        auto dmat = sim.InitialState(first_state);
        double time_n = sim.time_start();

        double max_prob = 0.;
        double max_time_n = 0.;

        for (int i = 0; i < sim.size(); ++i) {
          sim.RK4(dmat, time_n);
          const double prob = std::abs(dmat.getf(1, 1));
          if (prob > max_prob) {
            max_prob = prob;
            max_time_n = time_n;
          }
        }

        th.get("max_time_ns/D") = max_time_n / std::sqrt(chirp_rate) / u::ns;
        th.get("max_prob/D") = max_prob;

        // double correctly = std::abs(dmat.trace());
      }

      th.Fill();
      ++pb;
    }
  }

  auto file = rs::file::Create(
    result_dirpath / ("tree" + std::to_string(size_of_model) + ".root")
  );

  th.Write();
}
