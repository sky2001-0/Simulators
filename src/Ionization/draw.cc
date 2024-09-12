#include <iostream>
#include <string>

#include "../../lib/Simulators.hh"
#include "../../lib/include/globals.hh"
#include "../../lib/RootSupport/RootSupport.hh"

#include "TGraph.h"

#include "Calculator.hh"

using namespace simulators_support;
using namespace rs::utils;
using namespace Ionization;



void draw()
{
  CCheck();
  SetStyle();
  IgnoreWarning();

  auto func = [] (const int n, const int l) {
    const std::string state_name = to_string(n) + l_to_str(l);

    auto g = rs::graph::Create<TGraph>(
      100,
      state_name.c_str(),
      (state_name + " Ionization cross-section").c_str(),
      "light wavelength (nm)", "cross-section (Mb)"
    );

    Calculator c(ps_reduced_mass, n, l);

    auto ran = c.GetRange();
    for (int i = 0; i < 100; ++i) {
      const double light_lambda = ran[i];
      const double y = c.Calc(light_lambda);
      g->SetPoint(i, light_lambda / u::nm, y / (1e-22 * u::m * u::m));
    }

    g->SetMinimum(0.);

    rs::draw::FastSaveToFile(g, result_dirpath / (state_name + ".png"), "AL");
  };


  func(1, 0);
  func(2, 0);
  func(2, 1);

  Calculator c(ps_reduced_mass, 2, 1);
  std::cout << (
    c.Calc(243. * u::nm)
  ) / u::m / u::m << " m2" << std::endl;
}

