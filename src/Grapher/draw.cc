#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>

#include "../../lib/Simulators.h"
#include "../../lib/include/src/Range.h"
#include "../../lib/RootSupport/RootSupport.h"

using namespace simulators_support;
using namespace rs::utils;



auto ran = csp::Range<double>(-10., 10., 1000, true);


double func(const double x)
{
  return std::sqrt(std::pow(x, 2) + 1.);
}


void draw()
{
  CCheck();
  SetStyle();
  IgnoreWarning();

  auto g1 = new TGraph(ran.size());
  for (auto itr = ran.begin(); itr != ran.end(); ++itr) {
    const double x = *itr;
    g1->SetPoint(itr(), x, func(x));
  }

  auto mg = std::make_unique<TMultiGraph>();
  mg->Add(g1);

  rs::graph::SetLimit(mg.get(), {-11., 11.}, {-0.1, 1.6});

  auto c = std::make_unique<TCanvas>(mg->GetName(), mg->GetTitle());
  mg->Draw("AL");
  // rs::draw::DrawLineHorizontal(c.get(), 0.);
  c->SaveAs((result_dirpath / "graph.png").c_str());
}
