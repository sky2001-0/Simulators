#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"

#include "../lib/include/globals.hh"
#include "../lib/include/Support.hh"
#include "../lib/RootSupport/RootSupport.hh"
#include "../lib/Simulators.hh"

using namespace simulators_support;
using namespace rs::utils;


const double T_r = 610. * u::K;
const double T_z = 320. * u::K;
const double chirp_rate_1D = 0.49 * u::GHz / u::ns;

std::vector<double> finish_ratios = {1.5, 10., 20.};//{0., 1., 2., 3., 5.};
std::vector<std::pair<double, Color_t>> c_rate_ratios = {
  {1.0, kRed},
//   {0.9, kGreen},
  {0.8, kCyan},
//   {0.7, kBlue}
};

const std::string label = "3D"; // "1D" or "3D"



namespace ps
{
  constexpr double f_0 = (
    ortho2P0_binding_energy - ortho1S_binding_energy
  ) / u::h;
  constexpr double f_1 = (
    ortho2P1_binding_energy - ortho1S_binding_energy
  ) / u::h;
  constexpr double f_2 = (
    ortho2P2_binding_energy - ortho1S_binding_energy
  ) / u::h;

  constexpr double lambda_0 = u::c / f_0;
  constexpr double lambda_1 = u::c / f_1;
  constexpr double lambda_2 = u::c / f_2;

  constexpr double recoil_velocity = recoil / m_1S;
}


double calc_factor(const double velocity, const double temperature) {
  return (
    ps::m_1S / std::sqrt(1. - std::pow(velocity / u::c, 2)) - ps::m_1S
  ) * std::pow(u::c, 2) / (u::k_B * temperature);
}


double cooling_amount(
  const double cool_finish,
  const double chirp_rate_ratio_by_1D,
  const double time
)
{
  const double chirp_rate = chirp_rate_ratio_by_1D * chirp_rate_1D;

  const double v_finish = cool_finish * ps::recoil_velocity;
  const double v_start = v_finish + chirp_rate * time / ps::f_1 * u::c;

  const double ann = std::exp(-ps::ortho1S_annihilation_rate * time);

  if (label == "1D") {
    return ann * std::erf(std::sqrt(calc_factor(v_start, T_r)));
  } else if (label == "3D") {
    const double tmp_z = calc_factor(v_start, T_z);
    return (
      ann
      * (1. - std::exp(-calc_factor(v_start, T_r)))
      * (1. - (1. + tmp_z) * std::exp(-tmp_z))
    );
  } else {
    throw std::runtime_error("exception");
  }
}



void tmp_cooled_atom() {
  CCheck();
  SetStyle();
  IgnoreWarning();

  const std::string title = (
    "cooled_atom_study_" + label
    + "_Tr_" + csp::utils::ToString(T_r, 0) + "K"
    + "_Tz_" + csp::utils::ToString(T_z, 0) + "K"
  );

  const int n_points = 1000;
  const double time_max = 1. * u::us;
  const double time_step = time_max / n_points;
  const double y_max = (label == "1D") ? 0.15 : 0.003;

  auto out = rs::file::Create(result_dirpath / (title + ".root"));

  auto mg = std::make_unique<TMultiGraph>("multigraph", title.c_str());
  {
    mg->GetXaxis()->SetTitle("t (us)");
    mg->GetYaxis()->SetTitle("ratio");
    // mg->GetXaxis()->SetTitleSize(0.07); // x軸のタイトルのサイズ
    // mg->GetYaxis()->SetTitleSize(0.07); // y軸のタイトルのサイズ
    // mg->GetXaxis()->SetLabelSize(0.07); // x軸のラベルのサイズ
    // mg->GetYaxis()->SetLabelSize(0.07); // y軸のラベルのサイズ
    mg->GetXaxis()->CenterTitle(true); // x軸のタイトルを中央揃え
    mg->GetYaxis()->CenterTitle(true); // y軸のタイトルを中央
    // mg->GetXaxis()->SetTitleOffset(1.0); // x軸のタイトルオフセット
    // mg->GetYaxis()->SetTitleOffset(1.6); // y軸のタイトルオフセット
    rs::graph::SetLimitY(mg, {0., 0.012});
  }

  auto l = std::make_unique<TLegend>(0.7, 0.65, 0.9, 0.85);

  for (const auto& pair : c_rate_ratios) {
    const auto& c_rate = pair.first;
    const auto& color_base = pair.second;

    int color_index = 0;
    for (const double finish : finish_ratios) {
      const char* g_title = Form(
        "v_{0}/v_{r} = %.1f, R_{c} = %.1f", finish, c_rate
      );

      std::vector<double> x_vals;
      std::vector<double> y_vals;

      for (auto time : csp::Range(0., time_step, n_points)) {
        x_vals.push_back(time / u::us);
        y_vals.push_back(cooling_amount(finish, c_rate, time));
      }

      auto g = new TGraph(n_points, x_vals.data(), y_vals.data());
      g->SetTitle(g_title);
      g->SetLineWidth(2);
      g->SetLineColor(color_base + color_index);
      l->AddEntry(g, g_title, "l");

      mg->Add(g);

      out->cd();
      g->Write(Form("graph_c1_%.1f_c2_%.1f", finish, c_rate));

      color_index++;
    }
  }

  auto c = std::make_unique<TCanvas>("c", title.c_str(), 800, 600);
  {
    // c->SetLeftMargin(0.2);
    // c->SetRightMargin(0.08);
    // c->SetTopMargin(0.1);
    // c->SetBottomMargin(0.15);
    mg->Draw("AL");
    l->Draw();
  }

  c->SaveAs((result_dirpath / "cooled_atom_study_50K.png").c_str());
  c->SaveAs((result_dirpath / "cooled_atom_study_50K.pdf").c_str());

  out->cd();
  c->Write();
  mg->Write();

  out->Close();
}
