#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

#include "TCanvas.h"
#include "TEllipse.h"
#include "TTree.h"

#include "../../lib/include/Support.hh"
#include "../../lib/RootSupport/RootSupport.hh"

#include "Particle.hh"

using namespace ParticleSimulator;



namespace ls
{
	constexpr double T_fin = 10. * u::K;
	constexpr double f_rep = 236.25 * u::MHz;
	constexpr double pulse_fullwidth = 100e9; //8.9 * u::GHz;
	constexpr double chirp_rate = 0.; //4.9 * u::GHz / u::ns;
	constexpr double t_cool_start = 0. * u::ns;
	constexpr double t_cool_end = 100. * u::ns;
	const double f_end = std::max(
		std::sqrt(u::k_B * T_fin / ps::m_1S) / u::c * ps::trans_angfreq / u::tau,
		pulse_fullwidth / 2.
	);

	constexpr double trans_prob = 1.;

	constexpr double center_x = 0. * u::mm;
	constexpr double center_y = 0. * u::mm;
	constexpr double center_z = 2.5 * u::mm;
	constexpr double diameter = 5. * u::mm;
}

namespace calc
{
	constexpr double t_start = 0. * u::ns;
	constexpr double t_step = 2. * u::ns;
	constexpr double t_end = 150. * u::ns;
}

namespace display
{
	constexpr int xy_bins = 4e3;
	constexpr int z_bins = 4e3;
	constexpr double xy_min = -20. * u::mm;
	constexpr double xy_max = 20. * u::mm;
	constexpr double z_min = 0. * u::mm;
	constexpr double z_max = 30. * u::mm;
	constexpr double T_min = 1. * u::K;
	constexpr double T_max = 1e3 * u::K;
}



void simulate()
{
	CCheck();
	rs::utils::SetStyle();
	gStyle->SetOptStat(0);
	rs::utils::IgnoreWarning();

	// Derete png

	const std::string cmd1 = (
		"rm -f " + (result_dirpath / "distribution_*.png").string()
	);
	std::system(cmd1.c_str());

	// Laser Setting

	std::vector<double> laser_time;
	{
		for (double t = ls::t_cool_start; t < ls::t_cool_end; t += 1. / ls::f_rep) {
			laser_time.push_back(t);
		}
	}

	auto capt = [](const double t) {
		const double t_cool = ls::t_cool_end - t;
		const double f_min = ls::chirp_rate * t_cool + ls::f_end;
		const double f_max = f_min + ls::pulse_fullwidth;
		const double v_min = f_min * u::tau / ps::trans_angfreq * u::c;
		const double v_max = f_max * u::tau / ps::trans_angfreq * u::c;
		return std::make_pair(v_min, v_max);
	};

	PulseLaser laser(
		true, true, true,
		ls::center_x, ls::center_y, ls::center_z, ls::diameter / 2.,
		ls::trans_prob, laser_time, capt
	);

	// Input

	auto f_in = rs::file::Open(result_dirpath / filename_ini);
	auto tree_raw = rs::file::GetObj<TTree>("tree", f_in);
	auto tree = std::make_unique<rs::tree::TreeHelper>(std::move(tree_raw));

  const int n_particles = tree->GetEntries();
	std::vector<Particle> particles;
	auto pb1 = csp::ProgressBar(n_particles);
	for (int entry = 0; entry < n_particles; ++entry, ++pb1) {
		tree->GetEntry(entry);
		particles.emplace_back(
			tree->Get("rx"),
			tree->Get("ry"),
			tree->Get("rz"),
			tree->Get("vx"),
			tree->Get("vy"),
			tree->Get("vz"),
			tree->Get("lifetime")
		);
	}

	// For Drawing

	auto h_distribution_xz = std::make_unique<TH2F>(
		"hXZ", "",
		display::xy_bins, display::xy_min / u::mm, display::xy_max / u::mm,
		display::z_bins, display::z_min / u::mm, display::z_max / u::mm
	);
	{
		h_distribution_xz->GetXaxis()->SetTitle("position x (mm)");
		h_distribution_xz->GetYaxis()->SetTitle("position z (mm)");
		h_distribution_xz->GetZaxis()->SetTitle("equivalent temperature (K)");
		h_distribution_xz->SetMinimum(display::T_min);
		h_distribution_xz->SetMaximum(display::T_max);
	}

	auto circle = std::make_unique<TEllipse>(
		laser.center_x_ / u::mm, laser.center_z_ / u::mm,
		std::sqrt(laser.beam_radius2_) / u::mm
	);
	circle->SetLineWidth(1);
	circle->SetFillStyle(0);

	// Simulation

	std::vector<double> ltimes = laser.RadiationTime();

	auto pb2 = csp::ProgressBar<int>(
		std::ceil((calc::t_end - calc::t_start) / calc::t_step)
	);

	for (double t = calc::t_start; t < calc::t_end; t += calc::t_step, ++pb2) {
		h_distribution_xz->Reset();
		h_distribution_xz->SetTitle(
			("XZ : time = " + csp::utils::ToString(t / u::ns, 1) + " ns").c_str()
		);

		std::vector<double> pulse_times;
		for (const double ltime : ltimes) {
			if (ltime < t - calc::t_step) {
				continue;
			}
			if (t <= ltime) {
				break;
			}
			pulse_times.push_back(ltime);
		}

		for (auto& particle : particles) {
			if (t > particle.lifetime_) {
				continue;
			}

			if (t > calc::t_start) {
				double t_before = t - calc::t_step;
				for (const double pulse_time : pulse_times) {
					particle.Move(pulse_time - t_before);
					laser.CoolParicle(particle, pulse_time);
					t_before = pulse_time;
				}
				particle.Move(t - t_before);
			}

			const int bin = h_distribution_xz->FindBin(
				particle.rx_ / u::mm, particle.rz_ / u::mm
			);
    	if (h_distribution_xz->GetBinContent(bin) == 0.) {
				h_distribution_xz->Fill(
					particle.rx_ / u::mm,
					particle.rz_ / u::mm,
					ps::m_1S * std::pow(particle.Calc_v_abs(), 2) / 3. / u::k_B
					// ps::m_1S * std::pow(particle.vx_, 2) / u::k_B
				);
			}
		}

		auto c = std::make_unique<TCanvas>();
		{
			c->SetRightMargin(0.2);
			h_distribution_xz->Draw("COLZ");
			circle->Draw("same");
			c->SetLogz();
			c->SaveAs(
				(
					result_dirpath
					/ ("distribution_"	+ std::to_string(pb2.get_progress()) + ".png")
				).c_str()
			);
		}
	}

	// Create gif

	std::string cmd = (
		"ffmpeg -loglevel error -y -framerate 10 -i "
		+ (result_dirpath / "distribution_%d.png").string()
		+ " -loop -1 -filter_complex 'split[a],palettegen,[a]paletteuse' "
		+ (result_dirpath / "distribution.gif").string()
	);
	const int result = std::system(cmd.c_str());
	if (result != 0) {
		std::cerr << "This cmd failed." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

