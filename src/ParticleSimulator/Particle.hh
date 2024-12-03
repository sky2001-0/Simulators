#ifndef PARTICLE_HH
#define PARTICLE_HH

#include <cmath>
#include <filesystem>
#include <functional>
#include <string>
#include <utility>

#include "TMath.h"

#include "../../lib/include/src/globals.h"
#include "../../lib/include/src/Support.h"
#include "../../lib/Simulators.hh"



namespace ps
{
	constexpr double v_recoil = ps::recoil / ps::m_1S;
}


namespace ParticleSimulator
{

const double T_r = 610. * u::K;
const double T_z = 320. * u::K;
const double r_ini_max = 1. * u::mm;

const std::filesystem::path result_dirpath = (
	simulators_support::result_dirpath / "particle_simulation"
);

void CCheck()
{
	simulators_support::CCheck();
	std::filesystem::create_directory(result_dirpath);
}

const std::string filename_ini = "initial_particles.root";


struct Particle
{
	double rx_;
	double ry_;
	double rz_;
	double vx_;
	double vy_;
	double vz_;
	double lifetime_;

private:
	static double calc_rr_random() noexcept
	{
		// f(r) \propto 1 / r
		const double uni = csp::math::Rand(0., 1.);
		return std::sqrt(uni) * u::mm;
	}

	static double calc_theta_random() noexcept
	{
		// f(r) \propto 1 / r
		return csp::math::Rand(0., u::tau);
	}

	static double calc_vxy_random() noexcept
	{
		// f(x) \propto \exp(-x^2 / sigma^2)

		const double uni = csp::math::Rand();
		const double sign = csp::math::Sign(csp::math::Rand(-1., 1.));
		const double sigma = std::sqrt(2. / ps::m_1S * u::k_B * T_r);
		return sign * TMath::ErfcInverse(uni) * sigma;
	}

	static double calc_vz_random() noexcept
	{
		// f(x) \propto x^3 \exp(-x^2 / \sigma^2)
		const double uni = csp::math::Rand();
		const double sigma = std::sqrt(2. / ps::m_1S * u::k_B * T_z);
		return (
			std::sqrt(-ROOT::Math::lambert_Wm1((uni - 1.) / std::exp(1.)) - 1.) * sigma
		);
	}

	static double calc_lifetime_random() noexcept
	{
		// \int_x^\infty f(x) = \exp(-ps::ortho1S_annihilation_rate * x)
		const double uni = csp::math::Rand();
		return -1. / ps::ortho1S_annihilation_rate * std::log(1. - uni);
	}


public:
	explicit Particle()
	{
		const double r = calc_rr_random();
		const double theta = calc_theta_random();
		rx_ = r * std::cos(theta);
		ry_ = r * std::sin(theta);
		rz_ = 0.;
		vx_ = calc_vxy_random();
		vy_ = calc_vxy_random();
		vz_ = calc_vz_random();
		lifetime_ = calc_lifetime_random();
	}

	explicit Particle(
		double rx,
		double ry,
		double rz,
		double vx,
		double vy,
		double vz,
		double lifetime
	)
	: rx_(rx), ry_(ry), rz_(rz), vx_(vx), vy_(vy), vz_(vz), lifetime_(lifetime)
	{
	}

  double Calc_v_abs() const noexcept
  {
    return std::sqrt(vx_ * vx_ + vy_ * vy_ + vz_* vz_);
  }

  void Move(const double time) noexcept
  {
    rx_ += vx_ * time;
    ry_ += vy_ * time;
    rz_ += vz_ * time;
  }
};


struct PulseLaser
{
	using RangeFromTime = std::function<std::pair<double, double>(const double)>;

	const double is_x_;
	const double is_y_;
	const double is_z_;
	const double center_x_;
	const double center_y_;
	const double center_z_;
	const double beam_radius2_;

	const double transition_probability_;

	const std::vector<double> laser_times_;

	const RangeFromTime captured_velocity_at_time_;

public:
	PulseLaser() = delete;

	PulseLaser(
		const bool is_x, const bool is_y, const bool is_z,
		const double center_x, const double center_y, const double center_z,
		const double beam_radius, const double transition_probability,
		const std::vector<double> laser_times,
		const RangeFromTime captured_velocity_at_time
	)
	: is_x_(is_x), is_y_(is_y), is_z_(is_z),
		center_x_(center_x), center_y_(center_y), center_z_(center_z),
		beam_radius2_(std::pow(beam_radius, 2)),
		transition_probability_(transition_probability),
		laser_times_(laser_times),
		captured_velocity_at_time_(captured_velocity_at_time)
	{
	}

	const std::vector<double>& RadiationTime() noexcept
	{
		return laser_times_;
	}

	void CoolParicle(Particle& particle, const double t)
	{
		const double& x = particle.rx_;
		const double& y = particle.ry_;
		const double& z = particle.rz_;
		double& vx = particle.vx_;
		double& vy = particle.vy_;
		double& vz = particle.vz_;
		auto [v_min, v_max] = captured_velocity_at_time_(t);

		if (
			is_z_
			&& (
				std::pow(x - center_x_, 2) + std::pow(y - center_y_, 2)
				< beam_radius2_
			)
			&& v_min <= std::abs(vz) && std::abs(vz) < v_max
			&& csp::math::Rand() <= transition_probability_
		) {
			vz -= csp::math::Sign(vz) * ps::v_recoil;

		} else if (
			is_x_
			&& (
				std::pow(y - center_y_, 2) + std::pow(z - center_z_, 2)
				< beam_radius2_
			)
			&& v_min <= std::abs(vx) && std::abs(vx) < v_max
			&& csp::math::Rand() <= transition_probability_
		) {
			vx -= csp::math::Sign(vx) * ps::v_recoil;

		} else if (
			is_y_
			&& (
				std::pow(z - center_z_, 2) + std::pow(x - center_x_, 2)
				< beam_radius2_
			)
			&& v_min <= std::abs(vy) && std::abs(vy) < v_max
			&& csp::math::Rand() <= transition_probability_
		) {
			vy -= csp::math::Sign(vy) * ps::v_recoil;

		}
	}
};

}



#endif // PARTICLE_HH
