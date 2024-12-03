#include "../../lib/include/Support.hh"
#include "../../lib/RootSupport/RootSupport.hh"

#include "Particle.hh"

using namespace ParticleSimulator;



constexpr int n_particles = 1e6;

void generate()
{
	CCheck();
	rs::utils::SetStyle();
	rs::utils::IgnoreWarning();

  auto f = rs::file::Create(result_dirpath / filename_ini);

  auto tree = std::make_unique<rs::tree::TreeHelper>({
		"rx", "ry", "rz", "vx", "vy", "vz", "lifetime"
	});


	auto pb = csp::ProgressBar(n_particles);
	for (int i = 0; i < n_particles; ++i, ++pb) {
		auto particle = Particle();

		tree->Getf("rx") = particle.rx_;
		tree->Getf("ry") = particle.ry_;
		tree->Getf("rz") = particle.rz_;
		tree->Getf("vx") = particle.vx_;
		tree->Getf("vy") = particle.vy_;
		tree->Getf("vz") = particle.vz_;
		tree->Getf("lifetime") = particle.lifetime_;

		tree->Fill();
	}

	f->cd();
  tree->Write();
}
