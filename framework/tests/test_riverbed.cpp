#include "catch.hpp"
#include "../geometry/emitter.h"
#include "../geometry/io.h"
#include "../SPH/sph.h"
#include "../SPH/params.h"
#include "../PBF/pbf.h"
#include "../PBF/params.h"


// Check out https://github.com/catchorg/Catch2 for more information about how to use Catch2
/*TEST_CASE("Riverbed SPH", "[Riverbed SPH]") {
	WCSPH::Vector boundary_size = { 5.0, 5.0, 2.0 };
	WCSPH::Vector boundary_left = { -2.5, -2.5, -1 };
	std::vector<WCSPH::Vector> fluid_sizes;
	std::vector<WCSPH::Vector> fluid_lefts;

	SPHParameters params;
	MCParameters mcparams;
	params.dt = 0.00025;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.0075;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.00025;
	params.max_velocity_cap = 5;
	params.fluid_pressure_stiffness = 1000.0;
	params.fluid_viscosity = 0.0025;
	params.boundary_viscosity = 0.0;

	params.particle_diameter = 2 * params.particle_radius;
	params.fluid_sampling_distance = params.particle_diameter;
	params.boundary_sampling_distance = 0.8 * params.particle_diameter;
	params.smoothing_length = 1.2 * params.particle_diameter;
	params.smoothing_length_squared = params.smoothing_length * params.smoothing_length;
	params.compact_support = 2 * params.smoothing_length;
	params.emit_frequency = 1.2;
	params.max_num_particles = 1500000;
	params.grid_cell_size = 1.2 * params.particle_radius;
	params.export_type = params.EXPORT_WITH_SURFACE;

	Emitter::Emitter emitter1(0.1, { 0.0, 1.0, 0.0 }, 0.1, { 0.0, 0.0, 1.0 }, { -2, -1.4, 0.1 }, { 1.0, 0.0, 0.0 }, { 2.0, 0.0, 0.0 }, params.particle_diameter, params.particle_diameter);

	std::vector<CompactNSearch::Real> schedule;
	schedule.push_back(0.0);
	schedule.push_back(2.1);
	emitter1.set_schedule(schedule);
	std::vector<Emitter::Emitter> emitters;
	emitters.push_back(emitter1);

	WCSPH::SPH sph(false, false, true, "../res/riverbed_sph/riverbed_sph_", params, mcparams, emitters);
	sph.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	sph.simulate(4);
	sph.printStats("Riverbed SPH");
}

TEST_CASE("Riverbed PBF", "[Riverbed PBF]") {
	PBD::Vector boundary_size = { 5.0, 5.0, 2.0 };
	PBD::Vector boundary_left = { -2.5, -2.5, -1.0 };
	std::vector<PBD::Vector> fluid_sizes;
	std::vector<PBD::Vector> fluid_lefts;

	PBFParameters params;
	MCParameters mcparams;
	params.dt = 0.00025;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.01;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.00025;
	params.max_velocity_cap = 5;
	params.fluid_pressure_stiffness = 1000.0;
	params.fluid_viscosity = 0.0025;
	params.boundary_viscosity = 0.0;
	params.pbf_iterations = 5;

	params.particle_diameter = 2 * params.particle_radius;
	params.fluid_sampling_distance = params.particle_diameter;
	params.boundary_sampling_distance = 0.8 * params.particle_diameter;
	params.smoothing_length = 1.2 * params.particle_diameter;
	params.smoothing_length_squared = params.smoothing_length * params.smoothing_length;
	params.compact_support = 2 * params.smoothing_length;
	params.emit_frequency = 1.2;
	params.max_num_particles = 1500000;
	params.grid_cell_size = 1.2 * params.particle_radius;
	params.export_type = params.EXPORT_WITH_SURFACE;
	params.pbf_iterations = 10;

	Emitter::Emitter emitter1(0.1, { 0.0, 1.0, 0.0 }, 0.1, { 0.0, 0.0, 1.0 }, { -2, -1.4, 0.25 }, { 1.0, 0.0, 0.0 }, { 2, 0.0, 0.0 }, params.particle_diameter, params.particle_diameter);

	std::vector<CompactNSearch::Real> schedule;
	schedule.push_back(0.0);
	schedule.push_back(2.1);
	emitter1.set_schedule(schedule);
	std::vector<Emitter::Emitter> emitters;
	emitters.push_back(emitter1);

	PBD::PBF pbf(false, false, true, "../res/riverbed_pbf/riverbed_pbf_", params, mcparams, emitters);
	pbf.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	pbf.simulate(4);
	pbf.printStats("Riverbed PBF");
}*/