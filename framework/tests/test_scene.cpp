#include "catch.hpp"
#include "../geometry/emitter.h"
#include "../geometry/io.h"
#include "../SPH/sph.h"
#include "../SPH/params.h"
#include "../PBF/pbf.h"
#include "../PBF/params.h"
#include <cmath>

// Check out https://github.com/catchorg/Catch2 for more information about how to use Catch2
/*TEST_CASE("Test Emitter Fountain SPH Open Box", "[Test Emitter Fountain SPH Open Box]")
{
	//WCSPH::Vector boundary_size = { 3.0, 3.0, 0.5 }; // Use this in general
	WCSPH::Vector boundary_size = { 3.0, 3.0, 1.0 }; // Use this for with obstacle
	WCSPH::Vector boundary_left = { -0.015, -0.015, -0.015 };
	std::vector<WCSPH::Vector> fluid_sizes;

	std::vector<WCSPH::Vector> fluid_lefts;

	SPHParameters params;
	MCParameters mcparams;
	params.dt = 0.00025;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.005;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.00025;
	params.max_velocity_cap = 5;
	params.fluid_pressure_stiffness = 1000.0;
	params.fluid_viscosity = 0.005;
	params.boundary_viscosity = 0.0;

	params.particle_diameter = 2 * params.particle_radius;
	params.fluid_sampling_distance = params.particle_diameter;
	params.boundary_sampling_distance = 0.8 * params.particle_diameter;
	params.smoothing_length = 1.2 * params.particle_diameter;
	params.smoothing_length_squared = params.smoothing_length * params.smoothing_length;
	params.compact_support = 2 * params.smoothing_length;
	params.grid_cell_size = 1.2 * params.particle_radius;
	params.emit_frequency = 1.2;
	params.max_num_particles = 100000;

	Emitter::Emitter emitter(0.0625, { 1.0, 0.0, 0.0 }, 0.0625, { 0.0, 1.0, 0.0 }, { 1.5, 1.5, 0.1 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 5 }, params.particle_diameter, params.particle_diameter);

	SECTION("Fountain SPH Surface with surface tension") {
		params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
		mcparams.ours = false;
		mcparams.sparse = false;
		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.0);
		schedule.push_back(2.1);
		emitter.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter);

		WCSPH::SPH sph(false, false, true, "../res/emitter_fountain_sph_surface_ob_st/emitter_fountain_sph_surface_ob_st_", params, mcparams, emitters);
		sph.load_geometry(true, true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
	}

	Emitter::Emitter emitter1(0.0625, { -0.5, -1.0, -1.0 }, 0.0625, { -1.0, -0.5, 1.0 }, { 1.65, 1.35, 2.0 }, { -1.5, 1.5, -0.75 }, { -1.5, 1.5, -0.75 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter2(0.0625, { 0.5, 1.0, -1.0 }, 0.0625, { 1.0, 0.5, 1.0 }, { 1.35, 1.65, 2.0 }, { 1.5, -1.5, -0.75 }, { 1.5, -1.5, -0.75 }, params.particle_diameter, params.particle_diameter);
	WCSPH::Vector obstacle_bl = { 1.0, 1.0, 0.6 };
	WCSPH::Vector obstacle_size = { 1.0, 1.0, 0.25 };
	std::pair <WCSPH::Vector, WCSPH::Vector>& obstacle_box1 = std::pair <WCSPH::Vector, WCSPH::Vector>();
	obstacle_box1.first = obstacle_bl;
	obstacle_box1.second = obstacle_size;
	SECTION("Two hoses and an obstacle with surface tension in SPH") {
		params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
		mcparams.ours = false;
		mcparams.sparse = false;
		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.0);
		schedule.push_back(2.1);
		emitter1.set_schedule(schedule);
		emitter2.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter1);
		emitters.push_back(emitter2);

		WCSPH::SPH sph(false, false, true, "../res/emitter_2hose_1obstacle_sph_surface_ob_st/emitter_2hose_1obstacle_sph_surface_ob_st_", params, mcparams, emitters);
		sph.load_geometry(true, true, boundary_size, boundary_left, fluid_sizes, fluid_lefts, obstacle_box1);
		sph.simulate(1.8);
	}
	
	WCSPH::Vector boundary_size_v = { 3.0, 3.0, 0.75 };
	WCSPH::Vector boundary_left_v = { 0.0, 0.0, 0.0 };
	Emitter::Emitter emitter3(0.075, { 1.0, 0.0, 0.0 }, 0.075, { 0.0, 1.0, 0.0 }, { 1.5, 1.5, 1.85 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 3.0 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter4(0.075, { 0.0, 1.0, 0.0 }, 0.075, { 1.0, 0.0, 0.0 }, { 1.5, 1.5, 2.15 }, { 0.0, 0.0, -1.0 }, { 0.0, 0.0, -3.0 }, params.particle_diameter, params.particle_diameter);
	SECTION("Two vertical hoses with surface tension in SPH") {
		params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
		mcparams.ours = false;
		mcparams.sparse = false;
		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.0);
		schedule.push_back(2.1);
		emitter3.set_schedule(schedule);
		emitter4.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter3);
		emitters.push_back(emitter4);

		WCSPH::SPH sph(false, false, true, "../res/emitter_2vertical_hoses_sph_surface_ob_st/emitter_2vertical_hoses_sph_surface_ob_st_", params, mcparams, emitters);
		sph.load_geometry(true, true, boundary_size_v, boundary_left_v, fluid_sizes, fluid_lefts);
		sph.simulate(1.2);
	}
}

TEST_CASE("Test Emitter Fountain PBF", "[Test Emitter Fountain PBF]")
{
	PBD::Vector boundary_size = { 3.0, 3.0, 3.0 };
	PBD::Vector boundary_left = { -0.015, -0.015, -0.015 };
	std::vector<PBD::Vector> fluid_sizes;

	std::vector<PBD::Vector> fluid_lefts;

	PBFParameters params;
	MCParameters mcparams;
	params.dt = 0.0005;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.005;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.001;
	params.max_velocity_cap = 5;
	params.fluid_pressure_stiffness = 1000.0;
	params.fluid_viscosity = 0.005;
	params.boundary_viscosity = 0.0;

	params.particle_diameter = 2 * params.particle_radius;
	params.fluid_sampling_distance = params.particle_diameter;
	params.boundary_sampling_distance = 0.8 * params.particle_diameter;
	params.smoothing_length = 1.2 * params.particle_diameter;
	params.smoothing_length_squared = params.smoothing_length * params.smoothing_length;
	params.compact_support = 2 * params.smoothing_length;
	params.grid_cell_size = 1.2 * params.particle_radius;
	params.emit_frequency = 1.2;
	params.max_num_particles = 100000;

	params.export_type = PBFParameters::EXPORT_WITH_SURFACE;
	mcparams.ours = false;
	mcparams.sparse = false;

	//Emitter::Emitter emitter1(0.0625, { 1.0, 0.0, 0.0 }, 0.0625, { 0.0, 0.0, 1.0 }, { 1.5, 0.0, 1.5 }, { 0.0, 1.0, 0.0 }, { 0.0, 3.0, 0.0 }, params.particle_diameter, params.particle_diameter);
	//Emitter::Emitter emitter2(0.0625, { 0.0, 1.0, 0.0 }, 0.0625, { 0.0, 0.0, 1.0 }, { 0.0, 1.5, 1.5 }, { 1.0, 0.0, 0.0 }, { 3.0, 0.0, 0.0 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter(0.0625, { 1.0, 0.0, 0.0 }, 0.0625, { 0.0, 1.0, 0.0 }, { 1.5, 1.5, 0.1 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 5.0 }, params.particle_diameter, params.particle_diameter);

	std::vector<CompactNSearch::Real> schedule;
	schedule.push_back(0.0);
	schedule.push_back(2.1);
	emitter.set_schedule(schedule);
	std::vector<Emitter::Emitter> emitters;
	emitters.push_back(emitter);

	PBD::PBF pbf(false, false, true, "../res/emitter_fountain_pbf/emitter_fountain_pbf_", params, mcparams, emitters);
	pbf.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	pbf.simulate(2);
}*/

TEST_CASE("Italian Fountain SPH", "[Italian Fountain SPH]")
{
	std::cout << "Testing Italian Fountain in SPH" << std::endl;
	WCSPH::Vector fluid1_size = { 2.89497 - 0.06, 0.4632 - 0.06, 0.15 };
	WCSPH::Vector fluid1_left = { 6.03216 + 0.03 , -3.96769 + 0.03, 2.65319 + 0.03 };
	WCSPH::Vector fluid2_size = { 1.96858 - 0.06, 3.70556, 0.15 };
	WCSPH::Vector fluid2_left = { 6.95855 + 0.03, -3.50449, 2.65319 + 0.03 };
	WCSPH::Vector fluid3_size = { 2.89497 - 0.06, 0.463194 - 0.06, 0.15 };
	WCSPH::Vector fluid3_left = { 6.03216 + 0.03, 0.20107 + 0.03, 2.65319 + 0.03 };

	std::vector<WCSPH::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid1_size);
	fluid_sizes.push_back(fluid2_size);
	fluid_sizes.push_back(fluid3_size);

	std::vector<WCSPH::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid1_left);
	fluid_lefts.push_back(fluid2_left);
	fluid_lefts.push_back(fluid3_left);

	SPHParameters params;
	MCParameters mcparams;
	params.dt = 0.0001;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.01;
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
	params.grid_cell_size = 1.2 * params.particle_radius;

	params.cohesion_coefficient = 0.1;

	params.emit_frequency = 1.2;
	params.max_num_particles = 1000000;

	params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
	mcparams.ours = false;
	mcparams.sparse = true;

	Emitter::Emitter emitter_left(0.05, { 1.0, 0.0, 0.0 }, 0.05, { 0.0, 1.0, 0.0 }, { 6.49536, -3.0413, 2.71108 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 3.0 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter_right(0.05, { 1.0, 0.0, 0.0 }, 0.05, { 0.0, 1.0, 0.0 }, { 6.49536, -0.262125, 2.71108 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 3.0 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter_mid(0.015, { 0.0, 1.0, 0.0 }, 0.015, { 1.0, 0.0, 1.0 }, { 6.16818, -1.65, 4.11 }, { 1.0, 0.0, -1.0 }, { 1.0, 0.0, -1.0 }, params.particle_diameter, params.particle_diameter);

	std::vector<CompactNSearch::Real> schedule;
	schedule.push_back(0.0);
	schedule.push_back(2.1);
	emitter_left.set_schedule(schedule);
	emitter_right.set_schedule(schedule);
	emitter_mid.set_schedule(schedule);
	std::vector<Emitter::Emitter> emitters;
	emitters.push_back(emitter_left);
	emitters.push_back(emitter_right);
	emitters.push_back(emitter_mid);

	WCSPH::SPH sph(false, false, true, "../res/italian_fountain_sph/italian_fountain_sph_", params, mcparams, emitters);
	sph.load_geometry(fluid_sizes, fluid_lefts);
	sph.simulate(1);
	sph.printStats("Italian Fountain SPH");
}
