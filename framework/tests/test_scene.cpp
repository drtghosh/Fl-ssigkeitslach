#include "catch.hpp"
#include "../geometry/emitter.h"
#include "../geometry/io.h"
#include "../SPH/sph.h"
#include "../SPH/params.h"
#include "../PBF/pbf.h"
#include "../PBF/params.h"
#include <cmath>
#include<cstdlib>

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

/*TEST_CASE("Italian Fountain SPH", "[Italian Fountain SPH]")
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
	sph.load_geometry_from_obj(fluid_sizes, fluid_lefts, "../res/ItalianFountain.obj");
	sph.simulate(1);
	sph.printStats("Italian Fountain SPH");
}

TEST_CASE("Lotus Leave SPH", "[Lotus Leave SPH]")
{
	std::cout << "Testing Lotus Leave in SPH" << std::endl;
	WCSPH::Vector fluid1_size = { 0.5, 0.5, 0.5 };
	WCSPH::Vector fluid1_left = { -1.84 , 3.0, -3.1 };
	WCSPH::Vector fluid2_size = { 0.25, 0.25, 0.25 };
	WCSPH::Vector fluid2_left = { -1.7 , 3.12, -2.3 };

	std::vector<WCSPH::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid1_size);
	fluid_sizes.push_back(fluid2_size);

	std::vector<WCSPH::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid1_left);
	fluid_lefts.push_back(fluid2_left);

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

	params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
	mcparams.ours = false;
	mcparams.sparse = true;

	WCSPH::SPH sph(false, false, true, "../res/lotus_leave_sph/lotus_leave_sph_", params, mcparams);
	sph.load_geometry_from_obj(fluid_sizes, fluid_lefts, "../res/lotus.obj");
	sph.simulate(1.2);
	sph.printStats("Lotus Leave SPH");
}

TEST_CASE("Adhesive Lotus Leave SPH", "[Adhesive Lotus Leave SPH]")
{
	std::cout << "Testing Adhesive Lotus Leave in SPH" << std::endl;
	WCSPH::Vector fluid1_size = { 0.5, 0.5, 0.5 };
	WCSPH::Vector fluid1_left = { -1.84 , 3.0, -3.1 };
	WCSPH::Vector fluid2_size = { 0.25, 0.25, 0.25 };
	WCSPH::Vector fluid2_left = { -1.7 , 3.12, -2.3 };

	std::vector<WCSPH::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid1_size);
	fluid_sizes.push_back(fluid2_size);

	std::vector<WCSPH::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid1_left);
	fluid_lefts.push_back(fluid2_left);

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
	params.adhesion_coefficient = 2.0;

	params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
	mcparams.ours = false;
	mcparams.sparse = true;

	WCSPH::SPH sph(false, false, true, "../res/adhesive_lotus_leave_sph/adhesive_lotus_leave_sph_", params, mcparams);
	sph.load_geometry_from_obj(fluid_sizes, fluid_lefts, "../res/lotus.obj");
	sph.simulate(0.9);
	sph.printStats("Adhesive Lotus Leave SPH");
}

TEST_CASE("No Surface Tension Lotus Leave SPH", "[No Surface Tension Lotus Leave SPH]")
{
	std::cout << "Testing No Surface Tension Lotus Leave in SPH" << std::endl;
	WCSPH::Vector fluid1_size = { 0.5, 0.5, 0.5 };
	WCSPH::Vector fluid1_left = { -1.84 , 3.0, -3.1 };
	WCSPH::Vector fluid2_size = { 0.25, 0.25, 0.25 };
	WCSPH::Vector fluid2_left = { -1.7 , 3.12, -2.3 };

	std::vector<WCSPH::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid1_size);
	fluid_sizes.push_back(fluid2_size);

	std::vector<WCSPH::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid1_left);
	fluid_lefts.push_back(fluid2_left);

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

	params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
	mcparams.ours = false;
	mcparams.sparse = true;

	WCSPH::SPH sph(false, false, false, "../res/no_st_lotus_leave_sph/no_st_lotus_leave_sph_", params, mcparams);
	sph.load_geometry_from_obj(fluid_sizes, fluid_lefts, "../res/lotus.obj");
	sph.simulate(0.9);
	sph.printStats("No Surface Tension Lotus Leave SPH");
}

TEST_CASE("Lotus Leave Rain SPH", "[Lotus Leave Rain SPH]")
{
	std::cout << "Testing Lotus Leave Rain in SPH" << std::endl;

	std::vector<WCSPH::Vector> fluid_sizes;
	std::vector<WCSPH::Vector> fluid_lefts;
	WCSPH::Vector fluid_size = { 0.2, 0.2, 0.2 };
	WCSPH::Vector fluid_left1 = { -3.65 , 1.0, -3.12 };
	WCSPH::Vector fluid_left2 = { -2.05 , 0.90, -3.12 };
	WCSPH::Vector fluid_left3 = { -3.5 , 3.15, -3.12 };
	WCSPH::Vector fluid_left4 = { -0.5 , 1.1, -3.12 };

	srand((unsigned)time(NULL));
	for (int i = 2; i < 8; i++) {
		fluid_sizes.push_back(fluid_size);
		fluid_lefts.push_back(fluid_left1 + i * WCSPH::Vector({ 0.5, 0.5, 0.0 }) + WCSPH::Vector({ (rand() % 10) / 50.0, (rand() % 10) / 100.0, 0.0 }));
	}

	for (int i = 0; i < 6; i++) {
		fluid_sizes.push_back(fluid_size);
		fluid_lefts.push_back(fluid_left2 + i * WCSPH::Vector({ 0.5, 0.5, 0.0 }) + WCSPH::Vector({ -(rand() % 10) / 40.0, (rand() % 10) / 50.0, 0.0 }));
	}

	for (int i = 0; i < 4; i++) {
		fluid_sizes.push_back(fluid_size);
		int random = (rand() % 10) / 10;
		fluid_lefts.push_back(fluid_left3 + i * WCSPH::Vector({ 0.5, 0.5, 0.0 }) + WCSPH::Vector({ -(rand() % 10) / 50.0, -(rand() % 10) / 50.0, 0.0 }));
	}

	for (int i = 0; i < 2; i++) {
		fluid_sizes.push_back(fluid_size);
		fluid_lefts.push_back(fluid_left4 + i * WCSPH::Vector({ 0.5, 0.5, 0.0 }));
	}

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

	//params.cohesion_coefficient = 0.1;

	params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
	mcparams.ours = false;
	mcparams.sparse = true;

	WCSPH::SPH sph(false, false, true, "../res/lotus_leave_rain_sph/lotus_leave_rain_sph_", params, mcparams);
	sph.load_geometry_from_obj(fluid_sizes, fluid_lefts, "../res/lotus.obj");
	sph.simulate(0.3);
	sph.printStats("Lotus Leave Rain SPH");
}

TEST_CASE("Italian Fountain PBF", "[Italian Fountain PBF]")
{
	std::cout << "Testing Italian Fountain in PBF" << std::endl;

	CompactNSearch::Real particle_radius = 0.01;

	PBD::Vector fluid1_size = { 2.89497 - (6 * particle_radius), 0.4632 - (6 * particle_radius), 0.065};
	PBD::Vector fluid1_left = { 6.03216 + (3 * particle_radius) , -3.96769 + (3 * particle_radius), 2.65319 + (3 * particle_radius) };
	PBD::Vector fluid2_size = { 1.96858 - (6 * particle_radius), 3.70556, 0.065 };
	PBD::Vector fluid2_left = { 6.95855 + (3 * particle_radius), -3.50449, 2.65319 + (3 * particle_radius) };
	PBD::Vector fluid3_size = { 2.89497 - (6 * particle_radius), 0.463194 - (6 * particle_radius), 0.065 };
	PBD::Vector fluid3_left = { 6.03216 + (3 * particle_radius), 0.20107 + (3 * particle_radius), 2.65319 + (3 * particle_radius) };

	std::vector<PBD::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid1_size);
	fluid_sizes.push_back(fluid2_size);
	fluid_sizes.push_back(fluid3_size);

	std::vector<PBD::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid1_left);
	fluid_lefts.push_back(fluid2_left);
	fluid_lefts.push_back(fluid3_left);

	PBFParameters params;
	MCParameters mcparams;
	params.dt = 0.001;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.01;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.002;
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

	params.pbf_iterations = 10;

	params.grid_cell_size = 1.2 * params.particle_radius;

	params.cohesion_coefficient = 0.01;

	params.emit_frequency = 1.2;
	params.max_num_particles = 1000000;

	params.export_type = PBFParameters::EXPORT_WITH_SURFACE;
	mcparams.ours = false;
	mcparams.sparse = true;

	Emitter::Emitter emitter_left(0.05, { 1.0, 0.0, 0.0 }, 0.05, { 0.0, 1.0, 0.0 }, { 6.49536, -3.0413, 2.71108 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 5.0 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter_right(0.05, { 1.0, 0.0, 0.0 }, 0.05, { 0.0, 1.0, 0.0 }, { 6.49536, -0.262125, 2.71108 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 5.0 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter1_mid(0.015, { 0.0, 1.0, 0.0 }, 0.015, { 1.0, 0.0, 1.0 }, { 6.16818, -1.656, 4.11 }, { 1.0, 0.0, -1.0 }, { 1.0, 0.0, -1.0 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter2_mid(0.015, { 0.0, 1.0, 0.0 }, 0.015, { 1.0, 0.0, 1.0 }, { 6.16818, -1.656, 4.11 }, { 1.0, 0.0, -1.0 }, { 2.0, 0.0, -2.0 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter3_mid(0.015, { 0.0, 1.0, 0.0 }, 0.015, { 1.0, 0.0, 1.0 }, { 6.16818, -1.656, 4.11 }, { 1.0, 0.0, -1.0 }, { 0.0, 0.0, -0.5 }, params.particle_diameter, params.particle_diameter);
	PBD::Vector f1r1 = { -0.2957, -0.40005, -0.05022 };
	PBD::Vector f1r2 = { 0.14812, -0.1657095858017, 0.44789 };
	PBD::Vector f1n = f1r1.cross(f1r2);
	std::cout << "f1n: " << f1n << std::endl;
	Emitter::Emitter emitter_left_frog(0.015, f1r1, 0.015, f1r2, { 8.9213, -3.98846, 3.3062 }, f1n, 2.5 * f1n, params.particle_diameter, params.particle_diameter);
	PBD::Vector f2r1 = { 0.26112, -0.424693, -0.03812 };
	PBD::Vector f2r2 = { 0.2029537969822, 0.084476, 0.44908 };
	PBD::Vector f2n = f2r1.cross(f2r2);
	std::cout << "f2n: " << f2n << std::endl;
	Emitter::Emitter emitter_right_frog(0.015, f2r1, 0.015, f2r2, { 8.9413, 0.64, 3.3062 }, f2n, 2.5 * f2n, params.particle_diameter, params.particle_diameter);

	std::vector<CompactNSearch::Real> schedule_fountain;
	schedule_fountain.push_back(0.0);
	schedule_fountain.push_back(10.0);
	emitter_left.set_schedule(schedule_fountain);
	emitter_right.set_schedule(schedule_fountain);

	std::vector<CompactNSearch::Real> schedule1_tap;
	schedule1_tap.push_back(0.0);
	schedule1_tap.push_back(3.0);
	emitter1_mid.set_schedule(schedule1_tap);
	std::vector<CompactNSearch::Real> schedule2_tap;
	schedule2_tap.push_back(3.0);
	schedule2_tap.push_back(5.0);
	emitter2_mid.set_schedule(schedule2_tap);
	std::vector<CompactNSearch::Real> schedule3_tap;
	schedule3_tap.push_back(5.0);
	schedule3_tap.push_back(8.0);
	emitter3_mid.set_schedule(schedule3_tap);

	std::vector<CompactNSearch::Real> schedule_frogs;
	schedule_frogs.push_back(0.0);
	schedule_frogs.push_back(2.0);
	schedule_frogs.push_back(4.0);
	schedule_frogs.push_back(6.0);
	schedule_frogs.push_back(8.0);
	schedule_frogs.push_back(10.0);
	emitter_left_frog.set_schedule(schedule_frogs);
	emitter_right_frog.set_schedule(schedule_frogs);
	std::vector<Emitter::Emitter> emitters;
	emitters.push_back(emitter_left);
	emitters.push_back(emitter_right);
	emitters.push_back(emitter1_mid);
	emitters.push_back(emitter2_mid);
	emitters.push_back(emitter3_mid);
	emitters.push_back(emitter_left_frog);
	emitters.push_back(emitter_right_frog);

	PBD::PBF pbf(false, false, true, "../res/italian_fountain_pbf/italian_fountain_pbf_", params, mcparams, emitters);
	pbf.load_geometry_from_obj(fluid_sizes, fluid_lefts, "../res/ItalianFountain2.obj");
	pbf.simulate(10);
	pbf.printStats("Italian Fountain PBF");
}*/

TEST_CASE("Italian Fountain PBF", "[Italian Fountain PBF]")
{
	std::cout << "Testing Italian Fountain in PBF" << std::endl;

	CompactNSearch::Real particle_radius = 0.01;

	PBD::Vector fluid1_size = { (8.34947 - 6.03216) - (6 * particle_radius), (-2.69144 + 3.96769) - (6 * particle_radius), 0.075 };
	PBD::Vector fluid1_left = { 6.03216 + (3 * particle_radius) , -3.96769 + (3 * particle_radius), 2.65319 + (3 * particle_radius) };
	PBD::Vector fluid2_size = { (8.34947 - 6.95855) - (6 * particle_radius), (-0.607062 + 2.69144), 0.075 };
	PBD::Vector fluid2_left = { 6.95855 + (3 * particle_radius), -2.69144, 2.65319 + (3 * particle_radius) };
	PBD::Vector fluid3_size = { (8.34947 - 6.03216) - (6 * particle_radius), (0.664264 + 0.607062) - (6 * particle_radius), 0.075 };
	PBD::Vector fluid3_left = { 6.03216 + (3 * particle_radius), -0.607062 + (3 * particle_radius), 2.65319 + (3 * particle_radius) };
	PBD::Vector fluid4_size = { (6.15775 - 6.03216) - (6 * particle_radius), (-1.54176 + 1.76554), 0.0375 };
	PBD::Vector fluid4_left = { 6.03216 + (3 * particle_radius), -1.76554, 3.71416 + (3 * particle_radius) };
	PBD::Vector fluid5_size = { (6.16656 - 6.03216) - (6 * particle_radius), (-1.51732 + 1.7861), 0.0375 };
	PBD::Vector fluid5_left = { 6.03216 + (3 * particle_radius), -1.7861, 3.21407 + (3 * particle_radius) };
	PBD::Vector fluid6_size = { 0.25, (-1.16533 + 2.13809), 0.075 };
	PBD::Vector fluid6_left = { 6.05241 + (3 * particle_radius), -2.13809, 2.73681 + (3 * particle_radius) };

	std::vector<PBD::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid1_size);
	fluid_sizes.push_back(fluid2_size);
	fluid_sizes.push_back(fluid3_size);
	fluid_sizes.push_back(fluid4_size);
	fluid_sizes.push_back(fluid5_size);
	fluid_sizes.push_back(fluid6_size);

	std::vector<PBD::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid1_left);
	fluid_lefts.push_back(fluid2_left);
	fluid_lefts.push_back(fluid3_left);
	fluid_lefts.push_back(fluid4_left);
	fluid_lefts.push_back(fluid5_left);
	fluid_lefts.push_back(fluid6_left);

	PBFParameters params;
	MCParameters mcparams;
	params.dt = 0.001;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.01;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.002;
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

	params.pbf_iterations = 10;

	params.grid_cell_size = 1.2 * params.particle_radius;

	params.cohesion_coefficient = 0.01;

	params.emit_frequency = 1.2;
	params.max_num_particles = 1500000;

	params.export_type = PBFParameters::EXPORT_WITH_SURFACE;
	mcparams.ours = false;
	mcparams.sparse = true;

	Emitter::Emitter emitter_left(0.05, { 1.0, 0.0, 0.0 }, 0.05, { 0.0, 1.0, 0.0 }, { 6.49536, -3.0413, 2.65319 + (7 * particle_radius) + 0.075 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 5.0 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter_right(0.05, { 1.0, 0.0, 0.0 }, 0.05, { 0.0, 1.0, 0.0 }, { 6.49536, -0.262125, 2.65319 + (7 * particle_radius) + 0.075 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 5.0 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter_mid(0.015, { 0.0, 1.0, 0.0 }, 0.015, { 1.0, 0.0, 1.0 }, { 6.16818, -1.656, 4.11 }, { 1.0, 0.0, -1.0 }, { 0.0, 0.0, -1.0 }, params.particle_diameter, params.particle_diameter);
	PBD::Vector f1r1 = { -0.2957, -0.40005, -0.05022 };
	PBD::Vector f1r2 = { 0.14812, -0.1657095858017, 0.44789 };
	PBD::Vector f1n = f1r1.cross(f1r2);
	std::cout << "f1n: " << f1n << std::endl;
	Emitter::Emitter emitter_left_frog(0.02, f1r1, 0.02, f1r2, { 8.20923, -3.89947, 3.45903 }, f1n, 5 * f1n, params.particle_diameter, params.particle_diameter);
	PBD::Vector f2r1 = { 0.26112, -0.424693, -0.03812 };
	PBD::Vector f2r2 = { 0.2029537969822, 0.084476, 0.44908 };
	PBD::Vector f2n = f2r1.cross(f2r2);
	std::cout << "f2n: " << f2n << std::endl;
	Emitter::Emitter emitter_right_frog(0.02, f2r1, 0.02, f2r2, { 8.25884, 0.541795, 3.45823 }, f2n, 5 * f2n, params.particle_diameter, params.particle_diameter);

	std::vector<CompactNSearch::Real> schedule_fountain;
	schedule_fountain.push_back(0.0);
	schedule_fountain.push_back(4.0);
	schedule_fountain.push_back(6.0);
	schedule_fountain.push_back(10.0);
	schedule_fountain.push_back(12.0);
	schedule_fountain.push_back(16.0);
	emitter_left.set_schedule(schedule_fountain);
	emitter_right.set_schedule(schedule_fountain);

	std::vector<CompactNSearch::Real> schedule_tap;
	schedule_tap.push_back(0.0);
	schedule_tap.push_back(20.0);
	emitter_mid.set_schedule(schedule_tap);

	std::vector<CompactNSearch::Real> schedule_frogs;
	schedule_frogs.push_back(0.0);
	schedule_frogs.push_back(2.0);
	schedule_frogs.push_back(4.0);
	schedule_frogs.push_back(6.0);
	schedule_frogs.push_back(8.0);
	schedule_frogs.push_back(10.0);
	schedule_frogs.push_back(12.0);
	schedule_frogs.push_back(14.0);
	emitter_left_frog.set_schedule(schedule_frogs);
	emitter_right_frog.set_schedule(schedule_frogs);
	std::vector<Emitter::Emitter> emitters;
	emitters.push_back(emitter_left);
	emitters.push_back(emitter_right);
	emitters.push_back(emitter_mid);
	emitters.push_back(emitter_left_frog);
	emitters.push_back(emitter_right_frog);

	PBD::PBF pbf(false, false, true, "../res/italian_fountain_pbf/italian_fountain_pbf_", params, mcparams, emitters);
	pbf.load_geometry_from_obj(fluid_sizes, fluid_lefts, "../res/ItalianFountain3.obj");
	pbf.simulate(0.01);
	pbf.printStats("Italian Fountain PBF");
}

/*TEST_CASE("Lotus Leave PBF", "[Lotus Leave PBF]")
{
	std::cout << "Testing Lotus Leave in PBF" << std::endl;

	CompactNSearch::Real particle_radius = 0.01;
	PBD::Vector fluid1_size = { 0.5, 0.5, 0.5 };
	PBD::Vector fluid1_left = { -1.84 , 3.0, -3.1 };
	PBD::Vector fluid2_size = { 0.25, 0.25, 0.25 };
	PBD::Vector fluid2_left = { -1.7 , 3.12, -2.3 };
	//PBD::Vector fluid3_size = { 3.0 - (6 * particle_radius), 3.0 - (6 * particle_radius), 0.013 };
	//PBD::Vector fluid3_left = { -1.55 + (3 * particle_radius), 3.5 + (3 * particle_radius), -4.86 + (3 * particle_radius) };

	std::vector<PBD::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid1_size);
	fluid_sizes.push_back(fluid2_size);
	//fluid_sizes.push_back(fluid3_size);

	std::vector<PBD::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid1_left);
	fluid_lefts.push_back(fluid2_left);
	//fluid_lefts.push_back(fluid3_left);

	PBFParameters params;
	MCParameters mcparams;
	params.dt = 0.001;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.01;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.002;
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

	params.pbf_iterations = 10;

	params.export_type = PBFParameters::EXPORT_WITH_SURFACE;
	mcparams.ours = false;
	mcparams.sparse = true;

	PBD::PBF pbf(false, false, true, "../res/lotus_leave_pbf/lotus_leave_pbf_", params, mcparams);
	pbf.load_geometry_from_obj(fluid_sizes, fluid_lefts, "../res/lotus.obj");
	pbf.simulate(5.0);
	pbf.printStats("Lotus Leave PBF");
}

TEST_CASE("Adhesive Lotus Leave PBF", "[Adhesive Lotus Leave PBF]")
{
	std::cout << "Testing Adhesive Lotus Leave in PBF" << std::endl;

	CompactNSearch::Real particle_radius = 0.01;
	PBD::Vector fluid1_size = { 0.5, 0.5, 0.5 };
	PBD::Vector fluid1_left = { -1.84 , 3.0, -3.1 };
	PBD::Vector fluid2_size = { 0.25, 0.25, 0.25 };
	PBD::Vector fluid2_left = { -1.7 , 3.12, -2.3 };
	//PBD::Vector fluid3_size = { 10.0 - (6 * particle_radius), 10.0 - (6 * particle_radius), 0.013 };
	//PBD::Vector fluid3_left = { -5.05 + (3 * particle_radius), 0.0 + (3 * particle_radius), -4.86 + (3 * particle_radius) };

	std::vector<PBD::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid1_size);
	fluid_sizes.push_back(fluid2_size);
	//fluid_sizes.push_back(fluid3_size);

	std::vector<PBD::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid1_left);
	fluid_lefts.push_back(fluid2_left);
	//fluid_lefts.push_back(fluid3_left);

	PBFParameters params;
	MCParameters mcparams;
	params.dt = 0.001;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.01;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.002;
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

	params.pbf_iterations = 10;

	params.cohesion_coefficient = 0.1;
	params.adhesion_coefficient = 2.0;

	params.export_type = PBFParameters::EXPORT_WITH_SURFACE;
	mcparams.ours = false;
	mcparams.sparse = true;

	PBD::PBF pbf(false, false, true, "../res/adhesive_lotus_leave_pbf/adhesive_lotus_leave_pbf_", params, mcparams);
	pbf.load_geometry_from_obj(fluid_sizes, fluid_lefts, "../res/lotus.obj");
	pbf.simulate(5.0);
	pbf.printStats("Adhesive Lotus Leave PBF");
}*/