#include "catch.hpp"
#include <iostream>
#include <chrono>
#include "../PBF/pbf.h"
#include "../geometry/io.h"


TEST_CASE("Dam break PBF", "[Dam break PBF]")
{
	std::cout << "Testing dam break in PBF" << std::endl;
	PBD::Vector boundary_size = { 0.18, 0.8, 1.0 };
	PBD::Vector boundary_left = { -0.015, -0.015, -0.015 };
	PBD::Vector fluid_size = { 0.15, 0.25, 0.5 };
	PBD::Vector fluid_left = { 0.0, 0.0, 0.0 };

	std::vector<PBD::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);

	std::vector<PBD::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);

	PBFParameters params;
	MCParameters mcparams;
	params.dt = 0.0005;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.005;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.001;
	params.fluid_pressure_stiffness = 1000.0;
	params.fluid_viscosity = 0.0025;
	params.boundary_viscosity = 0.0;

	params.particle_diameter = 2 * params.particle_radius;
	params.fluid_sampling_distance = params.particle_diameter;
	params.boundary_sampling_distance = 0.8 * params.particle_diameter;
	params.smoothing_length = 1.2 * params.particle_diameter;
	params.smoothing_length_squared = params.smoothing_length * params.smoothing_length;
	params.compact_support = 2 * params.smoothing_length;
	params.export_type = PBFParameters::export_type::EXPORT;

	PBD::PBF pbf(false, false, "../res/dam_break_pbf/dam_break_pbf_", params, mcparams);
	pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	pbf.simulate(1);
	pbf.printStats("Dam break PBF");
}

/*TEST_CASE("Crazy stuff", "[Crazy stuff]")
{
	std::cout << "Testing crazy stuff" << std::endl;
	PBD::Vector boundary_size = { 2, 1, 2 };
	PBD::Vector boundary_left = { 0.0, 0.0, 0.0 };
	PBD::Vector fluid_size = { 0.5, 0.5, 0.5 };
	PBD::Vector fluid_left = { 0.01, 0.01, 1 };
	PBD::Vector fluid_size2 = { 0.5, 0.5, 0.5 };
	PBD::Vector fluid_left2 = { 1, 0.01, 1 };

	std::vector<PBD::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);
	fluid_sizes.push_back(fluid_size2);

	std::vector<PBD::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);
	fluid_lefts.push_back(fluid_left2);

	PBFParameters params;
	MCParameters mcparams;
	params.dt = 0.0001;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.005;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.00025;
	params.fluid_pressure_stiffness = 1000.0;
	params.fluid_viscosity = 0.0025;
	params.boundary_viscosity = 0.0;

	params.particle_diameter = 2 * params.particle_radius;
	params.fluid_sampling_distance = params.particle_diameter;
	params.boundary_sampling_distance = 0.8 * params.particle_diameter;
	params.smoothing_length = 1.2 * params.particle_diameter;
	params.smoothing_length_squared = params.smoothing_length * params.smoothing_length;
	params.compact_support = 2 * params.smoothing_length;

	PBD::PBF pbf(false, false, "../res/crazy_stuff_pbf/crazy_stuff_pbf_", params, mcparams);
	pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	pbf.simulate(0.1);
	pbf.printStats("Crazy stuff");
}

TEST_CASE("Crazy stuff2", "[Crazy stuff2]")
{
	std::cout << "Testing crazy stuff" << std::endl;
	PBD::Vector boundary_size = { 2, 1, 2 };
	PBD::Vector boundary_left = { 0.0, 0.0, 0.0 };
	PBD::Vector fluid_size = { 0.1, 0.1, 0.1 };
	PBD::Vector fluid_left = { 0.015, 0.015, 1 };
	PBD::Vector fluid_size2 = { 0.1, 0.1, 0.1 };
	PBD::Vector fluid_left2 = { 1, 0.015, 1 };

	std::vector<PBD::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);
	fluid_sizes.push_back(fluid_size2);

	std::vector<PBD::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);
	fluid_lefts.push_back(fluid_left2);

	PBFParameters params;
	MCParameters mcparams;
	params.dt = 0.0001;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.005;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.00025;
	params.fluid_pressure_stiffness = 1000.0;
	params.fluid_viscosity = 0.0025;
	params.boundary_viscosity = 0.0;

	params.particle_diameter = 2 * params.particle_radius;
	params.fluid_sampling_distance = params.particle_diameter;
	params.boundary_sampling_distance = 0.8 * params.particle_diameter;
	params.smoothing_length = 1.2 * params.particle_diameter;
	params.smoothing_length_squared = params.smoothing_length * params.smoothing_length;
	params.compact_support = 2 * params.smoothing_length;

	PBD::PBF pbf(false, false, "../res/crazy_stuff2_pbf/crazy_stuff2_pbf_", params, mcparams);
	pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	pbf.simulate(0.1);
	pbf.printStats("Crazy stuff 2");
}*/