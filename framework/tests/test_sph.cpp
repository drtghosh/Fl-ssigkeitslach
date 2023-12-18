#include "catch.hpp"
#include <iostream>
#include <chrono>
#include "../SPH/sph.h"
#include "../geometry/io.h"



// Check out https://github.com/catchorg/Catch2 for more information about how to use Catch2
/*TEST_CASE("Gravity only", "[Gravity only]")
{
	std::cout << "Testing gravity only" << std::endl;
	//const std::vector<learnSPH::TriMesh> meshes = learnSPH::read_tri_meshes_from_obj("../res/box.obj");
	//const learnSPH::TriMesh& box = meshes[0]; 
	WCSPH::Vector boundary_size = { 1.0, 1.0, 10.0 };
	WCSPH::Vector boundary_left = { 0.0, 0.0, 0.0 };
	WCSPH::Vector fluid_size = { 1.0, 1.0, 1.0 };
	WCSPH::Vector fluid_left = { 0.0, 0.0, 9.0 };

	std::vector<WCSPH::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);

	std::vector<WCSPH::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);

	SPHParameters params;
	MCParameters mcparams;
	params.dt = 0.01;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.03;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.1;

	WCSPH::SPH sph(true, false, "../res/gravity/gravity_", params, mcparams);
	sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	sph.simulate(5);
	sph.printStats("Gravity only");
}

TEST_CASE("Incremental SPH", "[Incremental SPH]")
{
	std::cout << "Incremental tests of SPH:" << std::endl;
	//const std::vector<learnSPH::TriMesh> meshes = learnSPH::read_tri_meshes_from_obj("../res/box.obj");
	//const learnSPH::TriMesh& box = meshes[0]; 
	WCSPH::Vector boundary_size;
	WCSPH::Vector boundary_left;
	WCSPH::Vector fluid_size = { 1.0, 1.0, 1.0 };
	WCSPH::Vector fluid_left = { 0.0, 0.0, 0.0 };

	std::vector<WCSPH::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);

	std::vector<WCSPH::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);


	SPHParameters params;
	MCParameters mcparams;
	params.dt = 0.001;
	params.dt_next_frame = 0.1;
	params.particle_radius = 0.03;
	params.fluid_rest_density = 1000;
	params.fluid_pressure_stiffness = 0.0;
	params.fluid_viscosity = 0.0;
	params.boundary_viscosity = 0.0;
	params.max_dt = 0.001;

	SECTION("No Boundary") {
		WCSPH::SPH sph(false, false, "../res/no_bounds/no_bounds_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(1);
		sph.printStats("No boundary");
	}

	params.fluid_pressure_stiffness = 1000.0;
	SECTION("Added stiffness") {
		WCSPH::SPH sph(false, false, "../res/no_bounds_fluid_stiffness/no_bounds_fluid_stiffness_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(1);
		sph.printStats("No boundary + stiffness");
	}

	params.fluid_viscosity = 0.1;
	SECTION("Added fluid viscosity") {
		WCSPH::SPH sph(false, false, "../res/no_bounds_fluid_visco/no_bounds_fluid_visco_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(1);
		sph.printStats("No boundary + stiffness + viscosity");
	}

	SECTION("Added another fluid box") {
		WCSPH::Vector second_fluid_size = { 1.0, 1.0, 1.0 };
		WCSPH::Vector second_fluid_left = { 1.0, 0.0, 0.0 };
		fluid_sizes.push_back(second_fluid_size);
		fluid_lefts.push_back(second_fluid_left);

		WCSPH::SPH sph(false, true, "../res/no_bounds_two_fluids/no_bounds_two_fluids_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(1);
		sph.printStats("Two fluid boxes");
	}

	SECTION("Turned on gravity") {
		WCSPH::SPH sph(false, false, "../res/no_bounds_gravity/no_bounds_gravity_", params, mcparams);
		sph.turn_on_gravity();
		sph.load_geometry(false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(1);
		sph.printStats("No boundary + gravity");
	}

	SECTION("Added a floor") {
		WCSPH::Vector boundary_size = { 21.0, 21.0, 11.0 };
		WCSPH::Vector boundary_left = { -10.0, -10.0, -1.0 };

		WCSPH::SPH sph(false, false, "../res/bounds_floor/bounds_floor_", params, mcparams);
		sph.turn_on_gravity();
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(1);
		sph.printStats("Floor + gravity");
	}

	params.boundary_viscosity = 1.0;
	SECTION("Added boundary viscosity") {
		WCSPH::Vector boundary_size = { 21.0, 21.0, 11.0 };
		WCSPH::Vector boundary_left = { -10.0, -10.0, -1.0 };

		WCSPH::SPH sph(false, false, "../res/bounds_floor_boundary_visco/bounds_floor_boundary_visco_", params, mcparams);
		sph.turn_on_gravity();
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(1);
		sph.printStats("Floor + gravity + boundary viscosity");
	}
}

TEST_CASE("Dam break SPH", "[Dam break SPH]")
{
	std::cout << "Testing dam break in SPH" << std::endl;
	WCSPH::Vector boundary_size = { 0.18, 0.8, 1.0 };
	WCSPH::Vector boundary_left = { -0.015, -0.015, -0.015 };
	WCSPH::Vector fluid_size = { 0.15, 0.25, 0.5 };
	WCSPH::Vector fluid_left = { 0.0, 0.0, 0.0 };

	std::vector<WCSPH::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);

	std::vector<WCSPH::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);

	SPHParameters params;
	MCParameters mcparams;
	params.dt = 0.00025;
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

	WCSPH::SPH sph(false, false, "../res/dam_break_sph/dam_break_sph_", params, mcparams);
	sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	sph.simulate(4);
	sph.printStats("Dam break SPH");
}

TEST_CASE("Crazy stuff", "[Crazy stuff]")
{
	std::cout << "Testing crazy stuff" << std::endl;
	WCSPH::Vector boundary_size = { 2, 1, 2 };
	WCSPH::Vector boundary_left = { 0.0, 0.0, 0.0 };
	WCSPH::Vector fluid_size = { 0.5, 0.5, 0.5 };
	WCSPH::Vector fluid_left = { 0.01, 0.01, 1 };
	WCSPH::Vector fluid_size2 = { 0.5, 0.5, 0.5 };
	WCSPH::Vector fluid_left2 = { 1, 0.01, 1 };

	std::vector<WCSPH::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);
	fluid_sizes.push_back(fluid_size2);

	std::vector<WCSPH::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);
	fluid_lefts.push_back(fluid_left2);

	SPHParameters params;
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

	WCSPH::SPH sph(false, false, "../res/crazy_stuff_sph/crazy_stuff_sph_", params, mcparams);
	sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	sph.simulate(0.1);
	sph.printStats("Crazy stuff");
}

TEST_CASE("Crazy stuff2", "[Crazy stuff2]")
{
	std::cout << "Testing crazy stuff" << std::endl;
	WCSPH::Vector boundary_size = { 2, 1, 2 };
	WCSPH::Vector boundary_left = { 0.0, 0.0, 0.0 };
	WCSPH::Vector fluid_size = { 0.1, 0.1, 0.1 };
	WCSPH::Vector fluid_left = { 0.015, 0.015, 1 };
	WCSPH::Vector fluid_size2 = { 0.1, 0.1, 0.1 };
	WCSPH::Vector fluid_left2 = { 1, 0.015, 1 };

	std::vector<WCSPH::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);
	fluid_sizes.push_back(fluid_size2);

	std::vector<WCSPH::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);
	fluid_lefts.push_back(fluid_left2);

	SPHParameters params;
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

	WCSPH::SPH sph(false, false, "../res/crazy_stuff2_sph/crazy_stuff2_sph_", params, mcparams);
	sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	sph.simulate(0.1);
	sph.printStats("Crazy stuff 2");
}*/