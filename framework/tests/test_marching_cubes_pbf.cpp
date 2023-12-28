#include "catch.hpp"
#include <iostream>
#include <chrono>
#include "../mcubes/marching_cubes.h"
#include "CompactNSearch/Config.h"
#include "../PBF/config.h"
#include "../geometry/io.h"
#include "../PBF/pbf.h"
#include "../PBF/params.h"
#include "../mcubes/mcparams.h"

// Check out https://github.com/catchorg/Catch2 for more information about how to use Catch2

/*TEST_CASE("Dam Break PBF with Surface Construction", "[Dam Break PBF with Surface Construction]")
{
	std::cout << "Testing dam break in PBF with surface construction" << std::endl;

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
	
	params.grid_cell_size = 1.0 * params.particle_radius;

	params.export_type = PBFParameters::export_type::EXPORT_WITH_SURFACE;

	SECTION("Regular Marching Cubes") {
		mcparams.ours = false;
		mcparams.sparse = false;
		PBD::PBF pbf(false, false, "../res/dam_break_pbf_mc/dam_break_pbf_mc_", params, mcparams);
		pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(4);
		pbf.printStats("Dam break in PBF with regular marching cubes");
	}
	SECTION("Sparse Marching Cubes") {
		mcparams.sparse = true;
		PBD::PBF pbf(false, false, "../res/dam_break_pbf_mc_sparse/dam_break_pbf_mc_sparse_", params, mcparams);
		pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(4);
		pbf.printStats("Dam break in PBF with sparse marching cubes");
	}
	SECTION("Our Marching Cubes") {
		mcparams.sparse = false;
		mcparams.ours = true;
		PBD::PBF pbf(false, false, "../res/dam_break_pbf_mc_ours/dam_break_pbf_mc_ours_", params, mcparams);
		pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(1);
		pbf.printStats("Dam break in PBF with our marching cubes");
	}
}

TEST_CASE("Dam Break with obstacles PBF with Surface Construction", "[Dam Break with obstacles PBF with Surface Construction]")
{
	std::cout << "Testing dam break with obstacles in PBF with surface construction" << std::endl;

	PBD::Vector boundary_size = { 0.18, 0.8, 1.0 };
	PBD::Vector boundary_left = { -0.015, -0.015, -0.015 };
	PBD::Vector fluid_size = { 0.15, 0.25, 0.5 };
	PBD::Vector fluid_left = { 0.0, 0.0, 0.0 };

	std::vector<PBD::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);

	std::vector<PBD::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);

	std::vector<std::array<PBD::Vector, 4>> obstacle_squares;
	std::array<PBD::Vector, 4> square1{ PBD::Vector{0.055,0.5,-0.015}, PBD::Vector{0.085,0.5,-0.015}, PBD::Vector{0.085,0.5,0.5}, PBD::Vector{0.055,0.5,0.5} };
	std::array<PBD::Vector, 4> square2{ PBD::Vector{0.055,0.6,-0.015}, PBD::Vector{0.085,0.6,-0.015}, PBD::Vector{0.085,0.6,0.5}, PBD::Vector{0.055,0.6,0.5} };
	std::array<PBD::Vector, 4> square3{ PBD::Vector{0.055,0.5,0.5}, PBD::Vector{0.085,0.5,0.5}, PBD::Vector{0.085,0.6,0.5}, PBD::Vector{0.055,0.6,0.5} };
	std::array<PBD::Vector, 4> square4{ PBD::Vector{0.055,0.5,-0.015}, PBD::Vector{0.055,0.5,0.5}, PBD::Vector{0.055,0.6,0.5}, PBD::Vector{0.055,0.6,-0.015} };
	std::array<PBD::Vector, 4> square5{ PBD::Vector{0.085,0.5,-0.015}, PBD::Vector{0.085,0.5,0.5}, PBD::Vector{0.085,0.6,0.5}, PBD::Vector{0.085,0.6,-0.015} };
	obstacle_squares.push_back(square1);
	obstacle_squares.push_back(square2);
	obstacle_squares.push_back(square3);
	obstacle_squares.push_back(square4);
	obstacle_squares.push_back(square5);

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

	params.grid_cell_size = 1.0 * params.particle_radius;

	params.export_type = PBFParameters::export_type::EXPORT_WITH_SURFACE;

	SECTION("Regular Marching Cubes") {
		mcparams.ours = false;
		mcparams.sparse = false;
		PBD::PBF pbf(false, false, "../res/dam_break_obstacle_pbf_mc/dam_break_obstacle_pbf_mc_", params, mcparams);
		pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts, obstacle_squares);
		pbf.simulate(4);
		pbf.printStats("Dam break with obstacles in PBF with regular marching cubes");
	}
	SECTION("Sparse Marching Cubes") {
		mcparams.sparse = true;
		PBD::PBF pbf(false, false, "../res/dam_break_obstacle_pbf_mc_sparse/dam_break_obstacle_pbf_mc_sparse_", params, mcparams);
		pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts, obstacle_squares);
		pbf.simulate(4);
		pbf.printStats("Dam break with obstacles in PBF with sparse marching cubes");
	}
	SECTION("Our Marching Cubes") {
		mcparams.sparse = false;
		mcparams.ours = true;
		PBD::PBF pbf(false, false, "../res/dam_break_pbf_mc_ours/dam_break_pbf_mc_ours_", params, mcparams);
		pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(4);
		pbf.printStats("Dam break in PBF with our marching cubes");
	}
}

TEST_CASE("Dam Break PBF with Surface Construction Time Step Tests 1", "[Dam Break PBF with Surface Construction Time Step Tests 1]")
{
	std::cout << "Testing dam break in PBF with surface construction with maxdt 0.0005" << std::endl;

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
	params.dt = 0.00025;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.005;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.0005;
	params.fluid_pressure_stiffness = 1000.0;
	params.fluid_viscosity = 0.0025;
	params.boundary_viscosity = 0.0;

	params.particle_diameter = 2 * params.particle_radius;
	params.fluid_sampling_distance = params.particle_diameter;
	params.boundary_sampling_distance = 0.8 * params.particle_diameter;
	params.smoothing_length = 1.2 * params.particle_diameter;
	params.smoothing_length_squared = params.smoothing_length * params.smoothing_length;
	params.compact_support = 2 * params.smoothing_length;

	params.grid_cell_size = 1.0 * params.particle_radius;

	params.export_type = PBFParameters::export_type::EXPORT_WITH_SURFACE;

	SECTION("Regular Marching Cubes") {
		mcparams.ours = false;
		mcparams.sparse = false;
		PBD::PBF pbf(false, false, "../res/dam_break_pbf_mc/dam_break_pbf_mc_", params, mcparams);
		pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(4);
		pbf.printStats("Dam break in PBF with regular marching cubes");
	}
	SECTION("Sparse Marching Cubes") {
		mcparams.sparse = true;
		PBD::PBF pbf(false, false, "../res/dam_break_pbf_mc_sparse_maxdt0005/dam_break_pbf_mc_sparse_maxdt0005_", params, mcparams);
		pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(2);
		pbf.printStats("Dam break in PBF with sparse marching cubes");
	}
}

TEST_CASE("Crazy Stuff PBF with Surface Construction", "[Crazy Stuff PBF with Surface Construction]")
{
	std::cout << "Testing crazy stuff in PBF with surface construction" << std::endl;
	PBD::Vector boundary_size = { 1.75, 1, 1.25 };
	PBD::Vector boundary_left = { 0.0, 0.0, 0.0 };
	PBD::Vector fluid_size = { 0.5, 0.5, 0.5 };
	PBD::Vector fluid_left = { 0.02, 0.02, 0.25 };
	PBD::Vector fluid_size2 = { 0.5, 0.5, 0.5 };
	PBD::Vector fluid_left2 = { 1, 0.02, 0.25 };

	std::vector<PBD::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);
	fluid_sizes.push_back(fluid_size2);

	std::vector<PBD::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);
	fluid_lefts.push_back(fluid_left2);

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

	params.grid_cell_size = 1.2 * params.particle_radius;

	params.export_type = PBFParameters::export_type::EXPORT_WITH_SURFACE;

	SECTION("Regular Marching Cubes") {
		mcparams.ours = false;
		mcparams.sparse = false;
		PBD::PBF pbf(false, false, "../res/crazy_stuff_pbf_mc/crazy_stuff_pbf_mc_", params, mcparams);
		pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(1.5);
		pbf.printStats("Crazy stuff in PBF with regular marching cubes");
	}
	SECTION("Sparse Marching Cubes") {
		mcparams.sparse = true;
		PBD::PBF pbf(false, false, "../res/crazy_stuff_pbf_mc_sparse/crazy_stuff_pbf_mc_sparse_", params, mcparams);
		pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(1.5);
		pbf.printStats("Crazy stuff in PBF with sparse marching cubes");
	}

	//SECTION("Surface Nets"){ //TODO for me maybe
	//}
}*/

/*TEST_CASE("Crazy Stuff PBF with Surface Construction2", "[Crazy Stuff PBF with Surface Construction2]")
{
	std::cout << "Testing crazy stuff in PBF with surface construction" << std::endl;
	PBD::Vector boundary_size = { 1.75, 1, 1.25 };
	PBD::Vector boundary_left = { 0.0, 0.0, 0.0 };
	PBD::Vector fluid_size = { 0.5, 0.5, 0.5 };
	PBD::Vector fluid_left = { 0.02, 0.02, 0.25 };
	PBD::Vector fluid_size2 = { 0.5, 0.5, 0.5 };
	PBD::Vector fluid_left2 = { 1, 0.02, 0.25 };

	std::vector<PBD::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);
	fluid_sizes.push_back(fluid_size2);

	std::vector<PBD::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);
	fluid_lefts.push_back(fluid_left2);

	PBFParameters params;
	MCParameters mcparams;
	params.dt = 0.00025;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.005;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.0005;
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

	params.export_type = PBFParameters::export_type::EXPORT_WITH_SURFACE;

	SECTION("Regular Marching Cubes") {
		mcparams.ours = false;
		mcparams.sparse = false;
		PBD::PBF pbf(false, false, "../res/crazy_stuff_pbf_mc2/crazy_stuff_pbf_mc2_", params, mcparams);
		pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(1.5);
		pbf.printStats("Crazy stuff in PBF with regular marching cubes");
	}
	SECTION("Sparse Marching Cubes") {
		mcparams.sparse = true;
		PBD::PBF pbf(false, false, "../res/crazy_stuff_pbf_mc_sparse2/crazy_stuff_pbf_mc_sparse2_", params, mcparams);
		pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(4);
		pbf.printStats("Crazy stuff in PBF with sparse marching cubes");
	}

	//SECTION("Surface Nets"){ //TODO for me maybe
	//}
}*/