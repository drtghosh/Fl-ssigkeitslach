#include "catch.hpp"
#include <iostream>
#include <chrono>
#include "../mcubes/marching_cubes.h"
#include "CompactNSearch/Config.h"
#include "../SPH/config.h"
#include "../geometry/io.h"
#include "../SPH/sph_out_bound.h"
#include "../SPH/params.h"
#include "../mcubes/mcparams.h"

// Check out https://github.com/catchorg/Catch2 for more information about how to use Catch2
/*TEST_CASE("Surface Tension with Gravity for unbounded fluid SPH with surface construction", "[Surface Tension with Gravity for unbounded fluid SPH with surface construction]")
{
	std::cout << "Testing surface tension with gravity for unbounded fluid SPH with surface construction: " << std::endl;
	WCSPH::Vector boundary_size = { 0.5, 0.5, 0.1 };
	WCSPH::Vector boundary_left = { -0.25, -0.25, -0.5 };
	WCSPH::Vector fluid_size = { 0.2, 0.2, 0.2 };
	WCSPH::Vector fluid_left = { -0.1, -0.1, 0.1 };

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

	params.grid_cell_size = 1.0 * params.particle_radius;

	params.export_type = SPHParameters::export_type::EXPORT_WITH_SURFACE;

	SECTION("Outer Boundary, No Surface Tension") {
		mcparams.ours = false;
		mcparams.sparse = true;
		WCSPH::SPHOB sph(false, false, false, "../res/out_bounds_no_surface_tension_sph/out_bounds_no_surface_tension_sph_", params, mcparams);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(1);
		sph.printStats("No surface tension with gravity in unbounded fluid SPH with regular marching cubes");
	}

	SECTION("Outer Boundary, With Surface Tension") {
		mcparams.ours = false;
		mcparams.sparse = true;
		WCSPH::SPHOB sph(false, false, true, "../res/out_bounds_surface_tension_sph/out_bounds_surface_tension_sph_", params, mcparams);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(1);
		sph.printStats("Surface tension with gravity in unbounded fluid SPH with regular marching cubes");
	}
}


TEST_CASE("Dam Break SPH with Surface Construction", "[Dam Break SPH with Surface Construction]")
{
	std::cout << "Testing dam break in SPH with surface construction" << std::endl;

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

	params.grid_cell_size = 1.0 * params.particle_radius;

	params.export_type = SPHParameters::export_type::EXPORT_WITH_SURFACE;

	SECTION("Regular Marching Cubes") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, false, "../res/dam_break_sph_mc/dam_break_sph_mc_", params, mcparams);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(4);
		sph.printStats("Dam break in SPH with regular marching cubes");
	}
	SECTION("Sparse Marching Cubes") {
		mcparams.sparse = true;
		WCSPH::SPH sph(false, false, false, "../res/dam_break_sph_mc_sparse/dam_break_sph_mc_sparse_", params, mcparams);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(0.5);
		sph.printStats("Dam break in SPH with sparse marching cubes");
	}
	SECTION("Our Marching Cubes") {
		mcparams.sparse = false;
		mcparams.ours = true;
		WCSPH::SPH sph(false, false, false, "../res/dam_break_sph_mc_ours/dam_break_sph_mc_ours_", params, mcparams);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(0.5);
		sph.printStats("Dam break in SPH with our marching cubes");
	}
	//SECTION("Surface Nets"){ //TODO for me maybe
	//}
}*/
