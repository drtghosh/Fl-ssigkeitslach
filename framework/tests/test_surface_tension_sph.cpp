#include "catch.hpp"
#include <iostream>
#include <chrono>
#include "../mcubes/marching_cubes.h"
#include "CompactNSearch/Config.h"
#include "../SPH/config.h"
#include "../geometry/io.h"
#include "../SPH/sph.h"
#include "../SPH/params.h"
#include "../mcubes/mcparams.h"

// Check out https://github.com/catchorg/Catch2 for more information about how to use Catch2
TEST_CASE("Surface Tension without Gravity for SPH with surface construction", "[Surface Tension without Gravity for SPH with surface construction]")
{
	std::cout << "Testing surface tension without gravity for SPH with surface construction: " << std::endl;
	WCSPH::Vector boundary_size = { 0.5, 0.5, 0.5};
	WCSPH::Vector boundary_left = { -0.25, -0.25, -0.25 };
	//WCSPH::Vector fluid_size = { 0.2, 0.2, 0.2 }; // Use this for the first two sections
	//WCSPH::Vector fluid_left = { -0.1, -0.1, -0.24 };
	WCSPH::Vector fluid_size = { 0.1, 0.1, 0.1 }; 
	WCSPH::Vector fluid_left = { -0.05, -0.05, -0.24 };

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

	params.export_type = SPHParameters::export_type::EXPORT_WITH_SURFACE;

	/*SECTION("No Boundary, No Surface Tension") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, false, "../res/no_bounds_no_surface_tension_sph/no_bounds_no_surface_tension_sph_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(4);
		sph.printStats("No boundary no surface tension without gravity in SPH with regular marching cubes");
	}

	SECTION("No Boundary, With Surface Tension") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, true, "../res/no_bounds_surface_tension_sph/no_bounds_surface_tension_sph_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(4);
		sph.printStats("No boundary with surface tension without gravity in SPH with regular marching cubes");
	}*/

	SECTION("Boundary, With Surface Tension, Beta 1.0, Gamma 1.0") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, true, "../res/bounds_surface_tension_sph_beta1_gamma1/bounds_surface_tension_sph_beta1_gamma1_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
		sph.printStats("Surface tension without gravity in bounded fluid SPH with regular marching cubes");
	}

	params.adhesion_coefficient = 1.0;
	params.cohesion_coefficient = 0.0;
	SECTION("Boundary, With Surface Tension, Beta 1.0, Gamma 0.0") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, true, "../res/bounds_surface_tension_sph_beta1_gamma0/bounds_surface_tension_sph_beta1_gamma0_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
		sph.printStats("Surface tension without gravity in bounded fluid SPH with regular marching cubes");
	}

	params.adhesion_coefficient = 1.0;
	params.cohesion_coefficient = 0.1;
	SECTION("Boundary, With Surface Tension, Beta 1.0, Gamma 0.1") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, true, "../res/bounds_surface_tension_sph_beta1_gamma01/bounds_surface_tension_sph_beta1_gamma01_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
		sph.printStats("Surface tension without gravity in bounded fluid SPH with regular marching cubes");
	}

	params.adhesion_coefficient = 1.0;
	params.cohesion_coefficient = 0.5;
	SECTION("Boundary, With Surface Tension, Beta 1.0, Gamma 0.5") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, true, "../res/bounds_surface_tension_sph_beta1_gamma05/bounds_surface_tension_sph_beta1_gamma05_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
		sph.printStats("Surface tension without gravity in bounded fluid SPH with regular marching cubes");
	}

	params.adhesion_coefficient = 2.0;
	params.cohesion_coefficient = 0.0;
	SECTION("Boundary, With Surface Tension, Beta 2.0, Gamma 0.0") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, true, "../res/bounds_surface_tension_sph_beta2_gamma0/bounds_surface_tension_sph_beta2_gamma0_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
		sph.printStats("Surface tension without gravity in bounded fluid SPH with regular marching cubes");
	}
	
	params.adhesion_coefficient = 2.0;
	params.cohesion_coefficient = 0.1;
	SECTION("Boundary, With Surface Tension, Beta 2.0, Gamma 0.1") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, true, "../res/bounds_surface_tension_sph_beta2_gamma01/bounds_surface_tension_sph_beta2_gamma01_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
		sph.printStats("Surface tension without gravity in bounded fluid SPH with regular marching cubes");
	}

	params.adhesion_coefficient = 2.0;
	params.cohesion_coefficient = 0.5;
	SECTION("Boundary, With Surface Tension, Beta 2.0, Gamma 0.5") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, true, "../res/bounds_surface_tension_sph_beta2_gamma05/bounds_surface_tension_sph_beta2_gamma05_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
		sph.printStats("Surface tension without gravity in bounded fluid SPH with regular marching cubes");
	}

	params.adhesion_coefficient = 2.0;
	params.cohesion_coefficient = 1.0;
	SECTION("Boundary, With Surface Tension, Beta 2.0, Gamma 1.0") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, true, "../res/bounds_surface_tension_sph_beta2_gamma1/bounds_surface_tension_sph_beta2_gamma1_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
		sph.printStats("Surface tension without gravity in bounded fluid SPH with regular marching cubes");
	}

	params.adhesion_coefficient = 0.001;
	params.cohesion_coefficient = 0.1;
	SECTION("Boundary, With Surface Tension, Beta 0.001, Gamma 0.1") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, true, "../res/bounds_surface_tension_sph_beta0_gamma01/bounds_surface_tension_sph_beta0_gamma01_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
		sph.printStats("Surface tension without gravity in bounded fluid SPH with regular marching cubes");
	}

	params.adhesion_coefficient = 0.001;
	params.cohesion_coefficient = 0.5;
	SECTION("Boundary, With Surface Tension, Beta 0.001, Gamma 0.5") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, true, "../res/bounds_surface_tension_sph_beta0_gamma05/bounds_surface_tension_sph_beta0_gamma05_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
		sph.printStats("Surface tension without gravity in bounded fluid SPH with regular marching cubes");
	}

	params.adhesion_coefficient = 0.001;
	params.cohesion_coefficient = 1.0;
	SECTION("Boundary, With Surface Tension, Beta 0.001, Gamma 1.0") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, true, "../res/bounds_surface_tension_sph_beta0_gamma1/bounds_surface_tension_sph_beta0_gamma1_", params, mcparams);
		sph.turn_off_gravity();
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
		sph.printStats("Surface tension without gravity in bounded fluid SPH with regular marching cubes");
	}
}

/*TEST_CASE("Surface Tension with Gravity for bounded fluid SPH with surface construction", "[Surface Tension with Gravity for bounded fluid SPH with surface construction]")
{
	std::cout << "Testing surface tension with gravity for bounded fluid SPH with surface construction: " << std::endl;
	WCSPH::Vector boundary_size = { 0.5, 0.5, 1.0};
	WCSPH::Vector boundary_left = { -0.25, -0.25, -0.25 };
	WCSPH::Vector fluid_size = { 0.2, 0.2, 0.2 };
	WCSPH::Vector fluid_left = { -0.1, -0.1, -0.1 };

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

	params.export_type = SPHParameters::export_type::EXPORT_WITH_SURFACE;

	SECTION("Boundary, No Surface Tension and No Adhesion") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, false, "../res/bounds_no_surface_tension_sph/bounds_no_surface_tension_sph_", params, mcparams);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(4);
		sph.printStats("No surface tension with gravity in bounded fluid SPH with regular marching cubes");
	}

	SECTION("Boundary, With Surface Tension and Adhesion") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, true, "../res/bounds_surface_tension_sph/bounds_surface_tension_sph_", params, mcparams);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(4);
		sph.printStats("Surface tension with gravity in bounded fluid SPH with regular marching cubes");
	}
}*/

/*TEST_CASE("Surface Tension with Gravity for bounded second fluid SPH with surface construction", "[Surface Tension with Gravity for bounded second fluid SPH with surface construction]")
{
	std::cout << "Testing surface tension with gravity for bounded second fluid SPH with surface construction: " << std::endl;
	WCSPH::Vector boundary_size = { 0.153, 0.153, 0.5 };
	WCSPH::Vector boundary_left = { -0.0765, -0.0765, -0.0765 };
	WCSPH::Vector fluid_size = { 0.02, 0.02, 0.016 };
	WCSPH::Vector fluid_left = { -0.01, -0.01, 0.044 };

	std::vector<WCSPH::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);

	std::vector<WCSPH::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);

	WCSPH::Vector second_fluid_size = { 0.15, 0.15, 0.04 };
	WCSPH::Vector second_fluid_left = { -0.075, -0.075, -0.075 };
	fluid_sizes.push_back(second_fluid_size);
	fluid_lefts.push_back(second_fluid_left);


	SPHParameters params;
	MCParameters mcparams;
	params.dt = 0.0001;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.0005;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.0001;
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

	SECTION("Boundary, No Surface Tension") {
		mcparams.ours = false;
		mcparams.sparse = true;
		WCSPH::SPH sph(false, false, false, "../res/drop_fluid_bounds_no_surface_tension_sph/drop_fluid_bounds_no_surface_tension_sph_", params, mcparams);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(0.1);
		sph.printStats("No surface tension with gravity in bounded second fluid SPH with regular marching cubes");
	}

	SECTION("Boundary, With Surface Tension") {
		mcparams.ours = false;
		mcparams.sparse = true;
		WCSPH::SPH sph(false, false, true, "../res/drop_fluid_bounds_surface_tension_sph/drop_fluid_bounds_surface_tension_sph_", params, mcparams);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(0.01);
		sph.printStats("Surface tension with gravity in bounded second fluid SPH with regular marching cubes");
	}
}*/


/*TEST_CASE("Dam Break SPH with Surface Construction", "[Dam Break SPH with Surface Construction]")
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
