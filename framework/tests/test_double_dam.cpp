#include "catch.hpp"
#include "../geometry/emitter.h"
#include "../geometry/io.h"
#include "../SPH/sph.h"
#include "../SPH/params.h"
#include "../PBF/pbf.h"
#include "../PBF/params.h"


// Check out https://github.com/catchorg/Catch2 for more information about how to use Catch2
/*TEST_CASE("Double Dam Break SPH", "[Double Dam Break SPH]") {
	WCSPH::Vector boundary_size = { 0.18, 1.4, 1.0 };
	WCSPH::Vector boundary_left = { -0.015, -0.015, -0.015 };
	WCSPH::Vector fluid_size = { 0.15, 0.25, 0.5 };
	WCSPH::Vector fluid_left = { 0.0, 0.0, 0.0 };
	WCSPH::Vector fluid_left2 = { 0.0, 1.12, 0.0 };

	std::vector<WCSPH::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);
	fluid_sizes.push_back(fluid_size);

	std::vector<WCSPH::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);
	fluid_lefts.push_back(fluid_left2);

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
	params.grid_cell_size = 1.2 * params.particle_radius;

	params.export_type = params.EXPORT_WITH_SURFACE;

	WCSPH::SPH sph(false, false, true, "../res/double_dam_break_sph/double_dam_break_sph_", params, mcparams);
	sph.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	sph.simulate(2);
	sph.printStats("Double Dam Break SPH");

	/*Data:
	Timesteps: 8000
	Boundary Particles: 74630
	Fluid Particles: 39738
	Runtime: 1437310 ms
	*/
/*}

TEST_CASE("Double Dam Break PBF", "[Double Dam Break PBF]") {
	PBD::Vector boundary_size = { 0.18, 1.4, 1.0 };
	PBD::Vector boundary_left = { -0.015, -0.015, -0.015 };
	PBD::Vector fluid_size = { 0.15, 0.25, 0.5 };
	PBD::Vector fluid_left = { 0.0, 0.0, 0.0 };
	PBD::Vector fluid_left2 = { 0.0, 1.12, 0.0 };

	std::vector<PBD::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);
	fluid_sizes.push_back(fluid_size);

	std::vector<PBD::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);
	fluid_lefts.push_back(fluid_left2);

	PBFParameters params;
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
	params.grid_cell_size = 1.2 * params.particle_radius;

	params.export_type = params.EXPORT_WITH_SURFACE;
	params.pbf_iterations = 5;

	PBD::PBF pbf(false, false, true, "../res/double_dam_break_pbf/double_dam_break_pbf_", params, mcparams);
	pbf.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	pbf.simulate(2);
	pbf.printStats("Double Dam Break PBF");

	/*Data:
	Timesteps: 8000
	Boundary Particles: 74630
	Fluid Particles: 39780
	Runtime: 3444449 ms
	*/
/*}

TEST_CASE("Double Dam Break SPH Particles", "[Double Dam Break SPH Particles]") {
	WCSPH::Vector boundary_size = { 0.18, 1.4, 1.0 };
	WCSPH::Vector boundary_left = { -0.015, -0.015, -0.015 };
	WCSPH::Vector fluid_size = { 0.15, 0.25, 0.5 };
	WCSPH::Vector fluid_left = { 0.0, 0.0, 0.0 };
	WCSPH::Vector fluid_left2 = { 0.0, 1.12, 0.0 };

	std::vector<WCSPH::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);
	fluid_sizes.push_back(fluid_size);

	std::vector<WCSPH::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);
	fluid_lefts.push_back(fluid_left2);

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
	params.grid_cell_size = 1.2 * params.particle_radius;

	WCSPH::SPH sph(false, false, true, "../res/double_dam_break_sph_particles/double_dam_break_sph_particles_", params, mcparams);
	sph.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	sph.simulate(2);
	sph.printStats("Double Dam Break SPH");

	/*Data:
	Timesteps: 8000
	Boundary Particles: 74630
	Fluid Particles: 39731
	Runtime: 1128054/1841819 ms
	*/
/*}

TEST_CASE("Double Dam Break PBF Particles", "[Double Dam Break Particles]") {
	PBD::Vector boundary_size = { 0.18, 1.4, 1.0 };
	PBD::Vector boundary_left = { -0.015, -0.015, -0.015 };
	PBD::Vector fluid_size = { 0.15, 0.25, 0.5 };
	PBD::Vector fluid_left = { 0.0, 0.0, 0.0 };
	PBD::Vector fluid_left2 = { 0.0, 1.12, 0.0 };

	std::vector<PBD::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);
	fluid_sizes.push_back(fluid_size);

	std::vector<PBD::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);
	fluid_lefts.push_back(fluid_left2);

	PBFParameters params;
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
	params.grid_cell_size = 1.2 * params.particle_radius;

	params.export_type = params.EXPORT_WITH_SURFACE;
	params.pbf_iterations = 5;

	PBD::PBF pbf(false, false, true, "../res/double_dam_break_pbf_particles/double_dam_break_pbf_particles_", params, mcparams);
	pbf.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	pbf.simulate(2);
	pbf.printStats("Double Dam Break PBF");

	/*Data:
	Timesteps: 8000
	Boundary Particles: 74630
	Fluid Particles: 39780
	Runtime: 3068366 ms
	*/
/*}

TEST_CASE("Double Dam Break SPH Max Dt", "[Double Dam Break SPH Max Dt]") {
	WCSPH::Vector boundary_size = { 0.18, 1.4, 1.0 };
	WCSPH::Vector boundary_left = { -0.015, -0.015, -0.015 };
	WCSPH::Vector fluid_size = { 0.15, 0.25, 0.5 };
	WCSPH::Vector fluid_left = { 0.0, 0.0, 0.0 };
	WCSPH::Vector fluid_left2 = { 0.0, 1.12, 0.0 };

	std::vector<WCSPH::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);
	fluid_sizes.push_back(fluid_size);

	std::vector<WCSPH::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);
	fluid_lefts.push_back(fluid_left2);

	SPHParameters params;
	MCParameters mcparams;
	params.dt = 0.0001;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.005;
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
	params.grid_cell_size = 1.2 * params.particle_radius;

	params.export_type = params.EXPORT_WITH_SURFACE;

	/*SECTION("Max DT: 0.1 ms") {
		WCSPH::SPH sph(false, false, true, "../res/double_dam_break_sph_dt0_1/double_dam_break_sph_dt0_1_", params, mcparams);
		sph.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
		sph.printStats("Double Dam Break SPH Max Dt");

		/*Data:
		Timesteps: 20001
		Boundary Particles: 74630
		Fluid Particles: 39734
		Runtime: 2783866 ms
		*/
	/*}

	SECTION("Max DT: 0.5 ms") {
		params.max_dt = 0.0005;
		WCSPH::SPH sph(false, false, true, "../res/double_dam_break_sph_dt0_5/double_dam_break_sph_dt0_5_", params, mcparams);
		sph.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
		sph.printStats("Double Dam Break SPH Max Dt");

		/*Data:
		Timesteps: 4001
		Boundary Particles: 74630
		Fluid Particles: 39731
		Runtime: 627802 ms
		*/
	/*}

	SECTION("Max DT: 1 ms") {
		params.max_dt = 0.001
			;
		WCSPH::SPH sph(false, false, true, "../res/double_dam_break_sph_dt1/double_dam_break_sph_dt1_", params, mcparams);
		sph.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
		sph.printStats("Double Dam Break SPH Max Dt");

		/*Data:
		Timesteps: 2001
		Boundary Particles: 74630
		Fluid Particles: 7097
		Runtime: 235096 ms
		*/
	/*}
	
}

TEST_CASE("Double Dam Break PBF Max Dt", "[Double Dam Break PBF Max Dt]") {
	PBD::Vector boundary_size = { 0.18, 1.4, 1.0 };
	PBD::Vector boundary_left = { -0.015, -0.015, -0.015 };
	PBD::Vector fluid_size = { 0.15, 0.25, 0.5 };
	PBD::Vector fluid_left = { 0.0, 0.0, 0.0 };
	PBD::Vector fluid_left2 = { 0.0, 1.12, 0.0 };

	std::vector<PBD::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);
	fluid_sizes.push_back(fluid_size);

	std::vector<PBD::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);
	fluid_lefts.push_back(fluid_left2);

	PBFParameters params;
	MCParameters mcparams;
	params.dt = 0.0001;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.005;
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
	params.grid_cell_size = 1.2 * params.particle_radius;

	params.export_type = params.EXPORT_WITH_SURFACE;
	params.pbf_iterations = 5;

	/*SECTION("Max DT: 0.1 ms") {
		PBD::PBF pbf(false, false, true, "../res/double_dam_break_pbf_dt0_1/double_dam_break_pbf_dt0_1_", params, mcparams);
		pbf.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(2);
		pbf.printStats("Double Dam Break PBF Max Dt");

		/*Data:
		Timesteps: 20001
		Boundary Particles: 74630
		Fluid Particles: 39780
		Runtime: 4440000 (74 min) ms
		*/
	/*}

	SECTION("Max DT: 0.5 ms") {
		params.max_dt = 0.0005;
		PBD::PBF pbf(false, false, true, "../res/double_dam_break_pbf_dt0_5/double_dam_break_pbf_dt0_5_", params, mcparams);
		pbf.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(2);
		pbf.printStats("Double Dam Break PBF Max Dt");

		/*Data:
		Timesteps: 4001
		Boundary Particles: 74630
		Fluid Particles: 39780
		Runtime: 1913284 ms
		*/
	/*}

	SECTION("Max DT: 1 ms") {
		params.max_dt = 0.001;
		PBD::PBF pbf(false, false, true, "../res/double_dam_break_pbf_dt1/double_dam_break_pbf_dt1_", params, mcparams);
		pbf.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(2);
		pbf.printStats("Double Dam Break PBF Max Dt");

		/*Data:
		Timesteps: 2001
		Boundary Particles: 74630
		Fluid Particles: 39775
		Runtime: 844643 ms
		*/
	/*}

}

/*TEST_CASE("Double Dam Break PBF Iter", "[Double Dam Break PBF Iter]") {
	PBD::Vector boundary_size = { 0.18, 1.4, 1.0 };
	PBD::Vector boundary_left = { -0.015, -0.015, -0.015 };
	PBD::Vector fluid_size = { 0.15, 0.25, 0.5 };
	PBD::Vector fluid_left = { 0.0, 0.0, 0.0 };
	PBD::Vector fluid_left2 = { 0.0, 1.12, 0.0 };

	std::vector<PBD::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);
	fluid_sizes.push_back(fluid_size);

	std::vector<PBD::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);
	fluid_lefts.push_back(fluid_left2);

	PBFParameters params;
	MCParameters mcparams;
	params.dt = 0.001;
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

	params.export_type = params.EXPORT_WITH_SURFACE;
	params.pbf_iterations = 2;

	SECTION("Iter: 2") {
		PBD::PBF pbf(false, false, true, "../res/double_dam_break_pbf_iter2/double_dam_break_pbf_iter2_", params, mcparams);
		pbf.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(2);
		pbf.printStats("Double Dam Break PBF Iter");

		/*Data:
		Timesteps: 2001
		Boundary Particles: 74630
		Fluid Particles: 28364/39780
		Runtime: 580648 ms
		*/
	/*}

	SECTION("Iter: 5") {
		params.pbf_iterations = 5;
		PBD::PBF pbf(false, false, true, "../res/double_dam_break_pbf_iter5/double_dam_break_pbf_iter5_", params, mcparams);
		pbf.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(2);
		pbf.printStats("Double Dam Break PBF Iter");

		/*Data:
		Timesteps: 2001
		Boundary Particles: 74630
		Fluid Particles: 39769
		Runtime: 879609 ms
		*/
	/*}

	SECTION("Iter: 10") {
		params.pbf_iterations = 10;
		PBD::PBF pbf(false, false, true, "../res/double_dam_break_pbf_iter10/double_dam_break_pbf_iter10_", params, mcparams);
		pbf.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(2);
		pbf.printStats("Double Dam Break PBF Iter");

		/*Data:
		Timesteps: 2001
		Boundary Particles: 74630
		Fluid Particles: 39772
		Runtime: 1466351 ms
		*/
	/*}

}*/

TEST_CASE("Double Dam Break PBF Max Dt Iter", "[Double Dam Break PBF Max Dt Iter]") {
	PBD::Vector boundary_size = { 0.18, 1.4, 1.0 };
	PBD::Vector boundary_left = { -0.015, -0.015, -0.015 };
	PBD::Vector fluid_size = { 0.15, 0.25, 0.5 };
	PBD::Vector fluid_left = { 0.0, 0.0, 0.0 };
	PBD::Vector fluid_left2 = { 0.0, 1.12, 0.0 };

	std::vector<PBD::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);
	fluid_sizes.push_back(fluid_size);

	std::vector<PBD::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);
	fluid_lefts.push_back(fluid_left2);

	PBFParameters params;
	MCParameters mcparams;
	params.dt = 0.0001;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.005;
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
	params.grid_cell_size = 1.2 * params.particle_radius;

	params.export_type = params.EXPORT_WITH_SURFACE;
	params.pbf_iterations = 5;

	SECTION("Max DT: 0.5 ms, Iter: 3") {
		params.max_dt = 0.0005;
		params.pbf_iterations = 3;
		PBD::PBF pbf(false, false, true, "../res/double_dam_break_pbf_dt0_5_iter3/double_dam_break_pbf_dt0_5_iter3_", params, mcparams);
		pbf.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(2);
		pbf.printStats("Double Dam Break PBF Max Dt");

		/*Data: (Iter 2)
		Timesteps: 4001
		Boundary Particles: 74630
		Fluid Particles: 23052
		Runtime: 779979 ms
		*/
		/*Data: (Iter 3)
		Timesteps: 4001
		Boundary Particles: 74630
		Fluid Particles: 39780
		Runtime: 1184690 ms
		*/
		}

	/*SECTION("Max DT: 1 ms, Iter: 5") {
		params.max_dt = 0.001;
		params.pbf_iterations = 5;
		PBD::PBF pbf(false, false, true, "../res/double_dam_break_pbf_dt1_iter5/double_dam_break_pbf_dt1_iter5_", params, mcparams);
		pbf.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(2);
		pbf.printStats("Double Dam Break PBF Max Dt");

		/*Data:
		Timesteps: 2001
		Boundary Particles: 74630
		Fluid Particles: 39774
		Runtime: 854354 ms
		*/
	/* }

	SECTION("Max DT: 5 ms, Iter: 10") {
		params.max_dt = 0.005;
		params.pbf_iterations = 10;
		PBD::PBF pbf(false, false, true, "../res/double_dam_break_pbf_dt5_iter10/double_dam_break_pbf_dt5_iter10_", params, mcparams);
		pbf.load_geometry(true, false, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		pbf.simulate(2);
		pbf.printStats("Double Dam Break PBF Max Dt");

		/*Data:
		Timesteps: 783
		Boundary Particles: 74630
		Fluid Particles: 39696
		Runtime: 589314 ms
		*/
	//}

}