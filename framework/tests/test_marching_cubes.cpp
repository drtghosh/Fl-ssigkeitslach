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
/*TEST_CASE("Marching Cubes", "[Marching Cubes]")
{
	std::cout << "Testing marching cubes" << std::endl;
	CompactNSearch::Real r_minor = 0.2;
	CompactNSearch::Real r_major = 0.7;

	auto torus = [&](WCSPH::Vector p) { //grid dim 2x2x2 [-1,+1]
		CompactNSearch::Real root = std::sqrt(p.x() * p.x() + p.y() * p.y());
		root -= r_major;

		return r_minor * r_minor - root * root - p.z() * p.z();
	};

	auto torus_gradient = [&](WCSPH::Vector p) {
		CompactNSearch::Real sqrt_xy = sqrt(p[0] * p[0] + p[1] * p[1]);
		WCSPH::Vector grad = { (sqrt_xy - r_major) * p[0] / sqrt_xy, (sqrt_xy - r_major) * p[1] / sqrt_xy, p[2] };
		grad.normalize();
		return grad;
	};
	CompactNSearch::Real cell_size = 0;
	std::string file = "../res/torus/";

	SECTION("Cell size: 0.1") {
		cell_size = 0.1;
		file += "torus_0.1.vtk";
	}
	SECTION("Cell size: 0.04") {
		cell_size = 0.04;
		file += "torus_0.04.vtk";
	}
	SECTION("Cell size: 0.01") {
		cell_size = 0.01;
		file += "torus_0.01.vtk";
	}
	SECTION("Cell size: 0.005") {
		cell_size = 0.005;
		file += "torus_0.005.vtk";
	}

	std::vector<WCSPH::Vector> positions;
	std::vector<CompactNSearch::Real> values;
	for (CompactNSearch::Real x = -1; x <= 1.0001; x += cell_size) {
		for (CompactNSearch::Real y = -1; y <= 1.0001; y += cell_size) {
			for (CompactNSearch::Real z = -1; z <= 1.0001; z += cell_size) {
				WCSPH::Vector p = {x,y,z};
				CompactNSearch::Real res = torus(p);
				values.emplace_back(res);
				positions.emplace_back(p);
			}
		}
	}
	std::array<int, 3 > dimensions = { 2 / cell_size + 1 , 2 / cell_size + 1 , 2 / cell_size + 1 };

	std::cout << "Grid size: " << positions.size() << std::endl;

	std::vector<WCSPH::Vector> vertices;
	std::vector<std::array<int, 3>> triangles; 
	std::vector<WCSPH::Vector> normals;


	auto start = std::chrono::system_clock::now();

	Parameters params;
	MCParameters mcparams;
	mcparams.sparse = false;

	std::unordered_map<uint64_t, CompactNSearch::Real> value_map;
	MC::MarchingCubes marching_cubes(mcparams, 0.0);
	marching_cubes.calculate(positions, values, value_map, dimensions, vertices, triangles);
	
	auto end = std::chrono::system_clock::now();
	auto elapsed =
		std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "Vertices size: " << vertices.size() << std::endl;
	std::cout << "Marching cubes time: " << elapsed.count() << " ms" << std::endl;

	for (unsigned int i = 0; i < vertices.size();i++) {
		WCSPH::Vector p = vertices[i];
		WCSPH::Vector normal = torus_gradient(p);
		normals.emplace_back(normal);
	}

	geometry::write_tri_mesh_to_vtk(file, vertices, triangles, normals);
}

TEST_CASE("Dam Break with Surface Construction", "[Dam Break with Surface Construction]")
{
	std::cout << "Testing dam break with surface construction" << std::endl;

	WCSPH::Vector boundary_size = { 0.18, 0.8, 1.0 };
	WCSPH::Vector boundary_left = { -0.015, -0.015, -0.015 };
	WCSPH::Vector fluid_size = { 0.15, 0.25, 0.5 };
	WCSPH::Vector fluid_left = { 0.0, 0.0, 0.0 };

	std::vector<WCSPH::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);

	std::vector<WCSPH::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);

	Parameters params;
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

	params.export_type = Parameters::export_type::EXPORT_WITH_SURFACE;

	SECTION("Regular Marching Cubes") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, "../res/dam_break_mc/dam_break_mc_", params, mcparams);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(0.5);
		sph.printStats("Dam break with regular marching cubes");
	}
	SECTION("Sparse Marching Cubes") {
		mcparams.sparse = true;
		WCSPH::SPH sph(false, false, "../res/dam_break_mc_sparse/dam_break_mc_sparse_", params, mcparams);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(0.5);
		sph.printStats("Dam break with sparse marching cubes");
	}
	SECTION("Our Marching Cubes") {
		mcparams.sparse = false;
		mcparams.ours = true;
		WCSPH::SPH sph(false, false, "../res/dam_break_mc_ours/dam_break_mc_ours_", params, mcparams);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(0.5);
		sph.printStats("Dam break with our marching cubes");
	}
	//SECTION("Surface Nets"){ //TODO for me maybe
	//}
}

TEST_CASE("Crazy Stuff with Surface Construction", "[Crazy Stuff with Surface Construction]")
{
	std::cout << "Testing crazy stuff with surface construction" << std::endl;
	WCSPH::Vector boundary_size = { 1.75, 1, 1.25 };
	WCSPH::Vector boundary_left = { 0.0, 0.0, 0.0 };
	WCSPH::Vector fluid_size = { 0.5, 0.5, 0.5 };
	WCSPH::Vector fluid_left = { 0.02, 0.02, 0.25 };
	WCSPH::Vector fluid_size2 = { 0.5, 0.5, 0.5 };
	WCSPH::Vector fluid_left2 = { 1, 0.02, 0.25 };

	std::vector<WCSPH::Vector> fluid_sizes;
	fluid_sizes.push_back(fluid_size);
	fluid_sizes.push_back(fluid_size2);

	std::vector<WCSPH::Vector> fluid_lefts;
	fluid_lefts.push_back(fluid_left);
	fluid_lefts.push_back(fluid_left2);

	Parameters params;
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

	params.grid_cell_size = 1.2 * params.particle_radius;

	params.export_type = Parameters::export_type::EXPORT_WITH_SURFACE;

	SECTION("Regular Marching Cubes") {
		mcparams.ours = false;
		mcparams.sparse = false;
		WCSPH::SPH sph(false, false, "../res/crazy_stuff_mc/crazy_stuff_mc_", params, mcparams);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(0.1);
		sph.printStats("Crazy stuff with marching cubes");
	}
	SECTION("Sparse Marching Cubes") {
		mcparams.sparse = true;
		WCSPH::SPH sph(false, false, "../res/crazy_stuff_mc_sparse/crazy_stuff_mc_sparse_", params, mcparams);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(0.1);
		sph.printStats("Crazy stuff with sparse marching cubes");
	}

	//SECTION("Surface Nets"){ //TODO for me maybe
	//}
}*/
