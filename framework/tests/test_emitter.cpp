#include "catch.hpp"
#include "../geometry/emitter.h"
#include "../geometry/io.h"
#include "../SPH/sph.h"
#include "../SPH/params.h"
#include "../PBF/pbf.h"
#include "../PBF/params.h"


// Check out https://github.com/catchorg/Catch2 for more information about how to use Catch2
/*TEST_CASE("Test Emitter Position + Shield", "[Test Emitter Position + Shield]")
{
	std::string path = "../res/emitters/";
	SECTION("Emitter parallel xy") {
		path += "emitter_xy.vtk";
		Emitter::Emitter emitter(5, { 1.0, 0.0, 0.0 }, 3, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 0.0 }, {0.0, 0.0, 1.0}, { 0.0, 0.0, 1.0 }, 0.5, 0.5);
		auto vertices = emitter.get_shield_vertices();
		auto triangles = emitter.get_shield_triangles();

		geometry::write_tri_mesh_to_vtk(path, vertices, triangles);
	}

	SECTION("Emitter parallel xy not normalized") {
		path += "emitter_xy_not_normalized.vtk";
		Emitter::Emitter emitter(5, { 3.0, 0.0, 0.0 }, 3, { 0.0, 3.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 3.0 }, { 0.0, 0.0, 1.0 }, 0.5, 0.5);
		auto vertices = emitter.get_shield_vertices();
		auto triangles = emitter.get_shield_triangles();

		geometry::write_tri_mesh_to_vtk(path, vertices, triangles);
	}

	SECTION("Emitter parallel xz") {
		path += "emitter_xz.vtk";
		Emitter::Emitter emitter(5, { 1.0, 0.0, 0.0 }, 3, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 }, 0.5, 0.5);
		auto vertices = emitter.get_shield_vertices();
		auto triangles = emitter.get_shield_triangles();

		geometry::write_tri_mesh_to_vtk(path, vertices, triangles);
	}

	SECTION("Emitter parallel 45 degrees") {
		path += "emitter_45.vtk";
		Emitter::Emitter emitter(5, { 1.0, 1.0, 0.0 }, 3, { -1.0, 1.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0}, 0.5, 0.5);
		auto vertices = emitter.get_shield_vertices();
		auto triangles = emitter.get_shield_triangles();

		geometry::write_tri_mesh_to_vtk(path, vertices, triangles);
	}

	SECTION("Emitter parallel 45 degrees z") {
		path += "emitter_45z.vtk";
		Emitter::Emitter emitter(5, { 1.0, 1.0, 1.0 }, 3, { 0.5, -1.0, 0.5 }, { 0.0, 0.0, 0.0 }, { -1.5, 0.0, 1.5 }, { -1.0, 0.0, 1.0 }, 0.5, 0.5);
		auto vertices = emitter.get_shield_vertices();
		auto triangles = emitter.get_shield_triangles();

		geometry::write_tri_mesh_to_vtk(path, vertices, triangles);
	}

	SECTION("Emitter parallel 45 degrees z + transform") {
		path += "emitter_45z_transform.vtk";
		Emitter::Emitter emitter(5, { 1.0, 1.0, 1.0 }, 3, { 0.5, -1.0, 0.5 }, { 0.0, 2.0, 3.0 }, { -1.5, 0.0, 1.5 }, { -1.0, 0.0, 1.0 }, 0.5, 0.5);
		auto vertices = emitter.get_shield_vertices();
		auto triangles = emitter.get_shield_triangles();

		geometry::write_tri_mesh_to_vtk(path, vertices, triangles);
	}
}

TEST_CASE("Emitter Emit", "[Emitter Emit]") {
	std::string path = "../res/emitters/";
	SECTION("Emitter parallel xy") {
		path += "emitter_xy_emit.vtk";
		Emitter::Emitter emitter(5, { 1.0, 0.0, 0.0 }, 3, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 }, 0.5, 0.5);
		std::vector<geometry::Vector> particle_positions;
		std::vector<geometry::Vector> particle_velocities;
		emitter.emit_particles(particle_positions, particle_velocities);

		geometry::write_particles_to_vtk(path, particle_positions);
	}

	SECTION("Emitter parallel xy not normalized") {
		path += "emitter_xy_not_normalized_emit.vtk";
		Emitter::Emitter emitter(5, { 3.0, 0.0, 0.0 }, 3, { 0.0, 3.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 3.0 }, { 0.0, 0.0, 1.0 }, 0.5, 0.5);
		std::vector<geometry::Vector> particle_positions;
		std::vector<geometry::Vector> particle_velocities;
		emitter.emit_particles(particle_positions, particle_velocities);

		geometry::write_particles_to_vtk(path, particle_positions);
	}

	SECTION("Emitter parallel xz") {
		path += "emitter_xz_emit.vtk";
		Emitter::Emitter emitter(5, { 1.0, 0.0, 0.0 }, 3, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 }, 0.5, 0.5);
		std::vector<geometry::Vector> particle_positions;
		std::vector<geometry::Vector> particle_velocities;
		emitter.emit_particles(particle_positions, particle_velocities);

		geometry::write_particles_to_vtk(path, particle_positions);
	}

	SECTION("Emitter parallel 45 degrees") {
		path += "emitter_45_emit.vtk";
		Emitter::Emitter emitter(5, { 1.0, 1.0, 0.0 }, 3, { -1.0, 1.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 }, 0.5, 0.5);
		std::vector<geometry::Vector> particle_positions;
		std::vector<geometry::Vector> particle_velocities;
		emitter.emit_particles(particle_positions, particle_velocities);

		geometry::write_particles_to_vtk(path, particle_positions);
	}

	SECTION("Emitter parallel 45 degrees z") {
		path += "emitter_45z_emit.vtk";
		Emitter::Emitter emitter(5, { 1.0, 1.0, 1.0 }, 3, { 0.5, -1.0, 0.5 }, { 0.0, 0.0, 0.0 }, { -1.5, 0.0, 1.5 }, { -1.0, 0.0, 1.0 }, 0.5, 0.5);
		std::vector<geometry::Vector> particle_positions;
		std::vector<geometry::Vector> particle_velocities;
		emitter.emit_particles(particle_positions, particle_velocities);

		geometry::write_particles_to_vtk(path, particle_positions);
	}

	SECTION("Emitter parallel 45 degrees z + transform") {
		path += "emitter_45z_transform_emit.vtk";
		Emitter::Emitter emitter(5, { 1.0, 1.0, 1.0 }, 3, { 0.5, -1.0, 0.5 }, { 0.0, 2.0, 3.0 }, { -1.5, 0.0, 1.5 }, { -1.0, 0.0, 1.0 }, 0.5, 0.5);
		std::vector<geometry::Vector> particle_positions;
		std::vector<geometry::Vector> particle_velocities;
		emitter.emit_particles(particle_positions, particle_velocities);

		geometry::write_particles_to_vtk(path, particle_positions);
	}
}*/

TEST_CASE("Test Emitter Fountain SPH", "[Test Emitter Fountain SPH]")
{
	std::string path = "../res/emitters/";

	WCSPH::Vector boundary_size = { 3.0, 3.0, 3.0 };
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
	params.fluid_viscosity = 0.0025;
	params.boundary_viscosity = 0.0;

	params.particle_diameter = 2 * params.particle_radius;
	params.fluid_sampling_distance = params.particle_diameter;
	params.boundary_sampling_distance = 0.8 * params.particle_diameter;
	params.smoothing_length = 1.2 * params.particle_diameter;
	params.smoothing_length_squared = params.smoothing_length * params.smoothing_length;
	params.compact_support = 2 * params.smoothing_length;
	params.emit_frequency = 1.2;
	params.max_num_emitted_particles = 100000;

	Emitter::Emitter emitter(0.0625, { 1.0, 0.0, 0.0 }, 0.0625, { 0.0, 1.0, 0.0 }, { 1.5, 1.5, 0.1 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 5 }, params.particle_diameter, params.particle_diameter);

	/*SECTION("Fountain SPH") {
		path += "fountain_sph.vtk";

		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.0);
		schedule.push_back(2.1);
		emitter.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter);

		WCSPH::SPH sph(false, false, "../res/emitter_fountain_sph/emitter_fountain_sph_", params, mcparams, emitters);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
	}*/

	SECTION("Fountain SPH Schedule") {
		path += "fountain_sph_schedule.vtk";

		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.2);
		schedule.push_back(0.5);
		schedule.push_back(0.8);
		schedule.push_back(1.2);
		schedule.push_back(1.8);
		schedule.push_back(2.1);
		emitter.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter);

		WCSPH::SPH sph(false, false, "../res/emitter_fountain_sph_schedule/emitter_fountain_sph_schedule_", params, mcparams, emitters);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
	}
}

TEST_CASE("Test Emitter Fountain PBF", "[Test Emitter Fountain PBF]")
{

}