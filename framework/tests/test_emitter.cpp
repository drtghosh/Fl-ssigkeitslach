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