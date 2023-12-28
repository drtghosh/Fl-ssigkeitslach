#include "catch.hpp"
#include "../geometry/emitter.h"
#include "../geometry/io.h"


// Check out https://github.com/catchorg/Catch2 for more information about how to use Catch2
TEST_CASE("Test Emitter", "[Test Emitter]")
{
	std::string path = "../res/emitters/";
	SECTION("Emitter parallel xy") {
		path += "emitter_xy.vtk";
		Emitter::Emitter emitter(5, { 1.0, 0.0, 0.0 }, 3, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 0.0 }, {0.0, 0.0, 1.0}, 0.0, 0.5);
		auto vertices = emitter.get_shield_vertices();
		auto triangles = emitter.get_shield_triangles();

		geometry::write_tri_mesh_to_vtk(path, vertices, triangles);
	}

	SECTION("Emitter parallel xy not normalized") {
		path += "emitter_xy_not_normalized.vtk";
		Emitter::Emitter emitter(5, { 3.0, 0.0, 0.0 }, 3, { 0.0, 3.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 3.0 }, 0.0, 0.5);
		auto vertices = emitter.get_shield_vertices();
		auto triangles = emitter.get_shield_triangles();

		geometry::write_tri_mesh_to_vtk(path, vertices, triangles);
	}

	SECTION("Emitter parallel xz") {
		path += "emitter_xz.vtk";
		Emitter::Emitter emitter(5, { 1.0, 0.0, 0.0 }, 3, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, 0.0, 0.5);
		auto vertices = emitter.get_shield_vertices();
		auto triangles = emitter.get_shield_triangles();

		geometry::write_tri_mesh_to_vtk(path, vertices, triangles);
	}

	SECTION("Emitter parallel 45 degrees") {
		path += "emitter_45.vtk";
		Emitter::Emitter emitter(5, { 1.0, 1.0, 0.0 }, 3, { -1.0, 1.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }, 0.0, 0.5);
		auto vertices = emitter.get_shield_vertices();
		auto triangles = emitter.get_shield_triangles();

		geometry::write_tri_mesh_to_vtk(path, vertices, triangles);
	}
	
}

TEST_CASE("Test Emitter Fountain", "[Test Emitter Fountain]")
{

}