#include "catch.hpp"
#include "../geometry/emitter.h"
#include "../geometry/io.h"
#include "../SPH/sph.h"
#include "../SPH/params.h"
#include "../PBF/pbf.h"
#include "../PBF/params.h"
#include <cmath>

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

/*TEST_CASE("Test Emitter Fountain SPH Open Box", "[Test Emitter Fountain SPH Open Box]")
{
	//WCSPH::Vector boundary_size = { 3.0, 3.0, 0.5 }; // Use this in general
	WCSPH::Vector boundary_size = { 3.0, 3.0, 1.0 }; // Use this for with obstacle
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
	params.fluid_viscosity = 0.005;
	params.boundary_viscosity = 0.0;

	params.particle_diameter = 2 * params.particle_radius;
	params.fluid_sampling_distance = params.particle_diameter;
	params.boundary_sampling_distance = 0.8 * params.particle_diameter;
	params.smoothing_length = 1.2 * params.particle_diameter;
	params.smoothing_length_squared = params.smoothing_length * params.smoothing_length;
	params.compact_support = 2 * params.smoothing_length;
	params.grid_cell_size = 1.2 * params.particle_radius;
	params.emit_frequency = 1.2;
	params.max_num_particles = 100000;

	Emitter::Emitter emitter(0.0625, { 1.0, 0.0, 0.0 }, 0.0625, { 0.0, 1.0, 0.0 }, { 1.5, 1.5, 0.1 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 5 }, params.particle_diameter, params.particle_diameter);

	SECTION("Fountain SPH") {
		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.0);
		schedule.push_back(2.1);
		emitter.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter);

		WCSPH::SPH sph(false, false, true, "../res/emitter_fountain_sph/emitter_fountain_sph_", params, mcparams, emitters);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
	}

	SECTION("Fountain SPH Schedule") {
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

	/*SECTION("Fountain SPH Surface Max Particle") {
		params.max_num_particles = 5000;

		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.0);
		schedule.push_back(2.1);
		emitter.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter);

		WCSPH::SPH sph(false, false, "../res/emitter_fountain_sph_max/emitter_fountain_sph_max_", params, mcparams, emitters);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(0.5);
	}

	SECTION("Fountain SPH Surface without surface tension") {
		params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
		mcparams.ours = false;
		mcparams.sparse = false;
		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.0);
		schedule.push_back(2.1);
		emitter.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter);

		WCSPH::SPH sph(false, false, false, "../res/emitter_fountain_sph_surface_ob/emitter_fountain_sph_surface_ob_", params, mcparams, emitters);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
	}

	SECTION("Fountain SPH Surface with surface tension") {
		params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
		mcparams.ours = false;
		mcparams.sparse = false;
		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.0);
		schedule.push_back(2.1);
		emitter.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter);

		WCSPH::SPH sph(false, false, true, "../res/emitter_fountain_sph_surface_ob_st/emitter_fountain_sph_surface_ob_st_", params, mcparams, emitters);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
		sph.simulate(2);
	}

	Emitter::Emitter emitter1(0.0625, { -0.5, -1.0, -1.0 }, 0.0625, { -1.0, -0.5, 1.0 }, { 1.65, 1.35, 2.0 }, { -1.5, 1.5, -0.75 }, { -1.5, 1.5, -0.75 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter2(0.0625, { 0.5, 1.0, -1.0 }, 0.0625, { 1.0, 0.5, 1.0 }, { 1.35, 1.65, 2.0 }, { 1.5, -1.5, -0.75 }, { 1.5, -1.5, -0.75 }, params.particle_diameter, params.particle_diameter);
	WCSPH::Vector obstacle_bl = { 1.0, 1.0, 0.6 };
	WCSPH::Vector obstacle_size = { 1.0, 1.0, 0.25 };
	std::pair <WCSPH::Vector, WCSPH::Vector>& obstacle_box1 = std::pair <WCSPH::Vector, WCSPH::Vector>();
	obstacle_box1.first = obstacle_bl;
	obstacle_box1.second = obstacle_size;
	SECTION("Two hoses and an obstacle with surface tension in SPH") {
		params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
		mcparams.ours = false;
		mcparams.sparse = false;
		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.0);
		schedule.push_back(2.1);
		emitter1.set_schedule(schedule);
		emitter2.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter1);
		emitters.push_back(emitter2);

		WCSPH::SPH sph(false, false, true, "../res/emitter_2hose_1obstacle_sph_surface_ob_st/emitter_2hose_1obstacle_sph_surface_ob_st_", params, mcparams, emitters);
		sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts, obstacle_box1);
		sph.simulate(1.8);
	}
	
	WCSPH::Vector boundary_size_v = { 3.0, 3.0, 0.75 };
	WCSPH::Vector boundary_left_v = { 0.0, 0.0, 0.0 };
	Emitter::Emitter emitter3(0.075, { 1.0, 0.0, 0.0 }, 0.075, { 0.0, 1.0, 0.0 }, { 1.5, 1.5, 1.85 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 3.0 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter4(0.075, { 0.0, 1.0, 0.0 }, 0.075, { 1.0, 0.0, 0.0 }, { 1.5, 1.5, 2.15 }, { 0.0, 0.0, -1.0 }, { 0.0, 0.0, -3.0 }, params.particle_diameter, params.particle_diameter);
	SECTION("Two vertical hoses with surface tension in SPH") {
		params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
		mcparams.ours = false;
		mcparams.sparse = false;
		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.0);
		schedule.push_back(2.1);
		emitter3.set_schedule(schedule);
		emitter4.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter3);
		emitters.push_back(emitter4);

		WCSPH::SPH sph(false, false, true, "../res/emitter_2vertical_hoses_sph_surface_ob_st/emitter_2vertical_hoses_sph_surface_ob_st_", params, mcparams, emitters);
		sph.load_geometry(true, boundary_size_v, boundary_left_v, fluid_sizes, fluid_lefts);
		sph.simulate(1.2);
	}

	Emitter::Emitter emitter5(0.075, { 0.0, 1.0, 0.0 }, 0.075, { 0.0, 0.0, 1.0 }, { 1.35, 1.5, 2.0 }, { 1.0, 0.0, 0.0 }, { 3.0, 0.0, 0.0 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter6(0.075, { 0.0, 0.0, 1.0 }, 0.075, { 0.0, 1.0, 0.0 }, { 1.65, 1.5, 2.0 }, { -1.0, 0.0, 0.0 }, { -3.0, 0.0, 0.0 }, params.particle_diameter, params.particle_diameter);
	SECTION("Four orthogonal hoses with surface tension in SPH") {
		params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
		mcparams.ours = false;
		mcparams.sparse = false;
		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.0);
		schedule.push_back(1.5);
		emitter3.set_schedule(schedule);
		emitter4.set_schedule(schedule);
		emitter5.set_schedule(schedule);
		emitter6.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter3);
		emitters.push_back(emitter4);
		emitters.push_back(emitter5);
		emitters.push_back(emitter6);

		WCSPH::SPH sph(false, false, true, "../res/emitter_4orthogonal_hoses_sph_surface_ob_st/emitter_4orthogonal_hoses_sph_surface_ob_st_", params, mcparams, emitters);
		sph.load_geometry(true, boundary_size_v, boundary_left_v, fluid_sizes, fluid_lefts);
		sph.simulate(1.2);
	}

	//WCSPH::Vector boundary_size_s = { 4.0, 4.0, 0.75 };
	//WCSPH::Vector boundary_left_s = { -2.0, -2.0, -1.0 };
	//Emitter::Emitter emitter_sp(0.0625, { 0.0, 1.0, 0.0 }, 0.0625, { 1.0, 0.0, 0.0 }, { 0.0, 0.0, 0.2 }, { 0.0, 0.0, -1.0 }, { 0.0, 0.0, -0.5 }, params.particle_diameter, params.particle_diameter);
	SECTION("Sphere obstacle with surface tension in SPH") {
		params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
		params.fluid_viscosity = 0.02;
		mcparams.ours = false;
		mcparams.sparse = false;
		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.0);
		schedule.push_back(0.2);
		schedule.push_back(0.4);
		schedule.push_back(1.5);
		emitter_sp.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter_sp);

		WCSPH::SPH sph(false, false, true, "../res/emitter_sphere_obstacle_sph_surface_ob_st/emitter_sphere_obstacle_sph_surface_ob_st_", params, mcparams, emitters);
		sph.load_geometry(true, boundary_size_s, boundary_left_s, fluid_sizes, fluid_lefts, std::pair <WCSPH::Vector, WCSPH::Vector>(), 0.18);
		sph.simulate(2.0);
	}

	SECTION("Adhesive sphere obstacle with surface tension in SPH") {
		params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
		params.fluid_viscosity = 0.02;
		params.cohesion_coefficient = 0.5;
		params.adhesion_coefficient = 3.0;
		mcparams.ours = false;
		mcparams.sparse = false;
		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.0);
		schedule.push_back(1.5);
		schedule.push_back(1.6);
		schedule.push_back(2.0);
		emitter_sp.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter_sp);

		WCSPH::SPH sph(false, false, true, "../res/emitter_adhesive_sphere_obstacle_sph_surface_ob_st/emitter_adhesive_sphere_obstacle_sph_surface_ob_st_", params, mcparams, emitters);
		sph.load_geometry(true, boundary_size_s, boundary_left_s, fluid_sizes, fluid_lefts, std::pair <WCSPH::Vector, WCSPH::Vector>(), 0.18);
		sph.simulate(2.0);
	}

	SECTION("Adhesive sphere obstacle with surface tension in SPH") {
		params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
		params.fluid_viscosity = 0.005;
		params.cohesion_coefficient = 0.5;
		params.adhesion_coefficient = 3.0;
		mcparams.ours = false;
		mcparams.sparse = false;
		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.0);
		schedule.push_back(1.5);
		schedule.push_back(1.6);
		schedule.push_back(2.0);
		emitter_sp.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter_sp);

		WCSPH::SPH sph(false, false, true, "../res/emitter_adhesive_sphere_obstacle_less_viscous_sph_surface_ob_st/emitter_adhesive_sphere_obstacle_less_viscous_sph_surface_ob_st_", params, mcparams, emitters);
		sph.load_geometry(true, boundary_size_s, boundary_left_s, fluid_sizes, fluid_lefts, std::pair <WCSPH::Vector, WCSPH::Vector>(), 0.18);
		sph.simulate(2.0);
	}

	SECTION("Adhesive sphere obstacle with surface tension in SPH") {
		params.export_type = SPHParameters::EXPORT_WITH_SURFACE;
		params.fluid_viscosity = 0.005;
		mcparams.ours = false;
		mcparams.sparse = false;
		std::vector<CompactNSearch::Real> schedule;
		schedule.push_back(0.0);
		schedule.push_back(1.5);
		schedule.push_back(1.6);
		schedule.push_back(2.0);
		emitter_sp.set_schedule(schedule);
		std::vector<Emitter::Emitter> emitters;
		emitters.push_back(emitter_sp);

		WCSPH::SPH sph(false, false, true, "../res/emitter_sphere_obstacle_less_viscous_sph_surface_ob_st/emitter_sphere_obstacle_less_viscous_sph_surface_ob_st_", params, mcparams, emitters);
		sph.load_geometry(true, boundary_size_s, boundary_left_s, fluid_sizes, fluid_lefts, std::pair <WCSPH::Vector, WCSPH::Vector>(), 0.18);
		sph.simulate(2.0);
	}
}*/

/*TEST_CASE("Emitter Vertical", "[Emitter Vertical]") {
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
	params.emit_frequency = 1.0;
	params.max_num_particles = 100000;

	Emitter::Emitter emitter(0.0625, { 1.0, 0.0, 0.0 }, 0.0625, { 0.0, 0.0, 1.0 }, { 1.5, 0.0, 1.5 }, { 0.0, 1.0, 0.0 }, { 0.0, 3.0, 0.0 }, params.particle_diameter, params.particle_diameter);

	std::vector<CompactNSearch::Real> schedule;
	schedule.push_back(0.2);
	schedule.push_back(2.1);
	emitter.set_schedule(schedule);
	std::vector<Emitter::Emitter> emitters;
	emitters.push_back(emitter);

	WCSPH::SPH sph(false, false, "../res/emitter_sph_vertical/emitter_sph_vertical_", params, mcparams, emitters);
	sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	sph.simulate(1);
}*/

/*TEST_CASE("Multiple Emitters", "[Multiple Emitters]") {
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
	params.emit_frequency = 1.0;
	params.max_num_particles = 100000;

	Emitter::Emitter emitter1(0.0625, { 1.0, 0.0, 0.0 }, 0.0625, { 0.0, 0.0, 1.0 }, { 1.5, 0.0, 1.5 }, { 0.0, 1.0, 0.0 }, { 0.0, 3.0, 0.0 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter2(0.0625, { 0.0, 1.0, 0.0 }, 0.0625, { 0.0, 0.0, 1.0 }, { 0.0, 1.5, 1.5 }, { 1.0, 0.0, 0.0 }, { 3.0, 0.0, 0.0 }, params.particle_diameter, params.particle_diameter);

	std::vector<CompactNSearch::Real> schedule;
	schedule.push_back(0.2);
	schedule.push_back(2.1);
	emitter1.set_schedule(schedule);
	emitter2.set_schedule(schedule);
	std::vector<Emitter::Emitter> emitters;
	emitters.push_back(emitter1);
	emitters.push_back(emitter2);

	WCSPH::SPH sph(false, false, "../res/emitter_sph_multiple/emitter_sph_multiple_", params, mcparams, emitters);
	sph.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	sph.simulate(1);
}*/

TEST_CASE("Test Emitter Fountain PBF", "[Test Emitter Fountain PBF]")
{
	PBD::Vector boundary_size = { 3.0, 3.0, 3.0 };
	PBD::Vector boundary_left = { -0.015, -0.015, -0.015 };
	std::vector<PBD::Vector> fluid_sizes;

	std::vector<PBD::Vector> fluid_lefts;

	PBFParameters params;
	MCParameters mcparams;
	params.dt = 0.00025;
	params.dt_next_frame = 0.01;
	params.particle_radius = 0.005;
	params.fluid_rest_density = 1000;
	params.max_dt = 0.0005;
	params.max_velocity_cap = 5;
	params.fluid_pressure_stiffness = 1000.0;
	params.fluid_viscosity = 0.005;
	params.boundary_viscosity = 0.0;

	params.particle_diameter = 2 * params.particle_radius;
	params.fluid_sampling_distance = params.particle_diameter;
	params.boundary_sampling_distance = 0.8 * params.particle_diameter;
	params.smoothing_length = 1.2 * params.particle_diameter;
	params.smoothing_length_squared = params.smoothing_length * params.smoothing_length;
	params.compact_support = 2 * params.smoothing_length;
	params.grid_cell_size = 1.2 * params.particle_radius;
	params.emit_frequency = 1.2;
	params.max_num_particles = 100000;

	params.export_type = PBFParameters::EXPORT_WITH_SURFACE;
	mcparams.ours = false;
	mcparams.sparse = false;

	//Emitter::Emitter emitter1(0.0625, { 1.0, 0.0, 0.0 }, 0.0625, { 0.0, 0.0, 1.0 }, { 1.5, 0.0, 1.5 }, { 0.0, 1.0, 0.0 }, { 0.0, 3.0, 0.0 }, params.particle_diameter, params.particle_diameter);
	//Emitter::Emitter emitter2(0.0625, { 0.0, 1.0, 0.0 }, 0.0625, { 0.0, 0.0, 1.0 }, { 0.0, 1.5, 1.5 }, { 1.0, 0.0, 0.0 }, { 3.0, 0.0, 0.0 }, params.particle_diameter, params.particle_diameter);
	Emitter::Emitter emitter(0.0625, { 1.0, 0.0, 0.0 }, 0.0625, { 0.0, 1.0, 0.0 }, { 1.5, 1.5, 0.1 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 5.0 }, params.particle_diameter, params.particle_diameter);

	std::vector<CompactNSearch::Real> schedule;
	schedule.push_back(0.0);
	schedule.push_back(2.1);
	emitter.set_schedule(schedule);
	std::vector<Emitter::Emitter> emitters;
	emitters.push_back(emitter);

	PBD::PBF pbf(false, false, "../res/emitter_fountain_pbf/emitter_fountain_pbf_", params, mcparams, emitters);
	pbf.load_geometry(true, boundary_size, boundary_left, fluid_sizes, fluid_lefts);
	pbf.simulate(2);
}