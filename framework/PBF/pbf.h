#pragma once

#include "config.h"
#include <Eigen/Dense>
#include "../kernel/kernel.h"
#include "params.h"
#include <CompactNSearch/CompactNSearch.h>
#include "../mcubes/marching_cubes.h"
#include <array>
#include "../geometry/emitter.h"

namespace PBD
{
	class PBF
	{
	private:
		// Basic parameters
		PBFParameters parameters;
		MCParameters mcparameters;
		std::string result_path;
		kernel::Kernel kernel{ 1.0 };
		MC::MarchingCubes marching_cubes{ mcparameters };

		int timesteps;
		int duration = 0;
		CompactNSearch::NeighborhoodSearch cns_boundary{ 1.0 };
		CompactNSearch::NeighborhoodSearch cns{ 1.0 };
		bool gravity_only{ false };
		bool with_initial_velocity{ false };
		bool with_surface_tension{ false };
		std::vector<Emitter::Emitter> emitter;

		// Particle data for fluid
		std::vector<PBD::Vector> fluid_particles;
		std::vector<PBD::Vector> old_fluid_particles;
		std::vector<PBD::Vector> fluid_velocities;
		std::vector<PBD::Vector> fluid_accelerations;
		std::vector<PBD::Vector> fluid_normals;
		std::vector<PBD::Vector> pressure_accelerations;
		std::vector<PBD::Vector> velocity_accelerations;
		std::vector<CompactNSearch::Real> fluid_densities;
		std::vector<CompactNSearch::Real> fluid_pressures;
		std::vector<CompactNSearch::Real> fluid_particle_counts;
		PBD::Vector relative_velocity = { 2.0, 0.0, 0.0 };
		std::vector<bool> is_fluid_particle_active;

		CompactNSearch::Real particle_mass;

		// Boundary data
		PBD::Vector boundary_box_bottom = PBD::Vector(0.0, 0.0, 0.0);
		PBD::Vector boundary_box_top = PBD::Vector(1.0, 1.0, 1.0);
		bool has_boundary{ true };
		bool open_boundary{ false };

		// Obstacle data
		std::pair <PBD::Vector, PBD::Vector>& obstacle_data = std::pair <PBD::Vector, PBD::Vector>();
		bool has_obstacle{ false };
		CompactNSearch::Real obstacle_sphere_radius = 0.0;
		bool has_obstacle_sphere{ false };

		// Particle data for boundary
		CompactNSearch::Real boundary_particle_volume;
		CompactNSearch::Real boundary_particle_mass;
		std::vector<PBD::Vector> boundary_particles;
		std::vector<CompactNSearch::Real> boundary_masses;
		std::vector<CompactNSearch::Real> boundary_volumes;
		std::vector<PBD::Vector> boundary_mesh_vertices;
		std::vector<std::array<int, 3>> boundary_mesh_faces;

		//Export tri mesh
		std::vector<PBD::Vector> boundary_mesh_vertices_export;
		std::vector<std::array<int, 3>> boundary_mesh_faces_export;
		std::vector<PBD::Vector> obstacle_mesh_vertices_export;
		std::vector<std::array<int, 3>> obstacle_mesh_faces_export;
		std::vector<PBD::Vector>  obstacle_sphere_mesh_vertices_export;
		std::vector<std::array<int, 3>> obstacle_sphere_mesh_faces_export;

		std::vector<CompactNSearch::Real> boundary_densities;
		std::vector<CompactNSearch::Real> boundary_pressure;


		//Other parameters:
		PBD::Vector gravity_acceleration = { 0.0, 0.0, -9.81 }; // -z direction as in main.cpp

		// Grid points data
		std::vector<PBD::Vector> grid_points;
		std::array<int,3> grid_resolution;
		std::vector<CompactNSearch::Real> grid_values;
		std::vector<PBD::Vector> surface_vertices;
		std::vector<std::array<int, 3>> surface_triangles;
		std::vector<PBD::Vector> surface_normals;
		std::unordered_map<uint64_t, CompactNSearch::Real> grid_values_map_sparse;

		// Functions
		void initialize();
		void compute_boundary_mass(unsigned int boundary_id, CompactNSearch::PointSet const& ps_boundary);
		void calculate_particle_density(unsigned int fluid_id, CompactNSearch::PointSet const& ps_fluid, unsigned int boundary_id);
		void semi_implicit_euler();
		void calculate_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id);
		//void calculate_pressure(bool first_step);
		void calculate_fluid_normals(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid);
		//void calculate_pressure_acceleration(unsigned int fluid_id, CompactNSearch::PointSet const& ps_fluid, unsigned int boundary_id);
		void calculate_viscosity_acceleration(unsigned int fluid_id, CompactNSearch::PointSet const& ps_fluid, unsigned int boundary_id);
		void calculate_other_acceleration();
		void calculate_st_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id);
		void calculate_adhesion_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id);
		void update_time_step_size();
		void check_particle_positions();
		void update_particle_positions(unsigned int fluid_id, CompactNSearch::PointSet const& ps_fluid, unsigned int boundary_id);
		void update_particle_velocities();
		void export_data(unsigned int frame);
		void export_data_with_surface(unsigned int frame);
		void export_boundary();
		void export_obstacles();
		void emit_particles(CompactNSearch::Real t_sim);
		void update_active_status();

		void create_grid();
		void reset_grid_values();
		void update_grid_values(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid);
		void calculate_mc_normals();


	public:
		// Constructor
		PBF(bool gravity_only, bool with_initial_velocity, bool with_surface_tension, std::string result_path, PBFParameters params, MCParameters mcparams, std::vector<Emitter::Emitter> emitter = std::vector<Emitter::Emitter>());

		// Destructor
		virtual ~PBF() = default;

		// Functions
		void simulate(CompactNSearch::Real t_end);
		void load_geometry(bool has_boundary, bool open_boundary, PBD::Vector & boundary_size, PBD::Vector & bottom_left_boundary, std::vector<PBD::Vector>&fluid_sizes, std::vector<PBD::Vector>&bottom_lefts_fluid, std::pair <PBD::Vector, PBD::Vector>&obstacle_box = std::pair <PBD::Vector, PBD::Vector>({ 0,0,0 }, {0,0,0}), CompactNSearch::Real obstacle_sphere_radius = 0.0);
		void turn_off_gravity();
		void turn_on_gravity();
		void add_initial_velocity(std::vector<PBD::Vector> velocities);
		void add_relative_velocity(PBD::Vector relative_velocity);
		int num_boundary_particles();
		int num_fluid_particles();
		int num_timesteps();
		void printStats(std::string test);
	};
}