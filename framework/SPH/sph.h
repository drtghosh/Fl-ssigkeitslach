#pragma once

#include "config.h"
#include <Eigen/Dense>
#include "../kernel/kernel.h"
#include "params.h"
#include <CompactNSearch/CompactNSearch.h>
#include "../mcubes/marching_cubes.h"
#include <array>
#include "../geometry/emitter.h"

namespace WCSPH
{
	class SPH
	{
	private:
		// Basic parameters
		SPHParameters parameters;
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
		std::vector<WCSPH::Vector> fluid_particles;
		std::vector<WCSPH::Vector> fluid_velocities;
		std::vector<WCSPH::Vector> fluid_accelerations;
		std::vector<WCSPH::Vector> fluid_normals;
		std::vector<WCSPH::Vector> pressure_accelerations;
		std::vector<WCSPH::Vector> velocity_accelerations;
		std::vector<CompactNSearch::Real> fluid_densities;
		std::vector<CompactNSearch::Real> fluid_pressures;
		std::vector<CompactNSearch::Real> fluid_particle_counts;
		WCSPH::Vector relative_velocity = { 2.0, 0.0, 0.0 };
		std::vector<bool> is_fluid_particle_active;

		CompactNSearch::Real particle_mass;

		// Boundary data
		bool has_boundary{ true };

		// Particle data for boundary
		std::vector<WCSPH::Vector> boundary_particles;
		std::vector<CompactNSearch::Real> boundary_masses;
		std::vector<CompactNSearch::Real> boundary_volumes;
		std::vector<WCSPH::Vector> boundary_mesh_vertices;
		std::vector<std::array<int, 3>> boundary_mesh_faces;
		std::vector<CompactNSearch::Real> boundary_densities;
		std::vector<CompactNSearch::Real> boundary_pressure;

		//Export tri mesh
		std::vector<WCSPH::Vector> boundary_mesh_vertices_export;
		std::vector<std::array<int, 3>> boundary_mesh_faces_export;


		//Other parameters:
		WCSPH::Vector gravity_acceleration = { 0.0, 0.0, -9.81 }; // -z direction as in main.cpp

		// Grid points data
		std::vector<WCSPH::Vector> grid_points;
		std::array<int,3> grid_resolution;
		std::vector<CompactNSearch::Real> grid_values;
		std::vector<WCSPH::Vector> surface_vertices;
		std::vector<std::array<int, 3>> surface_triangles;
		std::vector<WCSPH::Vector> surface_normals;
		std::unordered_map<uint64_t, CompactNSearch::Real> grid_values_map_sparse;

		// Functions
		void initialize();
		void compute_boundary_mass(unsigned int boundary_id, CompactNSearch::PointSet const& ps_boundary);
		void calculate_particle_density(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id);
		void semi_implicit_euler();
		void calculate_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id);
		void calculate_pressure(bool first_step);
		void calculate_fluid_normals(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid);
		void calculate_pressure_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id);
		void calculate_viscosity_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id);
		void calculate_other_acceleration();
		void calculate_st_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id);
		void calculate_adhesion_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id);
		void update_time_step_size();
		void export_data(unsigned int frame);
		void export_data_with_surface(unsigned int frame);
		void export_boundary();
		void emit_particles(CompactNSearch::Real t_sim);
		void update_active_status();

		void create_grid();
		void reset_grid_values();
		void update_grid_values(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid);
		void calculate_mc_normals();


	public:
		// Constructor
		SPH(bool gravity_only, bool with_initial_velocity, bool with_surface_tension, std::string result_path, SPHParameters params, MCParameters mcparams, std::vector<Emitter::Emitter>& emitter = std::vector<Emitter::Emitter>());

		// Destructor
		virtual ~SPH() = default;

		// Functions
		void simulate(CompactNSearch::Real t_end);
		void load_geometry(std::vector<WCSPH::Vector>& fluid_sizes, std::vector<WCSPH::Vector>& bottom_lefts_fluid);
		void turn_off_gravity();
		void turn_on_gravity();
		void add_initial_velocity(std::vector<WCSPH::Vector> velocities);
		void add_relative_velocity(WCSPH::Vector relative_velocity);
		int num_boundary_particles();
		int num_fluid_particles();
		int num_timesteps();
		void printStats(std::string test);
	};
}