#pragma once

#include "config.h"
#include <Eigen/Dense>
#include "../kernel/kernel.h"
#include "params.h"
#include <CompactNSearch/CompactNSearch.h>
#include "../mcubes/marching_cubes.h"
#include <array>

namespace WCSPH
{
	class SPH
	{
	private:
		// Basic parameters
		Parameters parameters;
		MCParameters mcparameters;
		std::string result_path;
		kernel::Kernel kernel{ 1.0 };
		MC::MarchingCubes marching_cubes{ mcparameters };

		int timesteps;
		int duration = 0;
		//CompactNSearch::NeighborhoodSearch cns{ parameters.compact_support };
		CompactNSearch::NeighborhoodSearch cns_boundary{ 1.0 };
		CompactNSearch::NeighborhoodSearch cns{ 1.0 };
		bool gravity_only{ false };
		bool with_initial_velocity{ false };

		// Particle data for fluid
		std::vector<WCSPH::Vector> fluid_particles;
		std::vector<WCSPH::Vector> fluid_velocities;
		std::vector<WCSPH::Vector> fluid_accelerations;
		std::vector<WCSPH::Vector> pressure_accelerations;
		std::vector<WCSPH::Vector> velocity_accelerations;
		std::vector<CompactNSearch::Real> fluid_densities;
		std::vector<CompactNSearch::Real> fluid_pressures;
		std::vector<CompactNSearch::Real> fluid_particle_counts;
		WCSPH::Vector relative_velocity = { 2.0, 0.0, 0.0 };

		CompactNSearch::Real particle_mass;
		WCSPH::Vector boundary_box_bottom = WCSPH::Vector(0.0, 0.0, 0.0);
		WCSPH::Vector boundary_box_top = WCSPH::Vector(1.0, 1.0, 1.0);
		bool has_boundary{ true };

		// Particle data for boundary
		CompactNSearch::Real boundary_particle_volume;
		CompactNSearch::Real boundary_particle_mass;
		std::vector<WCSPH::Vector> boundary_particles;
		std::vector<CompactNSearch::Real> boundary_masses;
		std::vector<CompactNSearch::Real> boundary_volumes;
		std::vector<WCSPH::Vector> boundary_mesh_vertices;
		std::vector<std::array<int, 3>> boundary_mesh_faces;

		std::vector<CompactNSearch::Real> boundary_densities;
		std::vector<CompactNSearch::Real> boundary_pressure;


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
		void calculate_pressure_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id);
		void calculate_viscosity_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id);
		void calculate_other_acceleration();
		void update_time_step_size();
		void check_particle_positions();
		void export_data(unsigned int frame);
		void export_data_with_surface(unsigned int frame);

		void create_grid();
		void reset_grid_values();
		void update_grid_values(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid);
		void calculate_mc_normals();


	public:
		// Constructor
		SPH(bool gravity_only, bool with_initial_velocity, std::string result_path, Parameters params, MCParameters mcparams);

		// Destructor
		virtual ~SPH() = default;

		// Functions
		void simulate(CompactNSearch::Real t_end);
		void load_geometry(bool has_boundary, WCSPH::Vector& boundary_size, WCSPH::Vector& bottom_left_boundary, std::vector<WCSPH::Vector>& fluid_sizes, std::vector<WCSPH::Vector>& bottom_lefts_fluid);
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