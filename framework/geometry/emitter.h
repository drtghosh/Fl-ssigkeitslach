#pragma once

#include <Eigen/Dense>
#include "CompactNSearch/Config.h"
#include <array>
#include "config.h"

#define USE_DOUBLE //TODO sort out configs

namespace Emitter
{
	class Emitter
	{
	private:
		//Emitter shap params
		CompactNSearch::Real r1;
		geometry::Vector r1_dir;
		CompactNSearch::Real r2;
		geometry::Vector r2_dir;
		geometry::Vector centre_position;
		geometry::Vector normal;

		const CompactNSearch::Real distance;
		const CompactNSearch::Real velocity;

		std::vector<CompactNSearch::Real> schedule; //times when emitters should be turned on and off, starts at 0.0 with off automatically (first entry is time when turned on)
		std::vector<std::array<int, 3>> shield_triangles;
		std::vector<geometry::Vector> shield_vertices;
		void calculate_shield();

	public:
		Emitter(CompactNSearch::Real radius1, geometry::Vector radius1_dir, CompactNSearch::Real radius2, geometry::Vector radius2_dir, geometry::Vector centre, geometry::Vector normal_dir, const CompactNSearch::Real velocity, const CompactNSearch::Real shield_distance);

		void set_schedule(std::vector<CompactNSearch::Real> schedule_times);
		std::vector<std::array<int, 3 >> get_shield_triangles();
		std::vector<geometry::Vector> get_shield_vertices();
		void emit_particles(std::vector<geometry::Vector>& particle_positions, std::vector<geometry::Vector>& particle_velocities);
	};
}