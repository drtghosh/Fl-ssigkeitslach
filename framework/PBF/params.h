#pragma once
#include <CompactNSearch/Config.h>

class PBFParameters
{
public:
	enum export_type { EXPORT, EXPORT_WITH_SURFACE };

	unsigned int num_threads = 4;
	unsigned int pbf_iterations = 5;
	CompactNSearch::Real particle_radius = 0.03;
	CompactNSearch::Real fluid_rest_density = 1000;
	CompactNSearch::Real dt = 0.1;
	CompactNSearch::Real dt_next_frame = dt;
	CompactNSearch::Real max_velocity_cap = 1;
	CompactNSearch::Real max_dt = 0.01;
	CompactNSearch::Real fluid_viscosity = 0.0;
	CompactNSearch::Real fluid_pressure_stiffness = 0.0;
	CompactNSearch::Real boundary_viscosity = 0.0;

	CompactNSearch::Real particle_diameter = 2 * particle_radius;
	CompactNSearch::Real fluid_sampling_distance = particle_diameter;
	CompactNSearch::Real boundary_sampling_distance = 0.8 * particle_diameter;
	CompactNSearch::Real smoothing_length = 1.2 * particle_diameter;
	CompactNSearch::Real smoothing_length_squared = smoothing_length * smoothing_length;
	CompactNSearch::Real compact_support = 2 * smoothing_length;

	export_type export_type = export_type::EXPORT;
	bool sparse = true;
	bool ours = false;

	CompactNSearch::Real grid_cell_size = 1.2 * particle_radius;
	CompactNSearch::Real initial_grid_value = 0.55;
};
