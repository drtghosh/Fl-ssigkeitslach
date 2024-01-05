#include "sph.h"
#include "../kernel/kernel.h"
#include "../CompactNSearch/include/CompactNSearch/CompactNSearch.h"
#include "../geometry/io.h"
#include "../geometry/sampling.h"

#include <iostream>
#include <omp.h>
#include <chrono>


#define USE_Double

namespace WCSPH
{

	SPH::SPH(bool gravity_only, bool with_initial_velocity, std::string result_path, SPHParameters params, MCParameters mcparams, std::vector<Emitter::Emitter>& emitter) : gravity_only(gravity_only), with_initial_velocity(with_initial_velocity), result_path(result_path), parameters(params), mcparameters(mcparams), emitter(emitter) {
		auto start = std::chrono::system_clock::now();

		this->kernel = kernel::Kernel(parameters.smoothing_length);
		this->cns_boundary = CompactNSearch::NeighborhoodSearch(parameters.compact_support);
		this->cns = CompactNSearch::NeighborhoodSearch(parameters.compact_support);
		this->marching_cubes = MC::MarchingCubes(mcparams, 0.0);

		auto end = std::chrono::system_clock::now();
		auto elapsed =
			std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		duration += elapsed.count();
	}

	void SPH::simulate(CompactNSearch::Real t_end) {
		auto start = std::chrono::system_clock::now();

		initialize();

		// Add initial velocity if needed
		if (with_initial_velocity) {
			add_relative_velocity(this->relative_velocity);
		}

		// Compute boundary mass and add boundary particles to cns
		unsigned int boundary_particles_id = 0;
		if (has_boundary) {
			boundary_particles_id = this->cns_boundary.add_point_set(boundary_particles.front().data(), boundary_particles.size());
			this->cns_boundary.find_neighbors();

			CompactNSearch::PointSet const& pointset_boundary = this->cns_boundary.point_set(boundary_particles_id);
			compute_boundary_mass(boundary_particles_id, pointset_boundary);

			boundary_particles_id = this->cns.add_point_set(boundary_particles.front().data(), boundary_particles.size());
		}

		//Add fluid particles to cns
		unsigned int fluid_particles_id = -1;
		if (fluid_particles.size() > 0) {
			fluid_particles_id = this->cns.add_point_set(fluid_particles.front().data(), fluid_particles.size());
			CompactNSearch::PointSet& pointset_fluid = this->cns.point_set(fluid_particles_id);
			this->cns.find_neighbors();
		}

		// Needed later for particle density calculation
		//unsigned int boundary_particles_id_fluid = this->cns.add_point_set(boundary_particles.front().data(), boundary_particles.size());

		// Simulation loop
		CompactNSearch::Real t_sim = 0.0;
		CompactNSearch::Real t_next_frame = 0.0;
		unsigned int frame_idx = 0;
		CompactNSearch::Real cx = 0.0;
		timesteps = 0;

		export_boundary();
		export_obstacles();

		while (t_sim < t_end) {
			if (fluid_particles.size() <= parameters.max_num_particles) {
				emit_particles(t_sim);
			}

			if (this->fluid_particles.size() == 0) {
				if (t_sim >= t_next_frame) {
					switch (parameters.export_type) {
					case SPHParameters::export_type::EXPORT:
						export_data(frame_idx);
						break;
					case SPHParameters::export_type::EXPORT_WITH_SURFACE:
						export_data_with_surface(frame_idx);
						break;
					default:
						export_data(frame_idx);
						break;
					}
					t_next_frame += parameters.dt_next_frame;
					frame_idx++;
				}
				// Timestep update
				t_sim += parameters.dt;
				timesteps += 1;
				continue;
			}
			// Update dt
			update_time_step_size();

			// Neighborhood Search
			//CompactNSearch::NeighborhoodSearch cns{ parameters.compact_support };
			//unsigned int boundary_particles_id = 0;
			//if (has_boundary) {
				//boundary_particles_id = this->cns.add_point_set(boundary_particles.front().data(), boundary_particles.size());
			//}
			//unsigned int fluid_particles_id = this->cns.add_point_set(fluid_particles.front().data(), fluid_particles.size());
			if (fluid_particles_id == -1) {
				fluid_particles_id = this->cns.add_point_set(fluid_particles.front().data(), fluid_particles.size());
			}
			else {
				this->cns.resize_point_set(fluid_particles_id, fluid_particles.front().data(), fluid_particles.size());
			}
			CompactNSearch::PointSet& pointset_fluid = this->cns.point_set(fluid_particles_id);
			this->cns.find_neighbors();

			if (!gravity_only) {
				// Compute fluid density
				calculate_particle_density(fluid_particles_id, pointset_fluid, boundary_particles_id);

				// Compute pressure
				calculate_pressure(t_sim == 0.0);
			}
			// Compute accelerations
			calculate_acceleration(fluid_particles_id, pointset_fluid, boundary_particles_id);

			// Time integration
			semi_implicit_euler();

			update_active_status();

			// Check particle positions
			if (has_boundary) {
				check_particle_positions();
			}
			//Export VTK
			if (t_sim >= t_next_frame) {
				if (parameters.export_type == SPHParameters::export_type::EXPORT_WITH_SURFACE) {
					create_grid();
					reset_grid_values();
					this->cns.resize_point_set(fluid_particles_id, fluid_particles.front().data(), fluid_particles.size());
					pointset_fluid = this->cns.point_set(fluid_particles_id);
					this->cns.find_neighbors();
					update_grid_values(fluid_particles_id, pointset_fluid);
					surface_vertices.clear();
					surface_triangles.clear();
					surface_normals.clear();
					marching_cubes.calculate(grid_points, grid_values, grid_values_map_sparse, grid_resolution, surface_vertices, surface_triangles);
					/*unsigned int vertices_id_mc = this->cns.add_point_set(surface_vertices.front().data(), surface_vertices.size());
					CompactNSearch::PointSet& pointset_vertices = this->cns.point_set(vertices_id_mc);
					this->cns.resize_point_set(fluid_particles_id, fluid_particles.front().data(), fluid_particles.size());
					pointset_fluid = this->cns.point_set(fluid_particles_id);
					this->cns.find_neighbors();*/
					calculate_mc_normals();
				}

				switch(parameters.export_type) {
					case SPHParameters::export_type::EXPORT : 
						export_data(frame_idx);
						break;
					case SPHParameters::export_type::EXPORT_WITH_SURFACE:
						export_data_with_surface(frame_idx);
						break;
					default:
						export_data(frame_idx);
						break;
				}
				t_next_frame += parameters.dt_next_frame;
				frame_idx++;
			}

			// Timestep update
			t_sim += parameters.dt;
			timesteps += 1;
		}

		auto end = std::chrono::system_clock::now();
		auto elapsed =
			std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		duration += elapsed.count();

	}

	void SPH::initialize() {
		this->fluid_accelerations.resize(fluid_particles.size());
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < (int)this->fluid_accelerations.size(); i++) {
			this->fluid_accelerations[i] = { 0.0,0.0,0.0 };
		}
		this->is_fluid_particle_active.resize(fluid_particles.size());
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < (int)this->is_fluid_particle_active.size(); i++) {
			this->is_fluid_particle_active[i] = true;
		}
		this->fluid_velocities.resize(fluid_particles.size());
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < (int)this->fluid_velocities.size(); i++) {
			this->fluid_velocities[i] = { 0.0,0.0,0.0 };
		}
		this->fluid_densities.resize(fluid_particles.size());
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < (int)this->fluid_densities.size(); i++) {
			this->fluid_densities[i] = 0.0;
		}
		this->fluid_pressures.resize(fluid_particles.size());
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < (int)this->fluid_pressures.size(); i++) {
			this->fluid_pressures[i] = 0.0;
		}
	}

	void SPH::load_geometry(bool has_boundary, WCSPH::Vector& boundary_size, WCSPH::Vector& bottom_left_boundary, std::vector<WCSPH::Vector>& fluid_sizes, std::vector<WCSPH::Vector>& bottom_lefts_fluid, std::vector<std::array<WCSPH::Vector, 4>>& obstacle_squares) {
		auto start = std::chrono::system_clock::now();
		
		this->has_boundary = has_boundary;
		if (has_boundary) {
			this->boundary_box_bottom = bottom_left_boundary;
			this->boundary_box_top = bottom_left_boundary + boundary_size;
			std::vector<WCSPH::Vector> vertices = { bottom_left_boundary, bottom_left_boundary + WCSPH::Vector(0.0, 0.0, boundary_size[2]), bottom_left_boundary + WCSPH::Vector(0.0, boundary_size[1], 0.0),
				bottom_left_boundary + WCSPH::Vector(0.0, boundary_size[1], boundary_size[2]), bottom_left_boundary + WCSPH::Vector(boundary_size[0], 0.0, 0.0), bottom_left_boundary + WCSPH::Vector(boundary_size[0], 0.0, boundary_size[2]),
				bottom_left_boundary + WCSPH::Vector(boundary_size[0], boundary_size[1], 0.0), bottom_left_boundary + boundary_size };

			std::vector<std::array<int, 3>> triangles = { {1, 2, 0}, {3, 6, 2}, {7, 4, 6}, {5, 0, 4}, {6, 0, 2}, {3, 5, 7},
														{1, 3, 2}, {3, 7, 6}, {7, 5, 4}, {5, 1, 0}, {6, 4, 0}, {3, 1, 5} };

			int vertices_before = vertices.size();

			boundary_mesh_faces_export = triangles;
			boundary_mesh_vertices_export = vertices;

			if (emitter.size() > 0) {
				for (unsigned int i = 0; i < emitter.size();i++) {
					std::vector<WCSPH::Vector> shield_vertices = emitter[i].get_shield_vertices();
					std::vector<std::array<int, 3>> shield_triangles = emitter[i].get_shield_triangles();
					for (unsigned int j = 0; j < shield_vertices.size(); j++) {
						vertices.push_back(shield_vertices[j]);
					}
					for (unsigned int j = 0; j < shield_triangles.size(); j++) {
						std::array<int, 3> triangle = shield_triangles[j];
						triangle[0] += vertices_before;
						triangle[1] += vertices_before;
						triangle[2] += vertices_before;
						triangles.push_back(triangle);
					}
				}
				
			}

			vertices_before = vertices.size();

			if (obstacle_squares.size() > 0) {
				std::vector<WCSPH::Vector> obstacle_vertices;
				std::vector<std::array<int, 3>> obstacle_triangles;

				for (int i = 0; i < obstacle_squares.size(); i++) {
					std::array<WCSPH::Vector, 4> square = obstacle_squares[i];
					for (int j = 0; j < square.size(); j++) {
						vertices.push_back(square[j]);
						obstacle_vertices.push_back(square[j]);
					}
					triangles.push_back({ vertices_before , vertices_before + 1 ,  vertices_before + 2 });
					triangles.push_back({ vertices_before + 2 , vertices_before + 3 ,  vertices_before });
					obstacle_triangles.push_back({ vertices_before , vertices_before + 1 ,  vertices_before + 2 });
					obstacle_triangles.push_back({ vertices_before + 2 , vertices_before + 3 ,  vertices_before });
					vertices_before += 4;
				}
				obstacle_mesh_vertices_export = obstacle_vertices;
				obstacle_mesh_faces_export = obstacle_triangles;
			}

			this->boundary_mesh_vertices = vertices;
			this->boundary_mesh_faces = triangles;
			geometry::sampling::triangle_mesh(this->boundary_particles, this->boundary_mesh_vertices, this->boundary_mesh_faces, this->parameters.boundary_sampling_distance);
			const CompactNSearch::Real boundary_volume = boundary_size[0] * boundary_size[1] * boundary_size[2];
			this->boundary_particle_volume = boundary_volume / this->boundary_particles.size();
			this->boundary_particle_mass = this->boundary_particle_volume * parameters.fluid_rest_density;
		}

		CompactNSearch::Real total_fluid_volume = 0.0;
		for (int i = 0; i < fluid_sizes.size(); ++i) {
			WCSPH::Vector top_right_fluid = bottom_lefts_fluid[i] + fluid_sizes[i];
			if (has_boundary) {
				WCSPH::Vector top_right_boundary = bottom_left_boundary + boundary_size;
				assert(bottom_lefts_fluid[i][0] >= bottom_left_boundary[0] && fluid_sizes[i][1] >= bottom_left_boundary[1] && fluid_sizes[i][2] >= bottom_left_boundary[2]);
				assert(top_right_fluid[0] <= top_right_boundary[0] && top_right_fluid[1] <= top_right_boundary[1] && top_right_fluid[2] <= top_right_boundary[2]);
			}
			geometry::sampling::fluid_box(this->fluid_particles, bottom_lefts_fluid[i], top_right_fluid, this->parameters.fluid_sampling_distance);
			total_fluid_volume += fluid_sizes[i][0] * fluid_sizes[i][1] * fluid_sizes[i][2];
			this->fluid_particle_counts.emplace_back(this->fluid_particles.size());
		}
		CompactNSearch::Real root = std::pow(parameters.max_num_particles, 1 / 3.) * parameters.particle_diameter;
		total_fluid_volume += root * root * root;
		this->particle_mass = (parameters.fluid_rest_density * total_fluid_volume) / (this->fluid_particles.size() + parameters.max_num_particles);

		auto end = std::chrono::system_clock::now();
		auto elapsed =
			std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		duration += elapsed.count();
	}

	void SPH::emit_particles(CompactNSearch::Real t_sim) {
		//Emit particles
		for (unsigned int i = 0; i < emitter.size(); i++) {
			std::vector<WCSPH::Vector> emitted_particle_positions;
			std::vector<WCSPH::Vector> emitted_particle_velocities;

			emitter[i].emit_particles(emitted_particle_positions, emitted_particle_velocities, t_sim);
			if (emitted_particle_positions.size() == 0) {
				continue;
			}
			bool emit = true;
			#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
			for (int j = 0; j < fluid_particles.size(); j++) {
				if (!is_fluid_particle_active[j]) {
					if ((fluid_particles[j] - emitted_particle_positions[0]).norm() <= parameters.emit_frequency * parameters.particle_diameter) {
						emit = false;
						break;
					}
				}
			}
			if (emit) {
				for (unsigned int j = 0; j < emitted_particle_positions.size(); j++) {
					fluid_particles.push_back(emitted_particle_positions[j]);
					fluid_velocities.push_back(emitted_particle_velocities[j]);
					fluid_accelerations.push_back({ 0.0,0.0,0.0 });
					is_fluid_particle_active.push_back(false);
				}
			}
		}
		this->fluid_densities.resize(fluid_particles.size());
		this->fluid_pressures.resize(fluid_particles.size());
	}

	void SPH::update_active_status() {
		//update active status
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < fluid_particles.size(); i++) {
			if (!is_fluid_particle_active[i]) {
				WCSPH::Vector emission_position;
				WCSPH::Vector dir = fluid_velocities[i].normalized();
				//Find origin position
				for (unsigned int j = 0; j < emitter.size(); j++) {
					std::vector<WCSPH::Vector> emission_positions = emitter[j].get_emition_positions();
					bool found = false;
					for (unsigned int k = 0; k < emission_positions.size(); k++) {
						WCSPH::Vector dist = (fluid_particles[i] - emission_positions[k]).normalized();
						found = (dist - dir).norm() < 0.0001;
						if (found) {
							emission_position = emission_positions[k];
							break;
						}
					}
					if (found) {
						break;
					}
				}
				bool is_active = (fluid_particles[i] - emission_position).norm() >= 3 * parameters.particle_diameter;
				if (is_active) {
					is_fluid_particle_active[i] = true;
				}
			}
		}
	}

	void SPH::compute_boundary_mass(unsigned int boundary_id, CompactNSearch::PointSet const& ps_boundary) {
		// std::vector<CompactNSearch::Real>& boundary_volumes not removed
		boundary_masses.resize(ps_boundary.n_points());
		boundary_volumes.resize(ps_boundary.n_points());

		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < ps_boundary.n_points(); ++i)
		{
			//Inititalize kernel sum
			CompactNSearch::Real kernel_sum = 0.0;
			CompactNSearch::Real density_sum = 0.0;
			// Get neighbors of particles
			for (int j = 0; j < ps_boundary.n_neighbors(boundary_id, i); ++j)
			{
				// Return the point id of the jth neighbor of the ith particle
				const unsigned int pid = ps_boundary.neighbor(boundary_id, i, j);
				const CompactNSearch::Real kernel_value = kernel.cubic_spline(boundary_particles.at(i), boundary_particles.at(pid));
				kernel_sum += kernel_value;
			}
			// Add contribution of particle itself
			kernel_sum += kernel.cubic_spline(boundary_particles.at(i), boundary_particles.at(i));
			boundary_masses[i] = this->parameters.fluid_rest_density / kernel_sum;
			boundary_volumes[i] = 1.0 / kernel_sum;
		}
	}

	void SPH::calculate_particle_density(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id) {
		// Loop over all fluid particles
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < ps_fluid.n_points(); ++i)
		{
			//Inititalize density
			CompactNSearch::Real density = 0.0;
			// Get neighbors of particles
			for (int j = 0; j < ps_fluid.n_neighbors(fluid_id, i); ++j)
			{
				// Return the point id of the jth neighbor of the ith particle
				const unsigned int fpid = ps_fluid.neighbor(fluid_id, i, j);
				density += kernel.cubic_spline(fluid_particles.at(i), fluid_particles.at(fpid));
			}
			// Add contribution of fluid particle itself
			density += kernel.cubic_spline(fluid_particles.at(i), fluid_particles.at(i));
			density *= this->particle_mass;
			// Get neighbors of boundary particles
			if (ps_fluid.n_neighbors(boundary_id, i) > 0 && has_boundary)
			{
				for (int j = 0; j < ps_fluid.n_neighbors(boundary_id, i); ++j)
				{
					// Return the point id of the jth neighbor of the ith particle
					const unsigned int bpid = ps_fluid.neighbor(boundary_id, i, j);
					density += (boundary_masses.at(bpid) * kernel.cubic_spline(fluid_particles.at(i), boundary_particles.at(bpid)));
				}
			}
			fluid_densities[i] = density;
		}
	}

	void SPH::semi_implicit_euler() {
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < (int)fluid_particles.size(); i++) {
				WCSPH::Vector new_v = fluid_velocities[i] + parameters.dt * fluid_accelerations[i];
				WCSPH::Vector new_x = fluid_particles[i] + parameters.dt * new_v;

				fluid_velocities[i] = new_v;
				fluid_particles[i] = new_x;
		}
	}

	void SPH::calculate_pressure(bool first_step) {
		CompactNSearch::Real subtractive_density = parameters.fluid_rest_density;
		if (first_step) {
			subtractive_density *= 1.02;
		}
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < fluid_pressures.size(); i++) {
			const CompactNSearch::Real pressure = std::max(parameters.fluid_pressure_stiffness * (fluid_densities.at(i) - subtractive_density), 0.0);
			fluid_pressures[i] = pressure;
		}
	}

	void SPH::calculate_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id) {
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < (int)fluid_accelerations.size(); i++) {
			fluid_accelerations[i] = { 0.0,0.0,0.0 };
		}

		if (!gravity_only) {
			calculate_pressure_acceleration(fluid_id, ps_fluid, boundary_id);
			calculate_viscosity_acceleration(fluid_id, ps_fluid, boundary_id);
		}

		calculate_other_acceleration();
	}

	void SPH::calculate_pressure_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id) {
		// Loop over all fluid particles
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < ps_fluid.n_points(); ++i) {
			if (is_fluid_particle_active[i]) {
				// Initialize acceleration
				WCSPH::Vector acceleration = { 0.0, 0.0, 0.0 };
				// Get fluid neighbors of particles
				for (int j = 0; j < ps_fluid.n_neighbors(fluid_id, i); ++j) {
					// Return the point id of the jth neighbor of the ith particle
					const unsigned int fpid = ps_fluid.neighbor(fluid_id, i, j);
					// Compute fluid-fluid interaction
					WCSPH::Vector ffa = -1.0 * ((fluid_pressures.at(i) / (fluid_densities.at(i) * fluid_densities.at(i))) + (fluid_pressures.at(fpid) / (fluid_densities.at(fpid) * fluid_densities.at(fpid)))) * kernel.cubic_grad_spline(fluid_particles.at(i), fluid_particles.at(fpid));
					acceleration += ffa;
				}
				acceleration *= this->particle_mass;
				// Get boundary neighbors of particles
				if (ps_fluid.n_neighbors(boundary_id, i) > 0 && has_boundary) {
					for (int j = 0; j < ps_fluid.n_neighbors(boundary_id, i); ++j) {
						// Return the point id of the jth neighbor of the ith particle
						const unsigned int bpid = ps_fluid.neighbor(boundary_id, i, j);
						// Compute fluid-boundary interaction
						WCSPH::Vector fba = -1.0 * parameters.fluid_rest_density * boundary_volumes.at(bpid) * (fluid_pressures.at(i) / (fluid_densities.at(i) * fluid_densities.at(i))) * kernel.cubic_grad_spline(fluid_particles.at(i), boundary_particles.at(bpid));
						acceleration += fba;
					}
				}
				fluid_accelerations[i] += acceleration;
			}
		}
	}

	void SPH::calculate_viscosity_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id) {
		// Loop over all fluid particles
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < ps_fluid.n_points(); ++i) {
			if (is_fluid_particle_active[i]) {
				// Initialize acceleration
				WCSPH::Vector acceleration = { 0.0, 0.0, 0.0 };
				// Get fluid neighbors of particles
				for (int j = 0; j < ps_fluid.n_neighbors(fluid_id, i); ++j) {
					// Return the point id of the jth neighbor of the ith particle
					const unsigned int fpid = ps_fluid.neighbor(fluid_id, i, j);
					// Compute fluid-fluid interaction
					const CompactNSearch::Real normff = (fluid_particles.at(i) - fluid_particles.at(fpid)).norm();
					WCSPH::Vector ffa = kernel.cubic_grad_spline(fluid_particles.at(i), fluid_particles.at(fpid)).dot(fluid_particles.at(i) - fluid_particles.at(fpid)) * (fluid_velocities.at(i) - fluid_velocities.at(fpid)) / (fluid_densities.at(fpid) * ((normff * normff) + (0.01 * parameters.smoothing_length_squared)));
					acceleration += ffa;
				}
				acceleration *= (2 * parameters.fluid_viscosity * this->particle_mass);
				// Get boundary neighbors of particles
				if (ps_fluid.n_neighbors(boundary_id, i) > 0 && has_boundary) {
					for (int j = 0; j < ps_fluid.n_neighbors(boundary_id, i); ++j) {
						// Return the point id of the jth neighbor of the ith particle
						const unsigned int bpid = ps_fluid.neighbor(boundary_id, i, j);
						// Compute fluid-boundary interaction
						const CompactNSearch::Real normfb = (fluid_particles.at(i) - boundary_particles.at(bpid)).norm();
						WCSPH::Vector fba = 2 * parameters.boundary_viscosity * boundary_volumes.at(bpid) * kernel.cubic_grad_spline(fluid_particles.at(i), boundary_particles.at(bpid)).dot(fluid_particles.at(i) - boundary_particles.at(bpid)) * fluid_velocities.at(i) / ((normfb * normfb) + (0.01 * parameters.smoothing_length_squared));
						acceleration += fba;
					}
				}
				fluid_accelerations[i] += acceleration;
			}
		}
	}

	void SPH::calculate_other_acceleration() {
		if (!gravity_only) {
			//TODO if we have other accelerations
			//if (is_fluid_particle_active[i]) {}
		}
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < (int)fluid_accelerations.size(); i++) {
			if (is_fluid_particle_active[i]) {
				fluid_accelerations[i] += gravity_acceleration;
			}
		}
	}

	void SPH::update_time_step_size() {
		CompactNSearch::Real max_norm = 0.0;
		for (size_t i = 0; i < fluid_velocities.size(); i++) {
			CompactNSearch::Real norm = fluid_velocities[i].norm();
			if (norm > max_norm) {
				max_norm = norm;
			}
			if (max_norm > parameters.max_velocity_cap) {
				max_norm = parameters.max_velocity_cap;
				break;
			}
		}

		CompactNSearch::Real dt_cfl = 0.5 * parameters.particle_radius / max_norm;
		parameters.dt = std::min(parameters.max_dt, dt_cfl);
	}

	void SPH::check_particle_positions() {
		int old_size = (int)fluid_particles.size();
		//learnSPH::Vector bottom_left= boundary_bounds.corner(boundary_bounds.BottomLeftFloor);
		//learnSPH::Vector top_right = boundary_bounds.corner(boundary_bounds.TopRightCeil); //TODO check
		std::vector<bool> new_is_fluid_particle_active;
		std::vector<WCSPH::Vector> new_positions;
		std::vector<WCSPH::Vector> new_velocities;
		std::vector<WCSPH::Vector> new_accelerations;
		std::vector<CompactNSearch::Real> new_densities;
		std::vector<CompactNSearch::Real> new_pressures;
		new_is_fluid_particle_active.reserve(old_size);
		new_positions.reserve(old_size);
		new_velocities.reserve(old_size);
		new_accelerations.reserve(old_size);
		new_densities.reserve(old_size);
		new_pressures.reserve(old_size);

		for (int i = 0; i < old_size; i++) {
			WCSPH::Vector pos = fluid_particles[i];
			bool is_bigger_BL = pos[0] >= this->boundary_box_bottom[0] && pos[1] >= this->boundary_box_bottom[1] && pos[2] >= this->boundary_box_bottom[2];
			bool is_smaller_TR = pos[0] <= this->boundary_box_top[0] && pos[1] <= this->boundary_box_top[1] && pos[2] <= this->boundary_box_top[2];
			bool is_inside = is_bigger_BL && is_smaller_TR;

			if (!is_inside) {
				continue;
			}
			new_is_fluid_particle_active.emplace_back(is_fluid_particle_active[i]);
			new_positions.emplace_back(pos);
			new_velocities.emplace_back(fluid_velocities[i]);
			new_accelerations.emplace_back(fluid_accelerations[i]);
			new_densities.emplace_back(fluid_densities[i]);
			new_pressures.emplace_back(fluid_pressures[i]);
		}
		this->is_fluid_particle_active = new_is_fluid_particle_active;
		this->fluid_particles = new_positions;
		this->fluid_velocities = new_velocities;
		this->fluid_accelerations = new_accelerations;
		this->fluid_densities = new_densities;
		this->fluid_pressures = new_pressures;

	}
	
	void SPH::turn_off_gravity() {
		this->gravity_acceleration = { 0.0, 0.0, 0.0 };
	}

	void SPH::turn_on_gravity() {
		this->gravity_acceleration = { 0.0, 0.0, -9.81 };
	}

	void SPH::add_initial_velocity(std::vector<WCSPH::Vector> velocities) {
		assert(velocities.size() == this->fluid_velocities.size());
		this->fluid_velocities = velocities;
	}

	void SPH::add_relative_velocity(WCSPH::Vector relative_velocity) {
		std::vector<WCSPH::Vector> velocities;
		for (int i = 0; i < (int)this->fluid_particle_counts[0]; i++) {
			velocities.emplace_back(0.5 * relative_velocity);
		}
		for (int i = (int)this->fluid_particle_counts[0]; i < (int)this->fluid_particle_counts[1]; i++) {
			velocities.emplace_back(-0.5 * relative_velocity);
		}
		assert(velocities.size() == this->fluid_velocities.size());
		this->fluid_velocities = velocities;
	}

	int SPH::num_boundary_particles() {
		return boundary_particles.size();
	}

	int SPH::num_fluid_particles() {
		int num_particles = 0;
		for (int i = 0;i < fluid_particle_counts.size();i++) {
			num_particles += fluid_particle_counts[i];
		}
		return num_particles;
	}

	int SPH::num_timesteps() {
		return timesteps;
	}

	void SPH::printStats(std::string test) {
		int num_boundary = boundary_particles.size();
		int num_fluid = 0;
		for (int i = 0;i < fluid_particle_counts.size();i++) {
			num_fluid += fluid_particle_counts[i];
		}

		std::cout << test << ": " << "(boundary: " << num_boundary << ", fluid: " << num_fluid << ", timesteps: " << timesteps << ")" << std::endl;
		std::cout << " => Time: " << duration << " ms" << std::endl;
	}

	void SPH::create_grid() {
		grid_points.clear();
		WCSPH::Vector bottom_left = { std::numeric_limits<CompactNSearch::Real>::max(), std::numeric_limits<CompactNSearch::Real>::max(), std::numeric_limits<CompactNSearch::Real>::max() };
		WCSPH::Vector top_right = { std::numeric_limits<CompactNSearch::Real>::min(), std::numeric_limits<CompactNSearch::Real>::min(), std::numeric_limits<CompactNSearch::Real>::min() };

		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < fluid_particles.size(); i++) {
			WCSPH::Vector pos = fluid_particles[i];
			for (int j = 0; j < 3; j++) {
				if (pos[j] < bottom_left[j]) {
					bottom_left[j] = pos[j];
				}
				if (pos[j] > top_right[j]) {
					top_right[j] = pos[j];
				}
			}
		}

		for (int k = 0; k < 3; k++) {
			bottom_left[k] -= parameters.compact_support;
			top_right[k] += parameters.compact_support;
		}

		for (int i = 0; i < 3; i++) {
			CompactNSearch::Real diff = top_right[i] - bottom_left[i];
			int res = (int)std::ceil(diff / parameters.grid_cell_size);
			grid_resolution[i]=res+1;
		}

		for (int i = 0; i < grid_resolution[0]; i++) {
			for (int j = 0; j < grid_resolution[1]; j++) {
				for (int k = 0; k < grid_resolution[2]; k++) {
					WCSPH::Vector new_point = { bottom_left[0] + i * parameters.grid_cell_size, bottom_left[1] + j * parameters.grid_cell_size, bottom_left[2] + k * parameters.grid_cell_size };
					grid_points.emplace_back(new_point);
				}
			}
		}
	}

	void SPH::reset_grid_values() {
		this->grid_values.resize(grid_points.size());
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < this->grid_values.size(); i++) {
			this->grid_values[i] = -1.0 * parameters.initial_grid_value;
		}
	}

	void SPH::update_grid_values(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid) {
		//CompactNSearch::NeighborhoodSearch cns_mc{ parameters.compact_support };
		//unsigned int fluid_id_mc = cns_mc.add_point_set(this->fluid_particles.front().data(), this->fluid_particles.size());
		//CompactNSearch::PointSet const& ps_fluid_mc = cns_mc.point_set(fluid_id_mc);
		//cns_mc.find_neighbors();
		// Loop over all fluid particles
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < ps_fluid.n_points(); ++i) {
			// Inititalize density
			CompactNSearch::Real normalized_density = 0.0;
			// Get neighbors among fluid particles
			for (int j = 0; j < ps_fluid.n_neighbors(fluid_id, i); ++j) {
				// Return the point id of the jth neighbor of the ith particle
				const unsigned int fpid = ps_fluid.neighbor(fluid_id, i, j);
				normalized_density += kernel.cubic_spline(fluid_particles.at(i), fluid_particles.at(fpid));
			}
			// Add contribution of fluid particle itself
			normalized_density += kernel.cubic_spline(fluid_particles.at(i), fluid_particles.at(i));
			// Get neighbors among grid 
			WCSPH::Vector pos = fluid_particles.at(i);
			WCSPH::Vector bottom_left = { pos[0] - parameters.compact_support, pos[1] - parameters.compact_support, pos[2] - parameters.compact_support };
			WCSPH::Vector top_right = { pos[0] + parameters.compact_support, pos[1] + parameters.compact_support, pos[2] + parameters.compact_support };
			WCSPH::Vector bottom_left_grid = this->grid_points[0];
			WCSPH::Vector starting_counter = { std::floor((bottom_left[0] - bottom_left_grid[0]) / parameters.grid_cell_size), std::floor((bottom_left[1] - bottom_left_grid[1]) / parameters.grid_cell_size), std::floor((bottom_left[2] - bottom_left_grid[2]) / parameters.grid_cell_size)};
			WCSPH::Vector ending_counter = { std::ceil((top_right[0] - bottom_left_grid[0]) / parameters.grid_cell_size), std::ceil((top_right[1] - bottom_left_grid[1]) / parameters.grid_cell_size), std::ceil((top_right[2] - bottom_left_grid[2]) / parameters.grid_cell_size) };
			for (int x = starting_counter[0]; x <= ending_counter[0]; x++) {
				for (int y = starting_counter[1]; y <= ending_counter[1]; y++) {
					for (int z = starting_counter[2]; z <= ending_counter[2]; z++) {
						WCSPH::Vector grid_point = { bottom_left_grid[0] + x * parameters.grid_cell_size, bottom_left_grid[1] + y * parameters.grid_cell_size, bottom_left_grid[2] + z * parameters.grid_cell_size };
						CompactNSearch::Real distance = (grid_point - pos).norm();
						if (distance <= parameters.compact_support) {
							CompactNSearch::Real value = kernel.cubic_spline(pos, grid_point) / normalized_density;
							this->grid_values[x * grid_resolution[1] * grid_resolution[2] + y * grid_resolution[2] + z] += value;
						}
					}
				}
			}
		}
		if (parameters.sparse) {
			this->grid_values_map_sparse.clear();
			for (int g = 0; g < grid_values.size(); g++) {
				if (std::abs(grid_values[g] + parameters.initial_grid_value) < std::numeric_limits<double>::epsilon()) {
					continue;
				}
				else {
					this->grid_values_map_sparse.try_emplace((uint64_t)(g), grid_values[g]);
				}
			}
		}
	}

	void SPH::calculate_mc_normals() {
		if (surface_vertices.size() > 0) {
			this->surface_normals.resize(surface_vertices.size());
			CompactNSearch::NeighborhoodSearch cns_mc{ parameters.compact_support };
			unsigned int fluid_id_mc = cns_mc.add_point_set(fluid_particles.front().data(), fluid_particles.size());
			unsigned int vertices_id_mc = cns_mc.add_point_set(surface_vertices.front().data(), surface_vertices.size());
			CompactNSearch::PointSet const& ps_fluid_mc = cns_mc.point_set(fluid_id_mc);
			CompactNSearch::PointSet const& ps_vertices_mc = cns_mc.point_set(vertices_id_mc);
			cns_mc.find_neighbors();
			#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
			for (int i = 0; i < ps_fluid_mc.n_points(); ++i) {
				// Inititalize density
				CompactNSearch::Real normalized_density = 0.0;
				// Get neighbors among fluid particles
				for (int j = 0; j < ps_fluid_mc.n_neighbors(fluid_id_mc, i); ++j) {
					// Return the point id of the jth neighbor of the ith particle
					const unsigned int fpid = ps_fluid_mc.neighbor(fluid_id_mc, i, j);
					normalized_density += kernel.cubic_spline(fluid_particles.at(i), fluid_particles.at(fpid));
				}
				// Add contribution of fluid particle itself
				normalized_density += kernel.cubic_spline(fluid_particles.at(i), fluid_particles.at(i));
				// Get neighbors among mc vertices
				for (int j = 0; j < ps_fluid_mc.n_neighbors(vertices_id_mc, i); ++j) {
					// Return the point id of the jth neighbor of the ith particle
					const unsigned int vpid = ps_fluid_mc.neighbor(vertices_id_mc, i, j);
					this->surface_normals[vpid] += (kernel.cubic_grad_spline(fluid_particles.at(i), surface_vertices.at(vpid)) / normalized_density);
				}
			}
			for (int i = 0; i < surface_normals.size(); i++) {
				surface_normals[i].normalize();
			}
		}
	}

	void SPH::export_data(unsigned int frame) {
		std::vector<double> particles_scalar_data(fluid_particles.size());

		// Scalar data will be the particle_id
		for (int i = 0; i < (int)fluid_particles.size(); i++) {
			particles_scalar_data[i] = i;
		}

		// Save output
		const std::string filename = this->result_path + std::to_string(frame) + ".vtk";
		geometry::write_particles_to_vtk(filename, fluid_particles, particles_scalar_data, fluid_velocities);
	}

	void SPH::export_data_with_surface(unsigned int frame) {
		// Save output
		const std::string filename = this->result_path + std::to_string(frame) + ".vtk";
		geometry::write_tri_mesh_to_vtk(filename, surface_vertices, surface_triangles, surface_normals);
	}

	void SPH::export_boundary() {
		if (has_boundary) {
			const std::string filename = this->result_path + "boundary/boundary.vtk";
			geometry::write_tri_mesh_to_vtk(filename, boundary_mesh_vertices_export, boundary_mesh_faces_export);
		}
	}

	void SPH::export_obstacles() {
		if (obstacle_mesh_vertices_export.size() > 0 && obstacle_mesh_faces_export.size() > 0) {
			const std::string filename = this->result_path + "obstacles/boundary.vtk";
			geometry::write_tri_mesh_to_vtk(filename, obstacle_mesh_vertices_export, obstacle_mesh_faces_export);
		}
	}

}