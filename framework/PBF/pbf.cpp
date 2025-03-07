#include "pbf.h"
#include "../kernel/kernel.h"
#include "../CompactNSearch/include/CompactNSearch/CompactNSearch.h"
#include "../geometry/io.h"
#include "../geometry/sampling.h"

#include <iostream>
#include <omp.h>
#include <chrono>


#define USE_Double

namespace PBD
{

	PBF::PBF(bool gravity_only, bool with_initial_velocity, bool with_surface_tension, std::string result_path, PBFParameters params, MCParameters mcparams, std::vector<Emitter::Emitter> emitter) : gravity_only(gravity_only), with_initial_velocity(with_initial_velocity), result_path(result_path), parameters(params), mcparameters(mcparams), emitter(emitter) {
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

	void PBF::simulate(CompactNSearch::Real t_end) {
		auto start = std::chrono::system_clock::now();

		initialize();

		// Add initial velocity if needed
		if (with_initial_velocity) {
			add_relative_velocity(this->relative_velocity);
		}

		// Compute boundary mass and add boundary particles to cns
		unsigned int boundary_particles_id = 0;
		if (has_boundary) {
			unsigned int boundary_particles_id = this->cns_boundary.add_point_set(boundary_particles.front().data(), boundary_particles.size());
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
					case PBFParameters::export_type::EXPORT:
						export_data(frame_idx);
						break;
					case PBFParameters::export_type::EXPORT_WITH_SURFACE:
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

			// Retain old particle positions
			this->old_fluid_particles = fluid_particles;

			if (fluid_particles_id == -1) {
				fluid_particles_id = this->cns.add_point_set(fluid_particles.front().data(), fluid_particles.size());
			}
			else {
				// Neighborhood Search
				this->cns.resize_point_set(fluid_particles_id, fluid_particles.front().data(), fluid_particles.size());
			}
			CompactNSearch::PointSet& pointset_fluid = cns.point_set(fluid_particles_id);
			this->cns.find_neighbors();

			if (!gravity_only) {
				// Compute fluid density
				calculate_particle_density(fluid_particles_id, pointset_fluid, boundary_particles_id);

				// Compute pressure
				//calculate_pressure(t_sim == 0.0);
			}

			// For surface tension compute fluid normals
			if (with_surface_tension) {
				calculate_fluid_normals(fluid_particles_id, pointset_fluid);
			}

			// Compute accelerations
			calculate_acceleration(fluid_particles_id, pointset_fluid, boundary_particles_id);

			// Time integration
			semi_implicit_euler();

			// Set up for PBF iterations
			cns.resize_point_set(fluid_particles_id, fluid_particles.front().data(), fluid_particles.size());
			pointset_fluid = cns.point_set(fluid_particles_id);
			cns.find_neighbors();

			for (int iter = 0; iter < parameters.pbf_iterations; iter++) {
				// do a pbf iteration
				update_particle_positions(fluid_particles_id, pointset_fluid, boundary_particles_id);
			}

			update_particle_velocities();

			update_active_status();
			
			// Check particle positions
			if (has_boundary) {
				check_particle_positions();
			}

			//Export VTK
			if (t_sim >= t_next_frame) {
				if (parameters.export_type == PBFParameters::export_type::EXPORT_WITH_SURFACE) {
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
					calculate_mc_normals();
				}

				switch(parameters.export_type) {
					case PBFParameters::export_type::EXPORT : 
						export_data(frame_idx);
						break;
					case PBFParameters::export_type::EXPORT_WITH_SURFACE:
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

	void PBF::initialize() {
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
		this->fluid_normals.resize(fluid_particles.size());
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < (int)this->fluid_normals.size(); i++) {
			this->fluid_normals[i] = { 0.0,0.0,0.0 };
		}
	}

	void PBF::load_geometry(bool has_boundary, bool open_boundary, PBD::Vector& boundary_size, PBD::Vector& bottom_left_boundary, std::vector<PBD::Vector>& fluid_sizes, std::vector<PBD::Vector>& bottom_lefts_fluid, std::pair <PBD::Vector, PBD::Vector>& obstacle_box, CompactNSearch::Real obstacle_sphere_radius) {
		auto start = std::chrono::system_clock::now();
		
		this->has_boundary = has_boundary;
		this->open_boundary = open_boundary;
		if (has_boundary) {
			this->boundary_box_bottom = bottom_left_boundary;
			this->boundary_box_top = bottom_left_boundary + boundary_size;
			std::vector<PBD::Vector> vertices = { bottom_left_boundary, bottom_left_boundary + PBD::Vector(0.0, 0.0, boundary_size[2]), bottom_left_boundary + PBD::Vector(0.0, boundary_size[1], 0.0),
				bottom_left_boundary + PBD::Vector(0.0, boundary_size[1], boundary_size[2]), bottom_left_boundary + PBD::Vector(boundary_size[0], 0.0, 0.0), bottom_left_boundary + PBD::Vector(boundary_size[0], 0.0, boundary_size[2]),
				bottom_left_boundary + PBD::Vector(boundary_size[0], boundary_size[1], 0.0), bottom_left_boundary + boundary_size };

			std::vector<std::array<int, 3>> triangles = { {1, 2, 0}, {3, 6, 2}, {7, 4, 6}, {5, 0, 4}, {6, 0, 2}, {3, 5, 7},
														{1, 3, 2}, {3, 7, 6}, {7, 5, 4}, {5, 1, 0}, {6, 4, 0}, {3, 1, 5} };
			if (open_boundary) {
				triangles = { {1, 2, 0}, {3, 6, 2}, {7, 4, 6}, {5, 0, 4}, {6, 0, 2},
							{1, 3, 2}, {3, 7, 6}, {7, 5, 4}, {5, 1, 0}, {6, 4, 0} };
			}

			int vertices_before = vertices.size();

			boundary_mesh_faces_export = triangles;
			boundary_mesh_vertices_export = vertices;

			if (emitter.size() > 0) {
				for (unsigned int i = 0; i < emitter.size(); i++) {
					std::vector<PBD::Vector> shield_vertices = emitter[i].get_shield_vertices();
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

			if (obstacle_box.second[0] > 0.001) {
				std::cout << "Obstacle box found:" << obstacle_box.second[0] << std::endl;
				this->obstacle_data = obstacle_box;
				this->has_obstacle = true;
				std::vector<PBD::Vector> obstacle_vertices;
				std::vector<std::array<int, 3>> obstacle_triangles;
				PBD::Vector obstacle_box_bottom = obstacle_box.first;
				PBD::Vector obstacle_size = obstacle_box.second;
				obstacle_vertices = { obstacle_box_bottom, obstacle_box_bottom + PBD::Vector(0.0, 0.0, obstacle_size[2]), obstacle_box_bottom + PBD::Vector(0.0, obstacle_size[1], 0.0),
					obstacle_box_bottom + PBD::Vector(0.0, obstacle_size[1], obstacle_size[2]), obstacle_box_bottom + PBD::Vector(obstacle_size[0], 0.0, 0.0), obstacle_box_bottom + PBD::Vector(obstacle_size[0], 0.0, obstacle_size[2]),
					obstacle_box_bottom + PBD::Vector(obstacle_size[0], obstacle_size[1], 0.0), obstacle_box_bottom + obstacle_size };

				obstacle_triangles = { {1, 2, 0}, {3, 6, 2}, {7, 4, 6}, {5, 0, 4}, {6, 0, 2}, {3, 5, 7},
									 {1, 3, 2}, {3, 7, 6}, {7, 5, 4}, {5, 1, 0}, {6, 4, 0}, {3, 1, 5} };

				for (int i = 0; i < obstacle_vertices.size(); i++) {
					vertices.push_back(obstacle_vertices[i]);
				}

				for (int j = 0; j < obstacle_triangles.size(); j++) {
					std::array<int, 3> triangle = obstacle_triangles[j];
					triangle[0] += vertices_before;
					triangle[1] += vertices_before;
					triangle[2] += vertices_before;
					triangles.push_back(triangle);
				}

				obstacle_mesh_vertices_export = obstacle_vertices;
				obstacle_mesh_faces_export = obstacle_triangles;
			}

			vertices_before = vertices.size();

			if (obstacle_sphere_radius > 0.0) {
				std::cout << "Obstacle sphere found!" << std::endl;
				this->has_obstacle_sphere = true;
				this->obstacle_sphere_radius = obstacle_sphere_radius;
				const std::vector<geometry::TriMesh> meshes = geometry::read_tri_meshes_from_obj("../res/sphere.obj");
				const geometry::TriMesh& sphere = meshes[0];
				std::vector<geometry::Vector> sphere_vertices = sphere.vertices;
				const std::vector<std::array<int, 3>> sphere_triangles = sphere.triangles;

				CompactNSearch::Real scale = 2.0 / obstacle_sphere_radius;
				for (int i = 0; i < sphere_vertices.size(); i++) {
					sphere_vertices[i] = sphere_vertices[i] / scale;
					vertices.push_back(sphere_vertices[i]);
				}

				for (int j = 0; j < sphere_triangles.size(); j++) {
					std::array<int, 3> triangle = sphere_triangles[j];
					triangle[0] += vertices_before;
					triangle[1] += vertices_before;
					triangle[2] += vertices_before;
					triangles.push_back(triangle);
				}

				obstacle_sphere_mesh_vertices_export = sphere_vertices;
				obstacle_sphere_mesh_faces_export = sphere_triangles;
			}

			if (false) {
				const std::vector<geometry::TriMesh> meshes = geometry::read_tri_meshes_from_obj("../res/meshes/riverbed.obj");
				for (int i = 0; i < (int)meshes.size(); i++) {
					std::cout << "mesh " << i << std::endl;
					for (int j = 0; j < (int)meshes[i].vertices.size(); j++) {
						vertices.push_back(meshes[i].vertices[j]);
					}
					for (int j = 0; j < (int)meshes[i].triangles.size(); j++) {
						if (meshes[i].triangles[j][0] > -1 && meshes[i].triangles[j][1] > -1 && meshes[i].triangles[j][2] > -1) {
							std::array<int, 3> new_triangle = meshes[i].triangles[j];
							for (int k = 0; k < 3; k++) {
								new_triangle[k] += vertices_before;
							}
							triangles.push_back(new_triangle);
						}
					}
					vertices_before = vertices.size();
				}
			}

			this->boundary_mesh_vertices = vertices;
			this->boundary_mesh_faces = triangles;
			geometry::sampling::triangle_mesh(this->boundary_particles, this->boundary_mesh_vertices, this->boundary_mesh_faces, this->parameters.boundary_sampling_distance);
		}

		//this->particle_mass = 0.0;
		CompactNSearch::Real total_fluid_volume = 0.0;
		for (int i = 0; i < fluid_sizes.size(); ++i) {
			PBD::Vector top_right_fluid = bottom_lefts_fluid[i] + fluid_sizes[i];
			if (has_boundary) {
				PBD::Vector top_right_boundary = bottom_left_boundary + boundary_size;
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

	void PBF::emit_particles(CompactNSearch::Real t_sim) {
		//Emit particles
		for (unsigned int i = 0; i < emitter.size(); i++) {
			std::vector<PBD::Vector> emitted_particle_positions;
			std::vector<PBD::Vector> emitted_particle_velocities;

			emitter[i].emit_particles(emitted_particle_positions, emitted_particle_velocities, t_sim);
			if (emitted_particle_positions.size() == 0) {
				continue;
			}
			bool emit = true;
			//#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
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
		this->fluid_normals.resize(fluid_particles.size());
	}

	void PBF::update_active_status() {
		//update active status
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < fluid_particles.size(); i++) {
			if (!is_fluid_particle_active[i]) {
				PBD::Vector emission_position;
				PBD::Vector dir = fluid_velocities[i].normalized();
				//Find origin position
				for (unsigned int j = 0; j < emitter.size(); j++) {
					std::vector<PBD::Vector> emission_positions = emitter[j].get_emition_positions();
					bool found = false;
					for (unsigned int k = 0; k < emission_positions.size(); k++) {
						PBD::Vector dist = (fluid_particles[i] - emission_positions[k]).normalized();
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

	void PBF::compute_boundary_mass(unsigned int boundary_id, CompactNSearch::PointSet const& ps_boundary) {
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

	void PBF::calculate_particle_density(unsigned int fluid_id, CompactNSearch::PointSet const& ps_fluid, unsigned int boundary_id) {
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

	void PBF::semi_implicit_euler() {
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < (int)fluid_particles.size(); i++) {
			PBD::Vector new_v = fluid_velocities[i] + parameters.dt * fluid_accelerations[i];
			PBD::Vector new_x = fluid_particles[i] + parameters.dt * new_v;

			//fluid_velocities[i] = new_v; // Do we need this?
			fluid_particles[i] = new_x;
		}
	}

	void PBF::calculate_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id) {
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < (int)fluid_accelerations.size(); i++) {
			fluid_accelerations[i] = { 0.0,0.0,0.0 };
		}

		if (!gravity_only) {
			calculate_viscosity_acceleration(fluid_id, ps_fluid, boundary_id);
		}

		calculate_other_acceleration();
		if (with_surface_tension) {
			calculate_st_acceleration(fluid_id, ps_fluid, boundary_id);
			if (has_boundary) {
				calculate_adhesion_acceleration(fluid_id, ps_fluid, boundary_id);
			}
		}
	}

	void PBF::calculate_viscosity_acceleration(unsigned int fluid_id, CompactNSearch::PointSet const& ps_fluid, unsigned int boundary_id) {
		// Loop over all fluid particles
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < ps_fluid.n_points(); ++i) {
			if (is_fluid_particle_active[i]) {
				// Initialize acceleration
				PBD::Vector acceleration = { 0.0, 0.0, 0.0 };
				// Get fluid neighbors of particles
				for (int j = 0; j < ps_fluid.n_neighbors(fluid_id, i); ++j) {
					// Return the point id of the jth neighbor of the ith particle
					const unsigned int fpid = ps_fluid.neighbor(fluid_id, i, j);
					// Compute fluid-fluid interaction
					const CompactNSearch::Real normff = (fluid_particles.at(i) - fluid_particles.at(fpid)).norm();
					PBD::Vector ffa = kernel.cubic_grad_spline(fluid_particles.at(i), fluid_particles.at(fpid)).dot(fluid_particles.at(i) - fluid_particles.at(fpid)) * (fluid_velocities.at(i) - fluid_velocities.at(fpid)) / (fluid_densities.at(fpid) * ((normff * normff) + (0.01 * parameters.smoothing_length_squared)));
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
						PBD::Vector fba = 2 * parameters.boundary_viscosity * boundary_volumes.at(bpid) * kernel.cubic_grad_spline(fluid_particles.at(i), boundary_particles.at(bpid)).dot(fluid_particles.at(i) - boundary_particles.at(bpid)) * fluid_velocities.at(i) / ((normfb * normfb) + (0.01 * parameters.smoothing_length_squared));
						acceleration += fba;
					}
				}
				fluid_accelerations[i] += acceleration;
			}
		}
	}

	void PBF::calculate_other_acceleration() {
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

	void PBF::calculate_fluid_normals(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid) {
		// Loop over all fluid particles
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < ps_fluid.n_points(); ++i) {
			// Initialize acceleration
			PBD::Vector normal = { 0.0, 0.0, 0.0 };
			// Get fluid neighbors of particles
			for (int j = 0; j < ps_fluid.n_neighbors(fluid_id, i); ++j) {
				// Return the point id of the jth neighbor of the ith particle
				const unsigned int fpid = ps_fluid.neighbor(fluid_id, i, j);
				// Compute fluid-fluid interaction
				PBD::Vector ffa = kernel.cubic_grad_spline(fluid_particles.at(i), fluid_particles.at(fpid)) / fluid_densities.at(fpid);
				normal += ffa;
			}
			normal *= (this->particle_mass * parameters.compact_support);
			fluid_normals[i] = normal;
		}
	}

	void PBF::calculate_st_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id) {
		// Loop over all fluid particles
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < ps_fluid.n_points(); ++i) {
			if (is_fluid_particle_active[i]) {
				// Initialize acceleration
				PBD::Vector acceleration = { 0.0, 0.0, 0.0 };
				// Get fluid neighbors of particles
				for (int j = 0; j < ps_fluid.n_neighbors(fluid_id, i); ++j) {
					// Return the point id of the jth neighbor of the ith particle
					const unsigned int fpid = ps_fluid.neighbor(fluid_id, i, j);
					// Compute fluid-fluid interaction
					const CompactNSearch::Real normff = (fluid_particles.at(i) - fluid_particles.at(fpid)).norm();
					// Compute cohesion related acceleration
					PBD::Vector co_acc = this->particle_mass * kernel.cohesion_kernel(fluid_particles.at(i), fluid_particles.at(fpid)) * (fluid_particles.at(i) - fluid_particles.at(fpid)) / normff;
					// Compute curvature related acceleration
					PBD::Vector curv_acc = this->fluid_normals.at(i) - this->fluid_normals.at(fpid);
					// Add both contributions to get surface tension acceleration
					acceleration += ((co_acc + curv_acc) / (fluid_densities.at(i) + fluid_densities.at(fpid)));
				}
				acceleration *= (-2.0 * parameters.fluid_rest_density * parameters.cohesion_coefficient);
				fluid_accelerations[i] += acceleration;
			}
		}
	}

	void PBF::calculate_adhesion_acceleration(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid, unsigned int boundary_id) {
		// Loop over all fluid particles
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < ps_fluid.n_points(); ++i) {
			if (is_fluid_particle_active[i]) {
				// Initialize acceleration
				PBD::Vector acceleration = { 0.0, 0.0, 0.0 };
				// Get boundary neighbors of particles
				if (ps_fluid.n_neighbors(boundary_id, i) > 0 && has_boundary) {
					for (int j = 0; j < ps_fluid.n_neighbors(boundary_id, i); ++j) {
						// Return the point id of the jth neighbor of the ith particle
						const unsigned int bpid = ps_fluid.neighbor(boundary_id, i, j);
						// Compute fluid-boundary interaction
						const CompactNSearch::Real normfb = (fluid_particles.at(i) - boundary_particles.at(bpid)).norm();
						acceleration += boundary_masses.at(bpid) * kernel.adhesion_kernel(fluid_particles.at(i), boundary_particles.at(bpid)) * (fluid_particles.at(i) - boundary_particles.at(bpid)) / normfb;
					}
				}
				acceleration *= -1.0 * parameters.adhesion_coefficient;
				fluid_accelerations[i] += acceleration;
			}
		}
	}

	void PBF::update_time_step_size() {
		// TOOOOOOOOO DOOOOOOOOOOOOOOO: use squared norm in stead of norm to check max velocity
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

	void PBF::update_particle_positions(unsigned int fluid_id, CompactNSearch::PointSet const& ps_fluid, unsigned int boundary_id) {
		// Compute fluid particles density
		calculate_particle_density(fluid_id, ps_fluid, boundary_id);
		// Inititalize lambda values
		std::vector<CompactNSearch::Real> lambdas(fluid_particles.size(), 0.0);
		// Store new positions
		std::vector<PBD::Vector> new_positions(fluid_particles.size());
		// Keep old positions // Not here
		//this->old_fluid_particles = fluid_particles;
		
		// Compute lambdas...............................................................................................................
		// Loop over all fluid particles
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < ps_fluid.n_points(); ++i) {
			// Calculate S(i)
			// Inititalize first term of S(i) as vector
			PBD::Vector first_vector = { 0.0, 0.0, 0.0 };
			// Inititalize second term of S(i)
			CompactNSearch::Real second_term = 0.0;
			// Get fluid neighbors of particles
			for (int j = 0; j < ps_fluid.n_neighbors(fluid_id, i); ++j) {
				// Return the point id of the jth neighbor of the ith particle
				const unsigned int fpid = ps_fluid.neighbor(fluid_id, i, j);
				PBD::Vector temp_vector = this->particle_mass * kernel.cubic_grad_spline(fluid_particles.at(i), fluid_particles.at(fpid)) / parameters.fluid_rest_density;
				CompactNSearch::Real norm_second = (-1.0 * temp_vector).norm();
				// Add vector to first term of S(i)
				first_vector += (temp_vector);
				// Add to second term of S(i)
				second_term += (norm_second * norm_second);
			}
			// Get boundary neighbors of particles
			if (ps_fluid.n_neighbors(boundary_id, i) > 0 && has_boundary) {
				for (int j = 0; j < ps_fluid.n_neighbors(boundary_id, i); ++j) {
					// Return the point id of the jth neighbor of the ith particle
					const unsigned int bpid = ps_fluid.neighbor(boundary_id, i, j);
					// Add vector to first term of S(i)
					first_vector += (this->boundary_volumes.at(bpid) * kernel.cubic_grad_spline(fluid_particles.at(i), boundary_particles.at(bpid)));
				}
			}
			CompactNSearch::Real norm_first = first_vector.norm();
			CompactNSearch::Real S_i = (norm_first * norm_first) + second_term;
			S_i /= this->particle_mass;
			// Calculate lambda(i)
			CompactNSearch::Real lambda_i = 0;
			CompactNSearch::Real eps = 0.0001;
			if (fluid_densities.at(i) > parameters.fluid_rest_density) {
				CompactNSearch::Real C_i = (fluid_densities.at(i) / parameters.fluid_rest_density) - 1.0;
				lambda_i = -1.0 * C_i / (S_i + eps);
			}
			lambdas.at(i) = lambda_i;
		}

		// Compute new positions...............................................................................................................
		// Loop over all fluid particles
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < ps_fluid.n_points(); ++i) {
			if (is_fluid_particle_active[i]) {
				// Inititalize position change vector
				PBD::Vector deltax_i = { 0.0, 0.0, 0.0 };
				// Get fluid neighbors of particles
				for (int j = 0; j < ps_fluid.n_neighbors(fluid_id, i); ++j) {
					// Return the point id of the jth neighbor of the ith particle
					const unsigned int fpid = ps_fluid.neighbor(fluid_id, i, j);
					deltax_i += ((lambdas.at(i) + lambdas.at(fpid)) * kernel.cubic_grad_spline(fluid_particles.at(i), fluid_particles.at(fpid)));
				}
				// Get boundary neighbors of particles
				if (ps_fluid.n_neighbors(boundary_id, i) > 0 && has_boundary) {
					for (int j = 0; j < ps_fluid.n_neighbors(boundary_id, i); ++j) {
						// Return the point id of the jth neighbor of the ith particle
						const unsigned int bpid = ps_fluid.neighbor(boundary_id, i, j);
						deltax_i += (boundary_masses.at(bpid) * lambdas.at(i) * kernel.cubic_grad_spline(fluid_particles.at(i), boundary_particles.at(bpid)) / this->particle_mass);
					}
				}
				deltax_i /= parameters.fluid_rest_density;
				new_positions.at(i) = fluid_particles.at(i) + deltax_i;
			}
			else {
				new_positions.at(i) = fluid_particles.at(i);
			}
		}

		// Update positions
		fluid_particles = new_positions;
	}

	void PBF::update_particle_velocities() {
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < (int)fluid_velocities.size(); i++) {
			fluid_velocities[i] = (fluid_particles[i] - old_fluid_particles[i]) / parameters.dt;
		}
	}

	void PBF::check_particle_positions() {
		int old_size = (int)fluid_particles.size();
		std::vector<bool> new_is_fluid_particle_active;
		std::vector<PBD::Vector> new_positions;
		std::vector<PBD::Vector> new_velocities;
		std::vector<PBD::Vector> new_accelerations;
		std::vector<PBD::Vector> new_normals;
		std::vector<CompactNSearch::Real> new_densities;
		std::vector<CompactNSearch::Real> new_pressures;
		new_is_fluid_particle_active.reserve(old_size);
		new_positions.reserve(old_size);
		new_velocities.reserve(old_size);
		new_accelerations.reserve(old_size);
		new_normals.reserve(old_size);
		new_densities.reserve(old_size);
		new_pressures.reserve(old_size);

		for (int i = 0; i < old_size; i++) {
			PBD::Vector pos = fluid_particles[i];
			bool is_bigger_BL = pos[0] >= this->boundary_box_bottom[0] && pos[1] >= this->boundary_box_bottom[1] && pos[2] >= this->boundary_box_bottom[2];
			bool is_smaller_TR = pos[0] <= this->boundary_box_top[0] && pos[1] <= this->boundary_box_top[1] && pos[2] <= this->boundary_box_top[2];
			bool keep = is_bigger_BL && is_smaller_TR;

			if (open_boundary) {
				bool is_out_top = pos[0] <= this->boundary_box_top[0] && pos[1] <= this->boundary_box_top[1] && pos[2] > this->boundary_box_top[2] && pos[2] <= (this->boundary_box_top[2] + 4 * (this->boundary_box_top[2] - this->boundary_box_bottom[2]));
				keep = keep || is_out_top;
			}

			if (has_obstacle) {
				PBD::Vector obstacle_bottom = this->obstacle_data.first;
				PBD::Vector obstacle_top = this->obstacle_data.first + this->obstacle_data.second;
				bool is_bigger_BL_obstacle = pos[0] >= obstacle_bottom[0] && pos[1] >= obstacle_bottom[1] && pos[2] >= obstacle_bottom[2];
				bool is_smaller_TR_obstacle = pos[0] <= obstacle_top[0] && pos[1] <= obstacle_top[1] && pos[2] <= obstacle_top[2];
				bool is_inside_obstacle = is_bigger_BL_obstacle && is_smaller_TR_obstacle;
				keep = keep && !is_inside_obstacle;
			}

			if (has_obstacle_sphere) {
				PBD::Vector obstacle_center = PBD::Vector(0.0, 0.0, 0.0);
				CompactNSearch::Real obstacle_radius = this->obstacle_sphere_radius;
				bool is_inside_obstacle_sphere = (pos - obstacle_center).norm() <= obstacle_radius;
				keep = keep && !is_inside_obstacle_sphere;
			}

			if (!keep) {
				continue;
			}

			new_is_fluid_particle_active.emplace_back(is_fluid_particle_active[i]);
			new_positions.emplace_back(pos);
			new_velocities.emplace_back(fluid_velocities[i]);
			new_accelerations.emplace_back(fluid_accelerations[i]);
			new_normals.emplace_back(fluid_normals[i]);
			new_densities.emplace_back(fluid_densities[i]);
			new_pressures.emplace_back(fluid_pressures[i]);
		}
		this->is_fluid_particle_active = new_is_fluid_particle_active;
		this->fluid_particles = new_positions;
		this->fluid_velocities = new_velocities;
		this->fluid_accelerations = new_accelerations;
		this->fluid_normals = new_normals;
		this->fluid_densities = new_densities;
		this->fluid_pressures = new_pressures;
	}
	
	void PBF::turn_off_gravity() {
		this->gravity_acceleration = { 0.0, 0.0, 0.0 };
	}

	void PBF::turn_on_gravity() {
		this->gravity_acceleration = { 0.0, 0.0, -9.81 };
	}

	void PBF::add_initial_velocity(std::vector<PBD::Vector> velocities) {
		assert(velocities.size() == this->fluid_velocities.size());
		this->fluid_velocities = velocities;
	}

	void PBF::add_relative_velocity(PBD::Vector relative_velocity) {
		std::vector<PBD::Vector> velocities;
		for (int i = 0; i < (int)this->fluid_particle_counts[0]; i++) {
			velocities.emplace_back(0.5 * relative_velocity);
		}
		for (int i = (int)this->fluid_particle_counts[0]; i < (int)this->fluid_particle_counts[1]; i++) {
			velocities.emplace_back(-0.5 * relative_velocity);
		}
		assert(velocities.size() == this->fluid_velocities.size());
		this->fluid_velocities = velocities;
	}

	int PBF::num_boundary_particles() {
		return boundary_particles.size();
	}

	int PBF::num_fluid_particles() {
		int num_particles = 0;
		for (int i = 0;i < fluid_particle_counts.size();i++) {
			num_particles += fluid_particle_counts[i];
		}
		return num_particles;
	}

	int PBF::num_timesteps() {
		return timesteps;
	}

	void PBF::printStats(std::string test) {
		int num_boundary = boundary_particles.size();
		//int num_fluid = 0;
		//for (int i = 0;i < fluid_particle_counts.size();i++) {
			//num_fluid += fluid_particle_counts[i];
		//}
		int num_fluid = fluid_particles.size();

		std::cout << test << ": " << "(boundary: " << num_boundary << ", fluid: " << num_fluid << ", timesteps: " << timesteps << ")" << std::endl;
		std::cout << " => Time: " << duration << " ms" << std::endl;
	}

	void PBF::create_grid() {
		grid_points.clear();
		PBD::Vector bottom_left = { std::numeric_limits<CompactNSearch::Real>::max(), std::numeric_limits<CompactNSearch::Real>::max(), std::numeric_limits<CompactNSearch::Real>::max() };
		PBD::Vector top_right = { std::numeric_limits<CompactNSearch::Real>::min(), std::numeric_limits<CompactNSearch::Real>::min(), std::numeric_limits<CompactNSearch::Real>::min() };

		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < fluid_particles.size(); i++) {
			PBD::Vector pos = fluid_particles[i];
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
					PBD::Vector new_point = { bottom_left[0] + i * parameters.grid_cell_size, bottom_left[1] + j * parameters.grid_cell_size, bottom_left[2] + k * parameters.grid_cell_size };
					grid_points.emplace_back(new_point);
				}
			}
		}
	}

	void PBF::reset_grid_values() {
		this->grid_values.resize(grid_points.size());
		#pragma omp parallel for num_threads(parameters.num_threads) schedule(static)
		for (int i = 0; i < this->grid_values.size(); i++) {
			this->grid_values[i] = -1.0 * parameters.initial_grid_value;
		}
	}

	void PBF::update_grid_values(unsigned int fluid_id, CompactNSearch::PointSet& ps_fluid) {
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
			PBD::Vector pos = fluid_particles.at(i);
			PBD::Vector bottom_left = { pos[0] - parameters.compact_support, pos[1] - parameters.compact_support, pos[2] - parameters.compact_support };
			PBD::Vector top_right = { pos[0] + parameters.compact_support, pos[1] + parameters.compact_support, pos[2] + parameters.compact_support };
			PBD::Vector bottom_left_grid = this->grid_points[0];
			PBD::Vector starting_counter = { std::floor((bottom_left[0] - bottom_left_grid[0]) / parameters.grid_cell_size), std::floor((bottom_left[1] - bottom_left_grid[1]) / parameters.grid_cell_size), std::floor((bottom_left[2] - bottom_left_grid[2]) / parameters.grid_cell_size)};
			PBD::Vector ending_counter = { std::ceil((top_right[0] - bottom_left_grid[0]) / parameters.grid_cell_size), std::ceil((top_right[1] - bottom_left_grid[1]) / parameters.grid_cell_size), std::ceil((top_right[2] - bottom_left_grid[2]) / parameters.grid_cell_size) };
			for (int x = starting_counter[0]; x <= ending_counter[0]; x++) {
				for (int y = starting_counter[1]; y <= ending_counter[1]; y++) {
					for (int z = starting_counter[2]; z <= ending_counter[2]; z++) {
						PBD::Vector grid_point = { bottom_left_grid[0] + x * parameters.grid_cell_size, bottom_left_grid[1] + y * parameters.grid_cell_size, bottom_left_grid[2] + z * parameters.grid_cell_size };
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

	void PBF::calculate_mc_normals() {
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

	void PBF::export_data(unsigned int frame) {
		std::vector<double> particles_scalar_data(fluid_particles.size());

		// Scalar data will be the particle_id
		for (int i = 0; i < (int)fluid_particles.size(); i++) {
			particles_scalar_data[i] = i;
		}

		// Save output
		const std::string filename = this->result_path + std::to_string(frame) + ".vtk";
		geometry::write_particles_to_vtk(filename, fluid_particles, particles_scalar_data, fluid_velocities);
	}

	void PBF::export_data_with_surface(unsigned int frame) {
		// Save output
		const std::string filename = this->result_path + std::to_string(frame) + ".vtk";
		geometry::write_tri_mesh_to_vtk(filename, surface_vertices, surface_triangles, surface_normals);
	}

	void PBF::export_boundary() {
		const std::string filename = this->result_path + "boundary/boundary.vtk";
		geometry::write_tri_mesh_to_vtk(filename, boundary_mesh_vertices_export, boundary_mesh_faces_export);
	}

	void PBF::export_obstacles() {
		if (obstacle_mesh_vertices_export.size() > 0 && obstacle_mesh_faces_export.size() > 0) {
			const std::string filename = this->result_path + "obstacles/boundary.vtk";
			geometry::write_tri_mesh_to_vtk(filename, obstacle_mesh_vertices_export, obstacle_mesh_faces_export);
		}
	}

}