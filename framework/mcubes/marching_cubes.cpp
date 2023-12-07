#include <iostream>
#include "marching_cubes_lut.h"
#include "marching_cubes.h"
#include <unordered_set>

namespace MC
{
	MarchingCubes::MarchingCubes(MCParameters mcparams, CompactNSearch::Real value) : mcparameters(mcparams), isovalue(value) {

	}

	void MarchingCubes::calculate(std::vector<MC::Vector>& positions, std::vector<CompactNSearch::Real>& values, std::unordered_map<uint64_t, CompactNSearch::Real>& value_map, std::array<int, 3> dimensions, std::vector<MC::Vector>& vertices, std::vector<std::array<int, 3>>& triangles) {
		assert(positions.size() == values.size());
		if (dimensions[0] > 1024 || dimensions[1] > 1024 || dimensions[2] > 1024) {
			std::cout << "WARNING GRID DIMENSION MIGHT CAUSE OVERFLOW" << std::endl;
		}
		if (mcparameters.sparse) {
			calculate_sparse(positions, value_map, dimensions, vertices, triangles);
		}
		else {
			if (mcparameters.ours) {
				calculate_ours(positions, values, dimensions, vertices, triangles);
			}
			else {
				calculate_regular(positions, values, dimensions, vertices, triangles);
			}
		}
	}

	MC::Vector MarchingCubes::calculate_interpolation(std::vector<MC::Vector>& positions, std::vector<CompactNSearch::Real>& values, int starting_index, int ending_index) {
		MC::Vector p;
		CompactNSearch::Real mu;
		CompactNSearch::Real v1 = values[starting_index];
		MC::Vector p1 = positions[starting_index];

		CompactNSearch::Real v2 = values[ending_index];
		MC::Vector p2 = positions[ending_index];

		if (std::abs(isovalue - v1) < 0.00001) {
			p = p1;
		}
		else if (std::abs(isovalue - v2) < 0.00001) {
			p = p2;
		}
		else if (std::abs(v1 - v2) < 0.00001) {
			p = p1;
		}
		else {
			mu = (isovalue - v1) / (v2 - v1);
			p = p1 + mu * (p2 - p1);
		}
		return p;
	}

	MC::Vector MarchingCubes::calculate_interpolation_sparse(std::vector<MC::Vector>& positions, CompactNSearch::Real v1, CompactNSearch::Real v2, int starting_index, int ending_index) {
		MC::Vector p;
		CompactNSearch::Real mu;
		MC::Vector p1 = positions[starting_index];
		MC::Vector p2 = positions[ending_index];

		if (std::abs(isovalue - v1) < 0.00001) {
			p = p1;
		}
		else if (std::abs(isovalue - v2) < 0.00001) {
			p = p2;
		}
		else if (std::abs(v1 - v2) < 0.00001) {
			p = p1;
		}
		else {
			mu = (isovalue - v1) / (v2 - v1);
			p = p1 + mu * (p2 - p1);
		}
		return p;
	}

	void MarchingCubes::calculate_regular(std::vector<MC::Vector>& positions, std::vector<CompactNSearch::Real>& values, std::array<int, 3> dimensions, std::vector<MC::Vector>& vertices, std::vector<std::array<int, 3>>& triangles){
		std::vector<MC::Vector> intersection_points;
		std::unordered_map<uint64_t, uint64_t> edge_map;
		// loop through the left and inner vertices of the grid
		// then loop through edges connected to the vertex
		// if the edge is not in the map, calculate the intersection point and add it to the map
		CompactNSearch::Real epsilon = 0.000001;
		std::vector<MC::Vector>::iterator it;
		for (int x = 0; x < dimensions[0]; x++) {
			for (int y = 0; y < dimensions[1]; y++) {
				for (int z = 0; z < dimensions[2]; z++) {
					int vertex_idx = (dimensions[2] * dimensions[1] * x) + (dimensions[2] * y) + z;
					CompactNSearch::Real current_value = values[vertex_idx];
					int dir_idx = 0;
					for (int dir = 0; dir < 3; dir++) {
						if (dir == 0) {
							dir_idx = vertex_idx + dimensions[2] * dimensions[1];
						}
						else if (dir == 1) {
							dir_idx = vertex_idx + dimensions[2];
						}
						else {
							dir_idx = vertex_idx + 1;
						}
						if (dir_idx < positions.size()) {
							CompactNSearch::Real dir_value = values[dir_idx];
							if (current_value * dir_value <= isovalue) {
								MC::Vector p = calculate_interpolation(positions, values, vertex_idx, dir_idx);
								uint64_t edge_key = (uint64_t)(3 * vertex_idx + dir);
								if (edge_map.find(edge_key) == edge_map.end()) {
									edge_map[edge_key] = intersection_points.size();
									intersection_points.push_back(p);
								}
							}
						}
					}
				}
			}
		}
		#pragma omp parallel for num_threads(mcparameters.num_threads) schedule(static)
		for (int x = 0; x < dimensions[0] - 1; x++) {
			#pragma omp parallel for num_threads(mcparameters.num_threads) schedule(static)
			for (int y = 0; y < dimensions[1] - 1; y++) {
				#pragma omp parallel for num_threads(mcparameters.num_threads) schedule(static)
				for (int z = 0; z < dimensions[2] - 1; z++) {
					int cube_mask = 0;
					//std::array<bool, 8> vertex_signs = { 0,0,0,0,0,0,0,0 };
					std::array<int, 8> cell_vertex_indices = { 0,0,0,0,0,0,0,0 };
					//Calculate cube mask
					for (unsigned int i = 0; i < CELL_VERTICES.size(); i++) {
						std::array<int, 3> cell_vertices = CELL_VERTICES[i];
						int idx = (dimensions[2] * dimensions[1] * (x + cell_vertices[0])) + (dimensions[2] * (y + cell_vertices[1])) + z + cell_vertices[2];
						cube_mask |= values[idx] < isovalue ? (1 << i) : 0;
						//vertex_signs[i] = values[idx] >= isovalue;
						cell_vertex_indices[i] = idx;
					}

					if (cube_mask == 0 || cube_mask == 255) {
						continue;
					}
					//Add triangle indices
					//std::array<std::array<int, 3>, 5> mc_table_data = get_marching_cubes_cell_triangulation(vertex_signs);
					for (unsigned int i = 0; i < 5; i++) {
						std::array<int, 3> tri = MARCHING_CUBES_TABLE[cube_mask][i];
						if (tri[0] != -1) {
							#pragma omp critical
							{
								std::array<int, 3> triangle = { 0, 0, 0 };
								for (int j = 0; j < 3; j++) {
									int edge_idx = tri[j];
									std::array<int, 2> edge = CELL_EDGES[edge_idx];
									int idx_1 = cell_vertex_indices[edge[0]];
									int idx_2 = cell_vertex_indices[edge[1]];
									int idx_small = std::min(idx_1, idx_2);
									uint64_t edge_key = (uint64_t)(3 * idx_small + CELL_EDGES_DIRECTION[edge_idx]);
									int vertex_idx = edge_map[edge_key];
									triangle[j] = vertex_idx;
								}
								triangles.push_back(triangle);
							}
						}
					}

				}
			}
		}

		vertices.resize(intersection_points.size());
		vertices = intersection_points;
	}

	void MarchingCubes::calculate_sparse(std::vector<MC::Vector>& positions, std::unordered_map<uint64_t, CompactNSearch::Real>& value_map, std::array<int, 3> dimensions, std::vector<MC::Vector>& vertices, std::vector<std::array<int, 3>>& triangles){
		std::vector<MC::Vector> intersection_points;
		std::unordered_map<uint64_t, uint64_t> edge_vertex_map;
		std::unordered_set<uint64_t> intersected_cells;
		//std::unordered_map<uint64_t, std::vector<uint64_t>> cell_edge_map;
		//std::unordered_map<uint64_t, CompactNSearch::Real>::iterator v;
		CompactNSearch::Real default_value = -1.0 * mcparameters.initial_grid_value;
		CompactNSearch::Real dir_value = default_value;
		CompactNSearch::Real current_value = default_value;
		
		for (std::unordered_map<uint64_t, CompactNSearch::Real>::iterator v = value_map.begin(); v != value_map.end(); v++) {
			int grid_point_idx = v->first;
			current_value =v->second;
			int x = grid_point_idx / (dimensions[2] * dimensions[1]);
			int y = (grid_point_idx - (x * dimensions[2] * dimensions[1])) / dimensions[2];
			int z = grid_point_idx - (x * dimensions[2] * dimensions[1]) - (y * dimensions[2]);
			int cell_idx = ((dimensions[2] - 1) * (dimensions[1] - 1) * x) + ((dimensions[2] - 1) * y) + z;
			
			if (current_value <= isovalue) {
				// positive z direction
				int dir_idx = grid_point_idx + 1;
				if (dir_idx < positions.size()) {
					if (value_map.find(dir_idx) != value_map.end()) {
						dir_value = value_map[dir_idx];
					}
					else {
						dir_value = default_value;
					}
					if (current_value * dir_value <= isovalue) {
						MC::Vector p = calculate_interpolation_sparse(positions, current_value, dir_value, grid_point_idx, dir_idx);
						uint64_t edge_key = (uint64_t)(3 * grid_point_idx + 2);
						if (edge_vertex_map.find(edge_key) == edge_vertex_map.end()) {
							edge_vertex_map[edge_key] = intersection_points.size();
							intersection_points.push_back(p);
						}
						intersected_cells.insert(cell_idx);
						if (cell_idx - (dimensions[2] - 1) >= 0) {
							intersected_cells.insert(cell_idx - (dimensions[2] - 1));
						}
						if ((cell_idx - ((dimensions[2] - 1) * (dimensions[1] - 1))) >= 0) {
							intersected_cells.insert(cell_idx - ((dimensions[2] - 1) * (dimensions[1] - 1)));
						}
						if ((cell_idx - ((dimensions[2] - 1) * (dimensions[1] - 1)) - (dimensions[2] - 1)) >= 0) {
							intersected_cells.insert(cell_idx - ((dimensions[2] - 1) * (dimensions[1] - 1)) - (dimensions[2] - 1));
						}
						//cell_edge_map[cell_idx].push_back(edge_key);
					}
				}
				// negative z direction
				dir_idx = grid_point_idx - 1;
				if (dir_idx >= 0) {
					if (value_map.find(dir_idx) != value_map.end()) {
						dir_value = value_map[dir_idx];
					}
					else {
						dir_value = default_value;
					}
					if (current_value * dir_value <= isovalue) {
						MC::Vector p = calculate_interpolation_sparse(positions, dir_value, current_value, dir_idx, grid_point_idx);
						uint64_t edge_key = (uint64_t)(3 * dir_idx + 2);
						if (edge_vertex_map.find(edge_key) == edge_vertex_map.end()) {
							edge_vertex_map[edge_key] = intersection_points.size();
							intersection_points.push_back(p);
						}
						intersected_cells.insert(cell_idx - 1);
						if (cell_idx - 1 - (dimensions[1] - 1) >= 0) {
							intersected_cells.insert(cell_idx - 1 - (dimensions[1] - 1));
						}
						if ((cell_idx - 1 - ((dimensions[2] - 1) * (dimensions[1] - 1))) >= 0) {
							intersected_cells.insert(cell_idx - 1 - ((dimensions[2] - 1) * (dimensions[1] - 1)));
						}
						if ((cell_idx - 1 - ((dimensions[2] - 1) * (dimensions[1] - 1)) - (dimensions[1] - 1)) >= 0) {
							intersected_cells.insert(cell_idx - 1 - ((dimensions[2] - 1) * (dimensions[1] - 1)) - (dimensions[1] - 1));
						}
						//cell_edge_map[cell_idx - 1].push_back(edge_key);
					}
				}
				// positive y direction
				dir_idx = grid_point_idx + dimensions[2];
				if (dir_idx < positions.size()) {
					if (value_map.find(dir_idx) != value_map.end()) {
						dir_value = value_map[dir_idx];
					}
					else {
						dir_value = default_value;
					}
					if (current_value * dir_value <= isovalue) {
						MC::Vector p = calculate_interpolation_sparse(positions, current_value, dir_value, grid_point_idx, dir_idx);
						uint64_t edge_key = (uint64_t)(3 * grid_point_idx + 1);
						if (edge_vertex_map.find(edge_key) == edge_vertex_map.end()) {
							edge_vertex_map[edge_key] = intersection_points.size();
							intersection_points.push_back(p);
						}
						intersected_cells.insert(cell_idx);
						if (cell_idx - 1 >= 0) {
							intersected_cells.insert(cell_idx - 1);
						}
						if ((cell_idx - ((dimensions[2] - 1) * (dimensions[1] - 1))) >= 0) {
							intersected_cells.insert(cell_idx - ((dimensions[2] - 1) * (dimensions[1] - 1)));
						}
						if ((cell_idx - 1 - ((dimensions[2] - 1) * (dimensions[1] - 1))) >= 0) {
							intersected_cells.insert(cell_idx - 1 - ((dimensions[2] - 1) * (dimensions[1] - 1)));
						}
						//cell_edge_map[cell_idx].push_back(edge_key);
					}
				}
				// negative y direction
				dir_idx = grid_point_idx - dimensions[2];
				if (dir_idx >= 0) {
					if (value_map.find(dir_idx) != value_map.end()) {
						dir_value = value_map[dir_idx];
					}
					else {
						dir_value = default_value;
					}
					if (current_value * dir_value <= isovalue) {
						MC::Vector p = calculate_interpolation_sparse(positions, dir_value, current_value, dir_idx, grid_point_idx);
						uint64_t edge_key = (uint64_t)(3 * dir_idx + 1);
						if (edge_vertex_map.find(edge_key) == edge_vertex_map.end()) {
							edge_vertex_map[edge_key] = intersection_points.size();
							intersection_points.push_back(p);
						}
						intersected_cells.insert(cell_idx - (dimensions[2] - 1));
						if (cell_idx - 1 - (dimensions[2] - 1) >= 0) {
							intersected_cells.insert(cell_idx - 1 - (dimensions[2] - 1));
						}
						if ((cell_idx - (dimensions[2] - 1) - ((dimensions[2] - 1) * (dimensions[1] - 1))) >= 0) {
							intersected_cells.insert(cell_idx - (dimensions[2] - 1) - ((dimensions[2] - 1) * (dimensions[1] - 1)));
						}
						if ((cell_idx - 1 - (dimensions[2] - 1) - ((dimensions[2] - 1) * (dimensions[1] - 1))) >= 0) {
							intersected_cells.insert(cell_idx - 1 - (dimensions[2] - 1) - ((dimensions[2] - 1) * (dimensions[1] - 1)));
						}
						//cell_edge_map[cell_idx - (dimensions[2] - 1)].push_back(edge_key);
					}
				}
				// positive x direction
				dir_idx = grid_point_idx + dimensions[2] * dimensions[1];
				if (dir_idx < positions.size()) {
					if (value_map.find(dir_idx) != value_map.end()) {
						dir_value = value_map[dir_idx];
					}
					else {
						dir_value = default_value;
					}
					if (current_value * dir_value <= isovalue) {
						MC::Vector p = calculate_interpolation_sparse(positions, current_value, dir_value, grid_point_idx, dir_idx);
						uint64_t edge_key = (uint64_t)(3 * grid_point_idx);
						if (edge_vertex_map.find(edge_key) == edge_vertex_map.end()) {
							edge_vertex_map[edge_key] = intersection_points.size();
							intersection_points.push_back(p);
						}
						intersected_cells.insert(cell_idx);
						if (cell_idx - 1 >= 0) {
							intersected_cells.insert(cell_idx - 1);
						}
						if ((cell_idx - (dimensions[2] - 1)) >= 0) {
							intersected_cells.insert(cell_idx - (dimensions[2] - 1));
						}
						if ((cell_idx - 1 - (dimensions[2] - 1)) >= 0) {
							intersected_cells.insert(cell_idx - 1 - (dimensions[2] - 1));
						}
						//cell_edge_map[cell_idx].push_back(edge_key);
					}
				}
				// negative x direction
				dir_idx = grid_point_idx - dimensions[2] * dimensions[1];
				if (dir_idx >= 0) {
					if (value_map.find(dir_idx) != value_map.end()) {
						dir_value = value_map[dir_idx];
					}
					else {
						dir_value = default_value;
					}
					if (current_value * dir_value <= isovalue) {
						MC::Vector p = calculate_interpolation_sparse(positions, dir_value, current_value, dir_idx, grid_point_idx);
						uint64_t edge_key = (uint64_t)(3 * dir_idx);
						if (edge_vertex_map.find(edge_key) == edge_vertex_map.end()) {
							edge_vertex_map[edge_key] = intersection_points.size();
							intersection_points.push_back(p);
						}
						intersected_cells.insert(cell_idx - ((dimensions[2] - 1) * (dimensions[1] - 1)));
						if (cell_idx - 1 - ((dimensions[2] - 1) * (dimensions[1] - 1)) >= 0) {
							intersected_cells.insert(cell_idx - 1 - ((dimensions[2] - 1) * (dimensions[1] - 1)));
						}
						if ((cell_idx - (dimensions[2] - 1) - ((dimensions[2] - 1) * (dimensions[1] - 1))) >= 0) {
							intersected_cells.insert(cell_idx - (dimensions[2] - 1) - ((dimensions[2] - 1) * (dimensions[1] - 1)));
						}
						if ((cell_idx - 1 - (dimensions[2] - 1) - ((dimensions[2] - 1) * (dimensions[1] - 1))) >= 0) {
							intersected_cells.insert(cell_idx - 1 - (dimensions[2] - 1) - ((dimensions[2] - 1) * (dimensions[1] - 1)));
						}
						//cell_edge_map[cell_idx - ((dimensions[2] - 1) * (dimensions[1] - 1))].push_back(edge_key);
					}
				}
			}
		}

		for (auto cell : intersected_cells) {
			int cell_idx = cell;
			int x = cell_idx / ((dimensions[2] - 1) * (dimensions[1] - 1));
			int y = (cell_idx - (x * (dimensions[2] - 1) * (dimensions[1] - 1))) / (dimensions[2] - 1);
			int z = cell_idx - (x * (dimensions[2] - 1) * (dimensions[1] - 1)) - (y * (dimensions[2] - 1));
			int vertex_idx = (dimensions[2] * dimensions[1] * x) + (dimensions[2] * y) + z;
			//std::vector<uint64_t> cell_edges = cell_edge_map[cell_idx];
			int cube_mask = 0;
			std::array<int, 8> cell_vertex_indices = { 0,0,0,0,0,0,0,0 };
			for (int v = 0; v < CELL_VERTICES.size(); v++) {
				std::array<int, 3> cell_vertices = CELL_VERTICES[v];
				int idx = (dimensions[2] * dimensions[1] * (x + cell_vertices[0])) + (dimensions[2] * (y + cell_vertices[1])) + z + cell_vertices[2];
				if (value_map.find(idx) != value_map.end()) {
					cube_mask |= value_map[idx] <= isovalue ? (1 << v) : 0;
				}
				else {
					cube_mask |= default_value <= isovalue ? (1 << v) : 0;
				}
				cell_vertex_indices[v] = idx;
			}
			if (cube_mask == 0 || cube_mask == 255) {
				continue;
			}
			for (unsigned int i = 0; i < 5; i++) {
				std::array<int, 3> tri = MARCHING_CUBES_TABLE[cube_mask][i];
				if (tri[0] != -1) {
					std::array<int, 3> triangle = { 0, 0, 0 };
					for (int j = 0; j < 3; j++) {
						int edge_idx = tri[j];
						std::array<int, 2> edge = CELL_EDGES[edge_idx];
						int idx_1 = cell_vertex_indices[edge[0]];
						int idx_2 = cell_vertex_indices[edge[1]];
						int idx_small = std::min(idx_1, idx_2);
						uint64_t edge_key = (uint64_t)(3 * idx_small + CELL_EDGES_DIRECTION[edge_idx]);
						int vertex_idx = edge_vertex_map.at(edge_key);
						triangle[j] = vertex_idx;
					}
					triangles.push_back(triangle);
				}
			}
		}
		vertices.resize(intersection_points.size());
		vertices = intersection_points;
	}

	void MarchingCubes::calculate_ours(std::vector<MC::Vector>& positions, std::vector<CompactNSearch::Real>& values, std::array<int, 3> dimensions, std::vector<MC::Vector>& vertices, std::vector<std::array<int, 3>>& triangles) {
		int vertex_idx = 0;

		#pragma omp parallel for num_threads(mcparameters.num_threads) schedule(static)
		for (int x = 0; x < dimensions[0] - 1; x++) {
			#pragma omp parallel for num_threads(mcparameters.num_threads) schedule(static)
			for (int y = 0; y < dimensions[1] - 1; y++) {
				#pragma omp parallel for num_threads(mcparameters.num_threads) schedule(static)
				for (int z = 0; z < dimensions[2] - 1; z++) {

					int cube_mask = 0;
					unsigned int cell[8] = { 0,0,0,0,0,0,0,0 };
					//Calculate cube mask
					for (unsigned int i = 0; i < CELL_VERTICES.size(); i++) {
						std::array<int, 3> cell_vertices = CELL_VERTICES[i];
						int idx = (dimensions[2] * dimensions[1] * (x + cell_vertices[0])) + (dimensions[2] * (y + cell_vertices[1])) + z + cell_vertices[2];
						cube_mask |= values[idx] < isovalue ? (1 << i) : 0;
						cell[i] = idx;
					}

					if (cube_mask == 0 || cube_mask == 255) {
						continue;
					}
					std::vector<MC::Vector> points;
					points.resize(15);

					//Interpolation
					int edge_mask = MARCHING_CUBES_EDGE_TABLE[cube_mask];
					for (unsigned int i = 0; i < CELL_EDGES.size(); i++) {
						if ((edge_mask & (1 << i)) > 0) {
							MC::Vector p;
							std::array<int, 2> edge = CELL_EDGES[i];
							CompactNSearch::Real v1 = values[cell[edge[0]]];
							MC::Vector p1 = positions[cell[edge[0]]];

							CompactNSearch::Real v2 = values[cell[edge[1]]];
							MC::Vector p2 = positions[cell[edge[1]]];

							CompactNSearch::Real mu;
							if (std::abs(isovalue - v1) < 0.00001) {
								p = p1;
							}
							else if (std::abs(isovalue - v2) < 0.00001) {
								p = p2;
							}
							else if (std::abs(v1 - v2) < 0.00001) {
								p = p1;
							}
							else {
								mu = (isovalue - v1) / (v2 - v1);
								p = p1 + mu * (p2 - p1);
							}
							points[i] = p;
						}
					}

					//Add vertices and triangle indices
					for (unsigned int i = 0; i < 5; i++) {
						std::array<int, 3> tri = MARCHING_CUBES_TABLE[cube_mask][i];
						if (tri[0] != -1) {
							MC::Vector v0 = points[tri[0]];
							MC::Vector v1 = points[tri[1]];
							MC::Vector v2 = points[tri[2]];

							#pragma omp critical
							{
								vertices.push_back(v0);
								vertices.push_back(v1);
								vertices.push_back(v2);

								std::array<int, 3> triangle = { vertex_idx, vertex_idx + 1,vertex_idx + 2 };
								triangles.push_back(triangle);
								vertex_idx += 3;
							}
						}
					}
				}
			}
		}
	}
}
