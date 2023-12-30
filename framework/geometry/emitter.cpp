#include "emitter.h"

namespace Emitter {

	Emitter::Emitter(CompactNSearch::Real radius1,geometry::Vector radius1_dir, CompactNSearch::Real radius2, geometry::Vector radius2_dir ,geometry::Vector centre, geometry::Vector normal_dir, const geometry::Vector velocity, const CompactNSearch::Real particle_dist, const CompactNSearch::Real shield_distance) : r1(radius1), r1_dir(radius1_dir), r2(radius2), r2_dir(radius2_dir), centre_position(centre), normal(normal_dir), velocity(velocity), particle_distance(particle_dist), distance(shield_distance) {
		assert(r1_dir.dot(r2_dir) == 0.0);
		assert(r1_dir.dot(normal) == 0.0);
		assert(r2_dir.dot(normal) == 0.0);

		r1_dir = r1_dir.normalized();
		r2_dir = r2_dir.normalized();
		normal = normal.normalized();

		calculate_shield();
	}

	void Emitter::set_schedule(std::vector<CompactNSearch::Real> schedule_times) {
		schedule = schedule_times;
	}

	void Emitter::calculate_shield() {
		geometry::Vector corner1 = centre_position + (r1 + distance) * r1_dir + (r2 + distance) * r2_dir + distance * -normal;
		geometry::Vector corner2 = centre_position + (r1 + distance) * -r1_dir + (r2 + distance) * r2_dir + distance * -normal;
		geometry::Vector corner3 = centre_position + (r1 + distance) * -r1_dir + (r2 + distance) * -r2_dir + distance * -normal;
		geometry::Vector corner4 = centre_position + (r1 + distance) * r1_dir + (r2 + distance) * -r2_dir + distance * -normal;

		shield_vertices.push_back(corner1);
		shield_vertices.push_back(corner2);
		shield_vertices.push_back(corner3);
		shield_vertices.push_back(corner4);

		shield_triangles.push_back({ 0,1,2 });
		shield_triangles.push_back({ 2,3,0 });

		corner3 = centre_position + (r1 + distance) * -r1_dir + (r2 + distance) * r2_dir + distance * normal;
		corner4 = centre_position + (r1 + distance) * r1_dir + (r2 + distance) * r2_dir + distance * normal;

		shield_vertices.push_back(corner3);
		shield_vertices.push_back(corner4);

		shield_triangles.push_back({ 0,1,4 });
		shield_triangles.push_back({ 4,5,0 });

		corner3 = centre_position + (r1 + distance) * -r1_dir + (r2 + distance) * -r2_dir + distance * normal;
		corner4 = centre_position + (r1 + distance) * -r1_dir + (r2 + distance) * r2_dir + distance * normal;

		shield_vertices.push_back(corner3);
		shield_vertices.push_back(corner4);

		shield_triangles.push_back({ 1,2,6 });
		shield_triangles.push_back({ 6,7,1 });

		corner3 = centre_position + (r1 + distance) * r1_dir + (r2 + distance) * -r2_dir + distance * normal;
		corner4 = centre_position + (r1 + distance) * -r1_dir + (r2 + distance) * -r2_dir + distance * normal;

		shield_vertices.push_back(corner3);
		shield_vertices.push_back(corner4);

		shield_triangles.push_back({ 2,3,8 });
		shield_triangles.push_back({ 8,9,2 });

		corner3 = centre_position + (r1 + distance) * r1_dir + (r2 + distance) * r2_dir + distance * normal;
		corner4 = centre_position + (r1 + distance) * r1_dir + (r2 + distance) * -r2_dir + distance * normal;

		shield_vertices.push_back(corner3);
		shield_vertices.push_back(corner4);

		shield_triangles.push_back({ 3,0,10 });
		shield_triangles.push_back({ 10,11,3 });
	}

	std::vector<std::array<int, 3 >> Emitter::get_shield_triangles() {
		return shield_triangles;
	}

	std::vector<geometry::Vector> Emitter::get_shield_vertices() {
		return shield_vertices;
	}

	void Emitter::emit_particles(std::vector<geometry::Vector>& particle_positions, std::vector<geometry::Vector>& particle_velocities, CompactNSearch::Real t_sim) {
		//Find points on xy plane that are inside ellipse on xy plane
		bool is_scheduled = false;
		for (unsigned int i = 0; i < schedule.size()-1; i+=2) {
			if (t_sim >= schedule[i] && t_sim < schedule[i + 1]) {
				is_scheduled = true;
			}
		}


		if (is_scheduled) {
			CompactNSearch::Real z = 0.0;
			for (CompactNSearch::Real x = -r1; x <= r1; x += particle_distance) {
				for (CompactNSearch::Real y = -r2; y <= r2; y += particle_distance) {
					bool is_inside_ellipse = (x * x) / (r1 * r1) + (y * y) / (r2 * r2) <= 1;
					if (is_inside_ellipse) {
						//Project onto 3d ellipse
						geometry::Vector point(x, y, z);
						geometry::Vector transformed_point(x, y, z);
						transformed_point[0] = point[0] * r1_dir[0] + point[1] * r2_dir[0] + point[2] * normal[0];
						transformed_point[1] = point[0] * r1_dir[1] + point[1] * r2_dir[1] + point[2] * normal[1];
						transformed_point[2] = point[0] * r1_dir[2] + point[1] * r2_dir[2] + point[2] * normal[2];
						transformed_point += centre_position;

						particle_positions.push_back(transformed_point);
						particle_velocities.push_back(velocity);
					}

				}
			}
		}
	}

	std::vector<geometry::Vector> Emitter::get_emition_positions() {
		std::vector<geometry::Vector> particle_positions;
		CompactNSearch::Real z = 0.0;
		for (CompactNSearch::Real x = -r1; x <= r1; x += particle_distance) {
			for (CompactNSearch::Real y = -r2; y <= r2; y += particle_distance) {
				bool is_inside_ellipse = (x * x) / (r1 * r1) + (y * y) / (r2 * r2) <= 1;
				if (is_inside_ellipse) {
					//Project onto 3d ellipse
					geometry::Vector point(x, y, z);
					geometry::Vector transformed_point(x, y, z);
					transformed_point[0] = point[0] * r1_dir[0] + point[1] * r2_dir[0] + point[2] * normal[0];
					transformed_point[1] = point[0] * r1_dir[1] + point[1] * r2_dir[1] + point[2] * normal[1];
					transformed_point[2] = point[0] * r1_dir[2] + point[1] * r2_dir[2] + point[2] * normal[2];
					transformed_point += centre_position;

					particle_positions.push_back(transformed_point);
				}
			}
		}
		return particle_positions;
	}
}