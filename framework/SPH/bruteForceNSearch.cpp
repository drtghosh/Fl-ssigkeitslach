#include "bruteForceNSearch.h"

#include <omp.h>

namespace BruteForceNSearch
{

	BruteForceNSearch::BruteForceNSearch(float radius) {
		this->radius = radius;
	}

	void BruteForceNSearch::add_points(std::vector<float> points) {
		if (points.size() % 3 != 0) {
			throw std::exception("Input not in correct format");
		}
		this->points.reserve(points.size());
		this->neighbor_indices.reserve(points.size()/3);
		for (unsigned int i = 0; i < points.size(); i++) {
			this->points.emplace_back(points.at(i));
			if (i % 3 == 0) {
				this->neighbor_indices.emplace_back(std::set<int>());
			}
		}
	}

	std::vector<float> BruteForceNSearch::get_points() {
		return this->points;
	}

	std::vector<std::set<int>> BruteForceNSearch::get_neighbor_indices() {
		return this->neighbor_indices;
	}

	void BruteForceNSearch::refresh_neighbor_indices() {
		this->neighbor_indices.clear();
		for (unsigned int i = 0; i < this->points.size(); i+=3) {
			this->neighbor_indices.emplace_back(std::set<int>());
		}
	}

	void BruteForceNSearch::find_neighbors() {
		refresh_neighbor_indices();

		for (unsigned int i = 0; i < this->points.size(); i+=3) {
			for (unsigned int j = 0; j < this->points.size(); j+=3) {
				float distance = calculate_distance(this->points.at(i), this->points.at(i+1), this->points.at(i+2), 
					this->points.at(j), this->points.at(j+1), this->points.at(j+2));

				if (distance > this->radius) {
				//if (distance >= this->radius) {
					continue;
				}
				int index_1 = i / 3;
				int index_2 = j / 3;

				this->neighbor_indices.at(index_1).insert(index_2);
			}
		}
	}

	void BruteForceNSearch::find_neighbors_2() {
		refresh_neighbor_indices();

		for (unsigned int i = 0; i < this->points.size(); i += 3) {
			for (unsigned int j = i; j < this->points.size(); j += 3) {
				float distance = calculate_distance(this->points.at(i), this->points.at(i + 1), this->points.at(i + 2),
					this->points.at(j), this->points.at(j + 1), this->points.at(j + 2));

				if (distance > this->radius) {
				//if (distance >= this->radius) {
					continue;
				}
				int index_1 = i / 3;
				int index_2 = j / 3;

				this->neighbor_indices.at(index_1).insert(index_2);
				this->neighbor_indices.at(index_2).insert(index_1);
			}
		}
	}

	void BruteForceNSearch::find_neighbors_parallel() {
		refresh_neighbor_indices();

		#pragma omp parallel for num_threads(10) schedule(static)
		for (int i = 0; i < this->points.size(); i += 3) {
			for (unsigned int j = i; j < this->points.size(); j += 3) {
				float distance = calculate_distance(this->points.at(i), this->points.at(i + 1), this->points.at(i + 2),
					this->points.at(j), this->points.at(j + 1), this->points.at(j + 2));

				if (distance > this->radius) {
					//if (distance >= this->radius) {
					continue;
				}
				int index_1 = i / 3;
				int index_2 = j / 3;

				this->neighbor_indices.at(index_1).insert(index_2);
				this->neighbor_indices.at(index_2).insert(index_1);
			}
		}
	}

	float BruteForceNSearch::calculate_distance(float x1, float y1, float z1, float x2, float y2, float z2) {
		return (float)sqrt(pow(x2 - x1,2)+ pow(y2 - y1, 2)+pow(z2 - z1, 2));
	}
}