#pragma once

#include <vector>
#include <set>
#include <math.h> 

namespace BruteForceNSearch
{
	class BruteForceNSearch
	{

	public:

		/**
		* Constructor.
		* Creates a new instance of the neighborhood search class.
		* @param r Search radius. If two points are closer to each other than a distance r they are considered neighbors.
		* @param erase_empty_cells If true. Empty cells in spatial hashing grid are erased if the points move.
		*/
		BruteForceNSearch(float radius);

		/**
		* Destructor.
		*/
		virtual ~BruteForceNSearch() = default;

		void add_points(std::vector<float> points);

		std::vector<float> get_points();

		std::vector<std::set<int>> get_neighbor_indices();

		void refresh_neighbor_indices();

		void find_neighbors();

		void find_neighbors_2();

		void find_neighbors_parallel();

	private:


		std::vector<float> points;

		std::vector<std::set<int>> neighbor_indices;

		float radius{ 0.0F };

		float calculate_distance(float x1, float y1, float z1, float x2, float y2, float z2);
	};
}
