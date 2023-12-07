#include "catch.hpp"
#include <array>
#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include "../SPH/bruteForceNSearch.h"
#include <CompactNSearch/include/CompactNSearch/CompactNSearch.h>


// Check out https://github.com/catchorg/Catch2 for more information about how to use Catch2
TEST_CASE("Performance Test", "[Performance]")
{
	//initialize variables
	srand(42);
	unsigned int point_number = 100;
	float max_dim = 10000;
	float radius = 2;
	std::vector<float> points;

	//create random point cloud with coords [0,max_dim]
	for (unsigned int i = 0; i < point_number*3; i++) {
		float r = ((float)rand() / (RAND_MAX));
		r = r * max_dim;

		points.emplace_back(r);
		//std::cout << "############ " << r << std::endl;
	}
	
	//Brute force performance
	BruteForceNSearch::BruteForceNSearch b_search(radius);
	b_search.add_points(points);

	auto start = std::chrono::system_clock::now();
	b_search.find_neighbors();
	auto end = std::chrono::system_clock::now();
	auto elapsed =
		std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "Brute Force time: " << elapsed.count() << " ms" << std::endl;

	const std::vector<std::set<int>> b_neighbors = b_search.get_neighbor_indices();

	start = std::chrono::system_clock::now();
	b_search.find_neighbors_2();
	end = std::chrono::system_clock::now();
	elapsed =
		std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "Brute Force 2 time: " << elapsed.count() << " ms" << std::endl;

	const std::vector<std::set<int>> b2_neighbors = b_search.get_neighbor_indices();

	/*for (unsigned int i = 0; i < point_number; i++) {
		for (auto iter = b_neighbors.at(i).begin(); iter != b_neighbors.at(i).end(); iter++) {
			std::cout << "############ " << i << " " << *iter << std::endl;
		}
		std::cout << "###############################" << std::endl;
		for (auto iter = b2_neighbors.at(i).begin(); iter != b2_neighbors.at(i).end(); iter++) {
			std::cout << "############ " << i << " " << *iter << std::endl;
		}
		std::cout << "###############################" << std::endl;
		std::cout << "###############################" << std::endl;
	}*/

	REQUIRE(b_neighbors == b2_neighbors);

	start = std::chrono::system_clock::now();
	b_search.find_neighbors_parallel();
	end = std::chrono::system_clock::now();
	elapsed =
		std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "Brute Force 2 parallel time: " << elapsed.count() << " ms" << std::endl;

	//CompactNSearch performance
	CompactNSearch::NeighborhoodSearch c_search(radius);
	std::vector<Eigen::Vector3d> point_set;

	for (unsigned int i = 0; i < point_number; i++) {
		Eigen::Vector3d point = { points.at(3 * i), points.at(3 * i + 1), points.at(3 * i + 2) };
		point_set.emplace_back(point);
		//std::cout << "point: " << point.at(0)<<", "<<point.at(1)<<", "<<point.at(2) << std::endl;
	}

	c_search.add_point_set(point_set.front().data(), point_number);

	//c_search.z_sort();
	//c_search.point_set(0).sort_field(points.data());

	start = std::chrono::system_clock::now();
	c_search.find_neighbors();
	end = std::chrono::system_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	std::cout << "CompactNSearch time: " << elapsed.count() << " ms" << std::endl;

}