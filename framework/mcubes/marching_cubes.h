#pragma once

#include <vector>
#include "config.h"
#include "CompactNSearch/Config.h"
#include "mcparams.h"

namespace MC
{
	class MarchingCubes
	{
	public:
		MarchingCubes(MCParameters mcparams, CompactNSearch::Real value = 0.0);
		virtual ~MarchingCubes() = default;
		void calculate(std::vector<MC::Vector>& positions, std::vector<CompactNSearch::Real>& values, std::unordered_map<uint64_t, CompactNSearch::Real>& value_map, std::array<int, 3> dimensions, std::vector<MC::Vector>& vertices, std::vector<std::array<int, 3>>& triangles);

	private:
		MCParameters mcparameters;
		CompactNSearch::Real isovalue;
		void calculate_regular(std::vector<MC::Vector>& positions, std::vector<CompactNSearch::Real>& values, std::array<int, 3> dimensions, std::vector<MC::Vector>& vertices, std::vector<std::array<int, 3>>& triangles);
		void calculate_sparse(std::vector<MC::Vector>& positions, std::unordered_map<uint64_t, CompactNSearch::Real>& value_map, std::array<int, 3> dimensions, std::vector<MC::Vector>& vertices, std::vector<std::array<int, 3>>& triangles);
		void calculate_ours(std::vector<MC::Vector>& positions, std::vector<CompactNSearch::Real>& values, std::array<int, 3> dimensions, std::vector<MC::Vector>& vertices, std::vector<std::array<int, 3>>& triangles);
		MC::Vector calculate_interpolation(std::vector<MC::Vector>& positions, std::vector<CompactNSearch::Real>& values, int starting_index, int ending_index);
		MC::Vector calculate_interpolation_sparse(std::vector<MC::Vector>& positions, CompactNSearch::Real v1, CompactNSearch::Real v2, int starting_index, int ending_index);
	};
};