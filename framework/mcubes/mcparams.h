#pragma once
#include <CompactNSearch/Config.h>

class MCParameters
{
public:
	unsigned int num_threads = 4;
	bool sparse = true;
	bool ours = false;
	CompactNSearch::Real initial_grid_value = 0.55;
};
