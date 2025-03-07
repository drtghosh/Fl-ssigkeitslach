#pragma once
#include <Eigen/Dense>

namespace PBD
{
#ifdef USE_DOUBLE
	using Vector = Eigen::Vector3d;
#else
	using Vector = Eigen::Vector3f;
#endif
}