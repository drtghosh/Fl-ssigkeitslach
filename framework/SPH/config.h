#pragma once
#include <Eigen/Dense>

namespace WCSPH
{
#ifdef USE_DOUBLE
	using Vector = Eigen::Vector3d;
#else
	using Real = Eigen::Vector3f;
#endif
}