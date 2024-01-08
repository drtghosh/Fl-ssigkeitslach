#pragma once
#include <cassert>
#include <Eigen/Dense>
#include <cmath>

namespace kernel
{
	class Kernel
	{
	public:
		Kernel(const double smoothing_length);
		virtual ~Kernel() = default;

		double cubic_spline(const Eigen::Vector3d x_i, const Eigen::Vector3d x_j);
		Eigen::Vector3d cubic_grad_spline(const Eigen::Vector3d x_i, const Eigen::Vector3d x_j);
		double cohesion_kernel(const Eigen::Vector3d x_i, const Eigen::Vector3d x_j);
		double adhesion_kernel(const Eigen::Vector3d x_i, const Eigen::Vector3d x_j);

	private:
		double smoothing_length;
		double inverse_h;
		double inverse_h3;
		double inverse_h4;
		double compact_support;
		double c3;
		double c6;
		double inverse_c9;
		double inverse_c325;
	};
};