#pragma once
#include <cassert>
#include <Eigen/Dense>

namespace kernel
{
	class Kernel
	{
	public:
		Kernel(const double smoothing_length);
		virtual ~Kernel() = default;

		double cubic_spline(const Eigen::Vector3d x_i, const Eigen::Vector3d x_j);
		Eigen::Vector3d cubic_grad_spline(const Eigen::Vector3d x_i, const Eigen::Vector3d x_j);

	private:
		double smoothing_length;
		double inverse_h;
		double inverse_h3;
		double inverse_h4;

	};
};