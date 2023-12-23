#include "kernel.h"

constexpr double PI = 3.14159265358979323846;
constexpr double alpha = 3.0 / (2.0 * PI);
constexpr double cohesion_alpha = 32.0 / PI;
const Eigen::Vector3d zero_grad = Eigen::Vector3d(0.0, 0.0, 0.0);

namespace kernel
{
	Kernel::Kernel(const double smoothing_length) {
		this->smoothing_length = smoothing_length;
		this->inverse_h = 1.0 / smoothing_length;
		this->inverse_h3 = 1.0 / (smoothing_length * smoothing_length * smoothing_length);
		this->inverse_h4 = 1.0 / (smoothing_length * smoothing_length * smoothing_length * smoothing_length);
		this->compact_support = 2.0 * smoothing_length;
		this->c3 = compact_support * compact_support * compact_support;
		this->inverse_c9 =	1.0 / (c3 * c3 * c3);
		this->c6 = c3 * c3;
		this->inverse_c325 = 1.0 / pow(compact_support, 3.25);
	}

	double Kernel::cubic_spline(const Eigen::Vector3d x_i, const Eigen::Vector3d x_j)
	{
		const double norm = (x_i - x_j).norm();
		//if (norm < std::numeric_limits<double>::epsilon()) {
			//return 0.0;
		//}
		const double q = norm * inverse_h;
		assert(q >= 0.0);

		if (q >= 2.0) {
			return 0.0;
		}
		if (q < 1.0) {
			return alpha * (2.0 / 3.0 - q * q + 0.5 * q * q * q) * inverse_h3;
		}
		else if (q < 2.0) {
			double a = 2.0 - q;
			return alpha * (a * a * a) * inverse_h3 / 6.0;
		}
	}

	Eigen::Vector3d Kernel::cubic_grad_spline(const Eigen::Vector3d x_i, const Eigen::Vector3d x_j)
	{
		const Eigen::Vector3d gradq = x_i - x_j;
		const double distance = (x_i - x_j).norm();

		if (distance <= std::numeric_limits<double>::epsilon()) {
			return zero_grad;
		}

		const double q = distance * inverse_h;
		assert(q >= 0.0);

		if (q >= 2.0) {
			return zero_grad;
		}
		if (q < 1.0) {
			return (alpha * (-2.0 * q + 1.5 * q * q) * inverse_h4 * gradq / distance) + zero_grad;
		}
		else if (q < 2.0) {
			double a = 2.0 - q;
			return (alpha * (-0.5 * a * a) * inverse_h4 * gradq / distance) + zero_grad;
		}
	}

	double Kernel::cohesion_kernel(const Eigen::Vector3d x_i, const Eigen::Vector3d x_j)
	{
		const double r = (x_i - x_j).norm();
		
		if (r >= compact_support) {
			return 0.0;
		}
		if (r < 0.5 * compact_support) {
			double diff = compact_support - r;
			return cohesion_alpha * inverse_c9 * ((2 * diff * diff * diff * r * r * r) - (c6 / 64));
		}
		else if (r < compact_support) {
			double diff = compact_support - r;
			return cohesion_alpha * inverse_c9 * (diff * diff * diff * r * r * r);
		}
	}

	double Kernel::adhesion_kernel(const Eigen::Vector3d x_i, const Eigen::Vector3d x_j)
	{
		const double r = (x_i - x_j).norm();

		if (r >= compact_support || r <= smoothing_length) {
			return 0.0;
		}
		else if (r < compact_support && r > smoothing_length) {
			double adhesiveness = 6.0 * r - 2.0 * compact_support - (4.0 * r * r / compact_support);
			if (adhesiveness < 0.0) {
				return 0.0;
			}
			else {
				return 0.007 * inverse_c325 * sqrt(sqrt(adhesiveness));
			}
		}
	}
}
