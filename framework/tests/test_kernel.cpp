#include "catch.hpp"
#include "../kernel/kernel.h"

using namespace kernel;
constexpr double PI = 3.14159265358979323846;

/*TEST_CASE("Cubic spline kernel", "[kernel]")
{
	const double h = 1.0;
	const double value_at_0 = 1.0 / (PI * h*h*h);
	const double value_at_1 = 1.0 / (4.0 * PI * h*h*h);
	const Eigen::Vector3d x_orig = Eigen::Vector3d(0.0, 0.0, 0.0);
	const Eigen::Vector3d e_x = Eigen::Vector3d(1.0, 0.0, 0.0);
	const Eigen::Vector3d e_y = Eigen::Vector3d(0.0, 1.0, 0.0);
	const Eigen::Vector3d e_z = Eigen::Vector3d(0.0, 0.0, 1.0);
	const Eigen::Vector3d x_far = Eigen::Vector3d(0.0, 0.0, 2.1 * h);
	Kernel kernel(h);

	srand((unsigned)time(NULL));
	const Eigen::Vector3d x_random = Eigen::Vector3d((rand() % 10)/10.0, (rand() % 10) / 10.0, (rand() % 10) / 10.0);

	SECTION("Kernel is symmetric") {
		REQUIRE(kernel.cubic_spline(x_orig, x_random) == kernel.cubic_spline(x_random, x_orig));
		REQUIRE(kernel.cubic_spline(x_orig, e_x) == kernel.cubic_spline(e_x, x_orig));
	}

	SECTION("Kernel is non-negative") {
		REQUIRE(kernel.cubic_spline(x_orig, x_random) >= 0.0);
		REQUIRE(kernel.cubic_spline(x_orig, e_y) >= 0.0);
	}

	SECTION("Kernel takes maximum value of 1 when h is not too small") {
		REQUIRE(kernel.cubic_spline(x_orig, x_random) <= 1.0);
		REQUIRE(kernel.cubic_spline(x_orig, e_z) <= 1.0);
		REQUIRE(kernel.cubic_spline(e_x, e_y) <= 1.0);
	}

	SECTION("Kernel is zero when distance is greater than 2h (compactness)") {
		REQUIRE(kernel.cubic_spline(x_orig, x_far) == 0.0);
	}

	SECTION("Kernel matching in critical domain values") {
		REQUIRE((kernel.cubic_spline(x_orig, x_orig) - value_at_0) <= std::numeric_limits<double>::epsilon());
		REQUIRE((kernel.cubic_spline(x_random, x_random) - value_at_0) <= std::numeric_limits<double>::epsilon());
		REQUIRE((kernel.cubic_spline(x_orig, e_x) - value_at_1) <= std::numeric_limits<double>::epsilon());
		REQUIRE((kernel.cubic_spline(x_orig, e_y) - value_at_1) <= std::numeric_limits<double>::epsilon());
		REQUIRE((kernel.cubic_spline(x_orig, e_z) - value_at_1) <= std::numeric_limits<double>::epsilon());
		REQUIRE(kernel.cubic_spline(x_orig, 2*e_x) == 0.0);
	}

	SECTION("Kernel shadows Dirac delta for small values") {
		const double small_h = 1e-6;
		Kernel small_h_kernel(small_h);
		REQUIRE(small_h_kernel.cubic_spline(x_orig, x_orig) > 1e12);
		REQUIRE(small_h_kernel.cubic_spline(x_orig, e_x) == 0.0);
		REQUIRE(small_h_kernel.cubic_spline(x_orig, x_far) == 0.0);
	}
}	


Eigen::Vector3d finite_diff_grad(const Eigen::Vector3d x_i, const Eigen::Vector3d x_j, const double h)
{
	const Eigen::Vector3d e_x = Eigen::Vector3d(1, 0, 0);
	const Eigen::Vector3d e_y = Eigen::Vector3d(0, 1, 0);
	const Eigen::Vector3d e_z = Eigen::Vector3d(0, 0, 1);

	const double epsilon = 1e-6;

	Kernel kernel(h);

	const double grad_x = (kernel.cubic_spline(x_i + epsilon * e_x, x_j) - kernel.cubic_spline(x_i - epsilon * e_x, x_j)) / (2.0 * epsilon);
	const double grad_y = (kernel.cubic_spline(x_i + epsilon * e_y, x_j) - kernel.cubic_spline(x_i - epsilon * e_y, x_j)) / (2.0 * epsilon);
	const double grad_z = (kernel.cubic_spline(x_i + epsilon * e_z, x_j) - kernel.cubic_spline(x_i - epsilon * e_z, x_j)) / (2.0 * epsilon);

	const Eigen::Vector3d grad = Eigen::Vector3d(grad_x, grad_y, grad_z);

	return grad;
}

TEST_CASE("Cubic spline kernel gradient", "[kernel]")
{
	const double h = 1.0;
	const Eigen::Vector3d x_orig = Eigen::Vector3d(0.0, 0.0, 0.0);
	const Eigen::Vector3d e_x = Eigen::Vector3d(1.0, 0.0, 0.0);
	const Eigen::Vector3d e_y = Eigen::Vector3d(0.0, 1.0, 0.0);
	const Eigen::Vector3d e_z = Eigen::Vector3d(0.0, 0.0, 1.0);
	const Eigen::Vector3d x_far = Eigen::Vector3d(0.0, 0.0, 2.1 * h);
	const Eigen::Vector3d e_x_grad_at_1 = 0.75 * e_x / (PI * h * h * h * h);
	Kernel kernel(h);

	srand((unsigned)time(NULL));
	const Eigen::Vector3d x_random = Eigen::Vector3d((rand() % 10) / 10.0, (rand() % 10) / 10.0, (rand() % 10) / 10.0);

	SECTION("Analytical gradient matches finite difference") {
		const Eigen::Vector3d grad_diff1 = kernel.cubic_grad_spline(x_orig, x_random) - finite_diff_grad(x_orig, x_random, h);
		const Eigen::Vector3d grad_diff2 = kernel.cubic_grad_spline(x_orig, e_x) - finite_diff_grad(x_orig, e_x, h);
		REQUIRE(std::abs(grad_diff1(0)) <= 1e-10);
		REQUIRE(std::abs(grad_diff1(1)) <= 1e-10);
		REQUIRE(std::abs(grad_diff1(2)) <= 1e-10);
		REQUIRE(std::abs(grad_diff2(0)) <= 1e-10);
		REQUIRE(std::abs(grad_diff2(1)) <= 1e-10);
		REQUIRE(std::abs(grad_diff2(2)) <= 1e-10);
	}

	SECTION("Gradient is zero when distance is greater than 2h") {
		REQUIRE(kernel.cubic_grad_spline(x_orig, x_far) == Eigen::Vector3d(0.0, 0.0, 0.0));
	}

	SECTION("Gradient is zero when distance is zero") {
		REQUIRE(kernel.cubic_grad_spline(x_orig, x_orig) == Eigen::Vector3d(0.0, 0.0, 0.0));
	}

	SECTION("Kernel gradient matching in critical domain values") {
		const Eigen::Vector3d grad_diff_at_1 = kernel.cubic_grad_spline(x_orig, e_x) - e_x_grad_at_1;
		REQUIRE(std::abs(grad_diff_at_1(0)) <= std::numeric_limits<double>::epsilon());
		REQUIRE(std::abs(grad_diff_at_1(1)) <= std::numeric_limits<double>::epsilon());
		REQUIRE(std::abs(grad_diff_at_1(2)) <= std::numeric_limits<double>::epsilon());
		REQUIRE(kernel.cubic_grad_spline(x_orig, 2 * e_x) == Eigen::Vector3d(0.0, 0.0, 0.0));
	}
}*/