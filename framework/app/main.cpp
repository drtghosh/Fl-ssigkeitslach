#include <stdlib.h>     // rand
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>    // std::max

#include <Eigen/Dense>

#include "../geometry/io.h"
#include "../geometry/sampling.h"

int main()
{
	std::cout << "Welcome to the learnSPH framework!!" << std::endl;
	std::cout << "Generating a sample scene...";

	const std::vector<geometry::TriMesh> meshes1 = geometry::read_tri_meshes_from_obj("../res/ItalianFountain.obj");
	for (int i = 0; i < (int)meshes1.size(); i++) {
		std::cout << "mesh " << i << std::endl;
		std::cout << "vertices: " << meshes1[i].vertices.size() << std::endl;
		std::cout << "triangles: " << meshes1[i].triangles.size() << std::endl;
	}
	const geometry::TriMesh& obj = meshes1[12];
	std::vector<Eigen::Vector3d> vertices = obj.vertices;
	const std::vector<std::array<int, 3>> triangles = obj.triangles;
	for (int i = 0; i < (int)vertices.size(); i++) {
		std::cout << vertices[i][2] << std::endl;
	}
	for (int i = 0; i < (int)triangles.size(); i++) {
		std::cout << triangles[i][0] << ' ' << triangles[i][1] << ' ' << triangles[i][2] << std::endl;
	}


	// Load a obj surface mesh
	const std::vector<geometry::TriMesh> meshes = geometry::read_tri_meshes_from_obj("../res/box.obj");
	const geometry::TriMesh& box = meshes[0];

	// Sample the mesh with particles
	const double sampling_distance = 0.05;
	std::vector<Eigen::Vector3d> particles;
	geometry::sampling::triangle_mesh(particles, box.vertices, box.triangles, sampling_distance);

	// Initialize data vectors
	std::vector<double> particles_scalar_data(particles.size());
	std::vector<Eigen::Vector3d> particles_vector_data(particles.size());

	// Scalar data will be the particle_id
	for (int i = 0; i < (int)particles.size(); i++) {
		particles_scalar_data[i] = i;
	}
	
	// Simulation loop
	for (int time_step = 0; time_step < 100; time_step++) {

		for (int particle_i = 0; particle_i < (int)particles.size(); particle_i++) {
			// Move particles a bit down in the Z direction
			particles[particle_i][2] -= 0.025;

			// Clamp the Z coord to the floor
			particles[particle_i][2] = std::max(particles[particle_i][2], -1.0);
			
			// Vector data is going to be the position
			particles_vector_data[particle_i] = particles[particle_i] - Eigen::Vector3d(-1, -1, -1);
		}

		// Save output
		const std::string filename = "../res/example_" + std::to_string(time_step) + ".vtk";
		geometry::write_particles_to_vtk(filename, particles, particles_scalar_data, particles_vector_data);
	}

	std::cout << "completed!" << std::endl;
	std::cout << "The scene files have been saved in the folder `<build_folder>/res`. \nYou can visualize them with Blender by using the Blender Sequence Loaded addon." << std::endl;

	return 0;
}