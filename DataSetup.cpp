//	DataSetup.cpp
//
//  Copyright (C) 2015 Bailin Deng <bldeng@gmail.com>
//
//  This file is part of LBC - Local Barycentric Coordinates.
//
//	LBC is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	LBC is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with LBC. If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////

#include "DataSetup.h"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace LBC{

void DataSetup::init(WeightingScheme scheme)
{
    normalize_sample_points();
	compute_boundary_values();
    construct_operator_matrices();
    compute_geodesic_distance_and_init_inner_points_coordinates();
	compute_solver_data(scheme);
}


void DataSetup::compute_projection_barycentric_coordinates(const DenseMatrix& base_points, const DenseMatrix& pts, DenseMatrix& proj_bc)
{
    assert(base_points.rows() == pts.rows());

    // Set up a rectangular linear system \sum_{i=1}^{n-1} x_i (b_i - b_0) = p,
    // where b_0, b_1, ..., b_{n-1} are the columns of matrix base_points,
    // and x_1, ..., x_{n-1} are the barycentric coordinate entries corresponding to b_1, ..., b_{n-1}.
    // This linear system is solved in a least squares manner,
    // to obtain the projection point from p onto the subspace spanned by b_1, ..., b_{n-1}
    int m = base_points.rows(), n = base_points.cols();
    DenseMatrix A = base_points.block(0, 1, m, n-1);
    A.colwise() -= base_points.col(0);

    // RHS of the linear systems
    DenseMatrix B = pts;
    B.colwise() -= base_points.col(0);

    Eigen::JacobiSVD< DenseMatrix, Eigen::FullPivHouseholderQRPreconditioner > jsvd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    // Solve the least squares system to obtain the coordinates for b_1, ..., b_{n-1}
    DenseMatrix x = jsvd.solve(B);
    proj_bc.resize(n, pts.cols());
    proj_bc.block(1, 0, n-1, pts.cols()) = x;

    // Compute the first coordinate by subtract the sum of other coordinates from 1
    proj_bc.row(0).fill(1.0);
    proj_bc.row(0) -= x.colwise().sum();
}

void DataSetup::compute_barycentric_coordinates(DataSetup::CageBoundaryFacetInfo &info)
{
	if(info.has_coordinates){
		return;
	}

	int n_face_vertices = info.facet_vertices.size(), n_boundary_points = info.boundary_points.size();

	DenseMatrix A(sample_points_.rows(), n_face_vertices);
	for(int i = 0; i < n_face_vertices; ++ i){
		A.col(i) = sample_points_.col(control_point_idx_(info.facet_vertices(i)));
	}

	DenseMatrix B(sample_points_.rows(), n_boundary_points);
	for(int i = 0; i < n_boundary_points; ++ i){
		B.col(i) = sample_points_.col(info.boundary_points(i));
	}

    DenseMatrix bc;
    compute_projection_barycentric_coordinates(A, B, bc);
    info.barycentric_coordinates = bc.transpose();
}


void DataSetup::compute_boundary_values()
{
	const int INNER_POINT = 0, BOUNDARY_POINT = 1, CONTROL_POINT = 2;

	int n_sample_points = sample_points_.cols(), n_control_points = control_point_idx_.size();
	IndexVector point_type = IndexVector::Constant(n_sample_points, INNER_POINT);

	init_coordinate_values_.setZero(n_sample_points, n_control_points);

	for(int i = 0; i < n_control_points; ++ i){
		point_type(control_point_idx_(i)) = CONTROL_POINT;
		init_coordinate_values_(control_point_idx_(i), i) = 1.0;
	}

	for(int i = 0; i < static_cast<int>(boundary_facet_info_.size()); ++ i){

		CageBoundaryFacetInfo &current_facet_info = boundary_facet_info_[i];
		compute_barycentric_coordinates(current_facet_info);

		int n_facet_vertices = current_facet_info.facet_vertices.size(),
				n_boundary_points = current_facet_info.boundary_points.size();

		DenseMatrix control_point_coordinate_values(n_facet_vertices, n_control_points);
		for(int j = 0; j < n_facet_vertices; ++ j){
			control_point_coordinate_values.row(j) =
					init_coordinate_values_.row(control_point_idx_(current_facet_info.facet_vertices(j)));
		}

		DenseMatrix boundary_point_coordinates = current_facet_info.barycentric_coordinates * control_point_coordinate_values;
		for(int j = 0; j < n_boundary_points; ++j){
			int cur_boundary_pt_idx = current_facet_info.boundary_points(j);
			init_coordinate_values_.row(cur_boundary_pt_idx) = boundary_point_coordinates.row(j);
			point_type(cur_boundary_pt_idx) = BOUNDARY_POINT;
		}
	}

	std::vector<TripletD> triplets_inner_point_mapping, triplets_control_point_mapping;

	int inner_pt_idx = 0;
	for(int i = 0; i < n_sample_points; ++ i){
		if(point_type(i) == INNER_POINT){
			triplets_inner_point_mapping.push_back(TripletD(i, inner_pt_idx ++, 1.0));
		}
	}
	inner_points_mapping_.resize(n_sample_points, inner_pt_idx);
	inner_points_mapping_.setFromTriplets(triplets_inner_point_mapping.begin(), triplets_inner_point_mapping.end());

	for(int i = 0; i < n_control_points; ++ i){
		triplets_control_point_mapping.push_back(TripletD(control_point_idx_(i), i, 1.0));
	}
	control_points_mapping_.resize(n_sample_points, n_control_points);
	control_points_mapping_.setFromTriplets(triplets_control_point_mapping.begin(), triplets_control_point_mapping.end());
}

void DataSetup::compute_solver_data(WeightingScheme scheme)
{
	control_points_ = sample_points_ * control_points_mapping_;
	inner_points_ = sample_points_ * inner_points_mapping_;
	grad_operator_for_inner_points_ = integrated_grad_operator_ * inner_points_mapping_;
	grad_const_ = integrated_grad_operator_ * init_coordinate_values_;

	int n_faces = cell_vertices_.cols();
	int n_control_points = control_points_.cols();
	grad_weights_.resize(n_faces, n_control_points);

	if(scheme == CONSTANT){
		grad_weights_.fill(1.0);
	}
	else{
		grad_weights_.fill(0.0);

		// Compute the geodesic distance at the centroid of each cell, by averaging the distance values at its vertices
		for(int i = 0; i < n_faces; ++ i){
			for(int j = 0; j < cell_vertices_.rows(); ++ j){
				grad_weights_.row(i) += geodesic_distance_.row(cell_vertices_(j, i));
			}
		}

        grad_weights_ /= cell_vertices_.rows();
        grad_weights_ /= geodesic_distance_.maxCoeff();

		// Normalize the distance values for each control point

		for(int i = 0; i < n_control_points; ++ i){

			for(int j = 0; j < n_faces; ++ j){
				double dist_value = grad_weights_(j, i);

				switch(scheme){
				case SQUARE:
					dist_value = dist_value * dist_value;
					break;
				case SQUAREROOT:
					dist_value = std::sqrt(dist_value);
					break;
				default:
					break;
				}

				grad_weights_(j, i) = dist_value;
			}
		}


	}
}

void DataSetup::normalize_sample_points()
{
    DenseVector centroid = sample_points_.rowwise().mean();
    sample_points_.colwise() -= centroid;
    double orig_diam = sample_points_.colwise().norm().maxCoeff();
    //std::cout << "Original diameter " << orig_diam << std::endl;
    sample_points_ *= 250.0 / std::max(1e-10, orig_diam);
    //std::cout << "new diameter " << sample_points_.colwise().norm().maxCoeff() << std::endl;
}

void DataSetup::construct_operator_matrices()
{
	int n_row = cell_vertices_.rows();
	assert( n_row == 3 || n_row == 4 );

	int dim = n_row - 1;
	int n_cells = cell_vertices_.cols();
	int n_pts = sample_points_.cols();

	cell_grad_operator_.resize(dim * n_cells, n_pts);
	integrated_grad_operator_.resize(dim * n_cells, n_pts);
	vertex_associted_measure_.resize(n_pts, n_pts);


	std::vector<TripletD> triplets_grad_operator, triplets_vertex_measure, triplets_integrated_grad_opeartor;

	for(int i = 0; i < n_cells; ++i){
		IndexVector current_idx = cell_vertices_.col(i);
		assert(dim + 1 == static_cast<int>(current_idx.size()));

		DenseMatrix grad_coef, integrated_grad_coef;
		double cell_measure;

		compute_gradient_coefficients(current_idx, grad_coef, integrated_grad_coef, cell_measure);

		for(int k = 0; k <= dim; ++ k){
			for(int j = 0; j < dim; ++ j){
                triplets_grad_operator.push_back(TripletD(dim * i + j, current_idx(k), grad_coef(j, k)));
                triplets_integrated_grad_opeartor.push_back(TripletD(dim * i + j, current_idx(k), integrated_grad_coef(j, k)));
			}

			// Here we compute the contribution to the vertex measure for each cell vertex from the current cell.
			// Eigen accumulates the values in triplets when creating a sparse matrix.
			triplets_vertex_measure.push_back(TripletD(current_idx(k), current_idx(k), cell_measure / double(current_idx.size())));
		}
	}


	cell_grad_operator_.setFromTriplets(triplets_grad_operator.begin(), triplets_grad_operator.end());
	integrated_grad_operator_.setFromTriplets(triplets_integrated_grad_opeartor.begin(), triplets_integrated_grad_opeartor.end());
	vertex_associted_measure_.setFromTriplets(triplets_vertex_measure.begin(), triplets_vertex_measure.end());
    integrated_divergence_operator_ = -integrated_grad_operator_.transpose();
    symmetric_laplacian_operator_ = integrated_divergence_operator_ * cell_grad_operator_;
}

double DataSetup::compute_geodesic_time_step()
{
	// Use the squared maximum distance between neighboring vertices as the time step
	int n_cells = cell_vertices_.cols(), n_face_vtx = cell_vertices_.rows();
	DenseVector face_max_edge_length;
	face_max_edge_length.setZero(n_cells);

	#pragma omp parallel for
	for(int i = 0; i < n_cells; ++i){
		for(int j = 0; j < n_face_vtx; ++ j){
			for(int k = j + 1; k < n_face_vtx; ++ k){
				face_max_edge_length(i) = std::max(
						face_max_edge_length(i),
						(sample_points_.col(cell_vertices_(j, i)) - sample_points_.col(cell_vertices_(k, i))).norm()
					);
			}
		}
	}

	double h = face_max_edge_length.maxCoeff();

	return h * h;
}

void DataSetup::compute_geodesic_distance_and_init_inner_points_coordinates()
{
	int n_row = cell_vertices_.rows();
	assert( n_row == 3 || n_row == 4 );
	int dim = n_row - 1;
	int n_cells = cell_vertices_.cols();
	int n_sample_points = sample_points_.cols();
	int n_control_points = control_point_idx_.size();

	geodesic_distance_.setZero(n_sample_points, n_control_points);

	double t = compute_geodesic_time_step();

	// To solve the Poisson equation, we set the value of the last variable to zero,
	// and solve for the remaining variables as a linear least squares problem.
	// The following SparseMatrix is used for evaluating the right hand side of the
	// corresponding linear system from the normalized gradient field.
	SparseMatrix rhs_op;


    // We will solve a linear least squares system for computing harmonic coordinate values for the inner sample points
    // The following matrix is used for computing the rhs of this system
    SparseMatrix rhs_op_inner_vals;
    inner_point_init_values_.setZero(inner_points_mapping_.cols(), n_control_points);



	#pragma omp parallel sections
	{
		#pragma omp section
		{
			SparseMatrix M_heatflow = vertex_associted_measure_ - symmetric_laplacian_operator_ * t;
			solver_heatflow_.compute(M_heatflow);
		}

		#pragma omp section
		{
			// Remove the last variable and set it to 0, such that the problem becomes over-determined and has a unique solution
			SparseMatrix sub_mat = symmetric_laplacian_operator_.block(0, 0, symmetric_laplacian_operator_.rows(), symmetric_laplacian_operator_.cols() - 1);
			SparseMatrix sub_mat_T = sub_mat.transpose();
			SparseMatrix M_geodesics = sub_mat_T * sub_mat;
			rhs_op = sub_mat_T * integrated_divergence_operator_;
			solver_geodesics_.compute(M_geodesics);
		}

        #pragma omp section
        {
            SparseMatrix reduced_grad_op = cell_grad_operator_ * inner_points_mapping_;
            SparseMatrix reduced_grad_op_T = reduced_grad_op.transpose();
            SparseMatrix M_inner_init_vals = reduced_grad_op_T * reduced_grad_op;
            solver_inner_init_vals_.compute(M_inner_init_vals);
            rhs_op_inner_vals = - reduced_grad_op_T * cell_grad_operator_;
        }
	}

	if(solver_heatflow_.info() != Eigen::Success || solver_geodesics_.info() != Eigen::Success){
		std::cerr << "Unable to factorize the linear system matrices" << std::endl;
		return;
	}

	#pragma omp parallel for
	for(int i = 0; i < n_control_points; ++i){
		// Integrate the heat flow
		DenseVector u0 = DenseVector::Zero(n_sample_points);
		u0(control_point_idx_(i)) = 1.0;
		DenseVector u = solver_heatflow_.solve(u0);

		// Normalize the gradient field
        DenseVector normalized_grad_u = - cell_grad_operator_ * u;
		for(int j = 0; j < n_cells; ++ j){
			double grad_norm = normalized_grad_u.segment(dim * j, dim).norm();
			normalized_grad_u.segment(dim * j, dim) *= 1.0 / grad_norm;
		}

		// Solve the Poisson equation
		geodesic_distance_.block(0, i, n_sample_points-1, 1) = solver_geodesics_.solve(rhs_op * normalized_grad_u);

		// Shift the values to achieve zero at the control point
		geodesic_distance_.col(i) -= DenseVector::Constant(n_sample_points, geodesic_distance_(control_point_idx_(i), i));

        // Solve a linear system to compute harmonic coordinate values for inner sample points
        inner_point_init_values_.col(i) = solver_inner_init_vals_.solve(rhs_op_inner_vals * init_coordinate_values_.col(i));
	}
}


void DataSetup::compute_gradient_coefficients(const IndexVector &point_idx, DenseMatrix &grad_coef,
		DenseMatrix &integrated_grad_coef, double &cell_measure)
{
	int dim = point_idx.size() - 1;
	assert(dim == 2 || dim == 3);

	integrated_grad_coef.setZero(dim, point_idx.size());

	// Compute the cell measure using the determinant of a matrix that stores edge vectors
	DenseMatrix edge_mat(dim, dim);
	for(int i = 0; i < dim; ++ i){
		edge_mat.col(i) = sample_points_.col(point_idx(i+1)) - sample_points_.col(point_idx(0));
	}
	cell_measure = std::fabs(edge_mat.determinant()) / (dim == 2 ? 2.0 : 6.0);


	if(dim == 2){
		for(int i = 0; i < 3; ++ i){
			Vector2d edge_vec = sample_points_.col(point_idx((i+1)%3)) - sample_points_.col(point_idx((i+2)%3));
			Vector2d rotated_edge_vec;
			rotated_edge_vec(0) = -edge_vec(1);
			rotated_edge_vec(1) = edge_vec(0);

			// Check the direction of the gradient vector
			Vector2d check_vec = sample_points_.col(point_idx(i)) - sample_points_.col(point_idx((i+1)%3));
			if(check_vec.dot(rotated_edge_vec) < 0.0){
				rotated_edge_vec *= -1.0;
			}

			integrated_grad_coef.col(i) = rotated_edge_vec * 0.5;
		}
	}
	else{
		for(int i = 0; i < 4; ++ i){
			Vector3d edge_vec1 = sample_points_.col(point_idx((i+2)%4)) - sample_points_.col(point_idx((i+1)%4)),
					edge_vec2 = sample_points_.col(point_idx((i+3)%4)) - sample_points_.col(point_idx((i+1)%4));
			Vector3d integrated_grad_vec = edge_vec1.cross(edge_vec2) / 6.0;

			Vector3d check_vec = sample_points_.col(point_idx((i+1)%4)) - sample_points_.col(point_idx(i));
			if(check_vec.dot(integrated_grad_vec) < 0.0){
				integrated_grad_vec *= -1.0;
			}

			integrated_grad_coef.col(i) = integrated_grad_vec;
		}
	}

	grad_coef = integrated_grad_coef / cell_measure;
}

}


