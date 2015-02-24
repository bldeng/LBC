//	DataSetup.h
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

#ifndef DATASETUP_H_
#define DATASETUP_H_

#include "common.h"
#include <vector>
#include <cassert>

#ifdef USE_CHOLMOD
#include <Eigen/CholmodSupport>
#endif

namespace LBC{

class DataSetup{
public:
	enum WeightingScheme{
		CONSTANT = 0,
		LINEAR = 1,
		SQUARE = 2,
		SQUAREROOT = 3
	};

	// Storage for information about boundary facets of the control cage
	struct CageBoundaryFacetInfo{

		// Control point array indices (in the range [0, n_control_points - 1]) of the cage vertices on this boundary facet
		// In 2D, this vector is two-dimensional, while in 3D it is three-dimensional
		IndexVector facet_vertices;

		// Sample point array indices (in the range [0, n_sample_points - 1]) of boundary sample points
		// that are not control points and lie on this facet
		IndexVector boundary_points;

		// The barycentric coordinates of the boundary points w.r.t. the facet vertices.
		// Each row corresponds to the barycentric coordinates of one boundary point.
		// They are either computed from the point positions, or provided in the constructor.
		DenseMatrix barycentric_coordinates;

		// Whether or not the coordinates have been provided
		bool has_coordinates;

		CageBoundaryFacetInfo(const IndexVector &_facet_vertices, const IndexVector &_boundary_points)
			:facet_vertices(_facet_vertices), boundary_points(_boundary_points), has_coordinates(false){}

		CageBoundaryFacetInfo(const IndexVector &_facet_vertices, const IndexVector &_boundary_points, const DenseMatrix &_coordinates)
			:facet_vertices(_facet_vertices), boundary_points(_boundary_points),
			 barycentric_coordinates(_coordinates), has_coordinates(true)
		{
			assert(barycentric_coordinates.rows() == boundary_points.size() && barycentric_coordinates.cols() == facet_vertices.size());
		}

	};

	DataSetup(const DenseMatrix &sample_points, const IndexVector &control_point_idx, const DenseIndexMatrix &cell_vertices,
			const std::vector< CageBoundaryFacetInfo > &boundary_facet_info, WeightingScheme scheme = SQUARE):
			sample_points_(sample_points), control_point_idx_(control_point_idx), cell_vertices_(cell_vertices),
			boundary_facet_info_(boundary_facet_info)
	{
		init(scheme);
	}

    const DenseMatrix& get_LBC_solver_control_points() const
    {
        return control_points_;
    }

    const DenseMatrix& get_LBC_solver_data_points() const
    {
        return inner_points_;
    }

    const SparseMatrix& get_LBC_solver_grad_operator() const
    {
        return grad_operator_for_inner_points_;
    }

    const DenseMatrix& get_LBC_solver_grad_const() const
    {
        return grad_const_;
    }

    const DenseMatrix& get_LBC_solver_grad_weights() const
    {
        return grad_weights_;
    }

    const DenseMatrix& get_geodesic_distance() const
    {
        return geodesic_distance_;
    }

    const DenseMatrix& get_inner_point_init_values() const
    {
        return inner_point_init_values_;
    }

	DenseMatrix get_full_coordinate_values(const DenseMatrix &inner_point_values) const
	{
		return inner_points_mapping_ * inner_point_values  + init_coordinate_values_;
	}    

    // Compute the projectoin point of given points onto the affine subspace spanned by a set of base points,
    // represented using the barycentric coordinates w.r.t. the base points.
    // The base points are stored as columns of the matrix base_poitns
    // The computed barycentric coordiantes are returned in vector proj_bc
    static void compute_projection_barycentric_coordinates(const DenseMatrix& base_points, const DenseMatrix& pts, DenseMatrix& proj_bc);


private:

	// Positions for all sample points, including control points and data points.
	// Each sample point corresponds to one column of this matrix.
	DenseMatrix sample_points_;

	// The indices of control points among the sample points,
	// i.e., the indices of columns in matrix sample_points_ that correspond to the control points
	IndexVector control_point_idx_;

	// Index matrix showing which sample points belong to the same cell.
	// Each column collects the indices of sample points within one cell.
	// In 2D, a cell is a triangle, so this matrix has three rows.
	// In 3D, a cell is a tetrahedron, meaning that there are four rows in this matrix.
	DenseIndexMatrix cell_vertices_;

	std::vector< CageBoundaryFacetInfo > boundary_facet_info_;

	// Matrix representation of the gradient operator for each cell,
	// The number of columns equals the number of sample points,
	// while the number of rows equals D * N_c, where D is the problem dimension and N_c is the number of cells.
	// integrated_grad_operator is the weighted version of the operator, with the gradient weighted by the cell measure
	SparseMatrix cell_grad_operator_, integrated_grad_operator_;

	// A diagonal matrix for evaluating the measures associated with each vertex
	SparseMatrix vertex_associted_measure_;

	// Matrix representing the divergence at each vertex for a vector field defined on faces.
	// It is computed by evaluating the outward flux from the associated cell area for the vertex
	SparseMatrix integrated_divergence_operator_;

	SparseMatrix symmetric_laplacian_operator_;

	// Data required by the LBC solver
	DenseMatrix control_points_, inner_points_;
	SparseMatrix grad_operator_for_inner_points_;
	DenseMatrix grad_const_, grad_weights_;


	// Matrices for separating inner points and other points from the sample point matrix
	// These matrices can be right-multiplied to the sample point matrix to obtain the corresponding point positions
	SparseMatrix inner_points_mapping_, control_points_mapping_;

	// Initial barycentric coordinate values for the sample points.
	// For inner points, the coordinate values are zero.
	// For boundary points (i.e., points whose coordinate values are not subject to optimization),
	// the values are computed according to the Dirichlet boundary conditions
	DenseMatrix init_coordinate_values_;

    // Initial barycentric coordiante values for the inner sample points
    DenseMatrix inner_point_init_values_;

	// Geodesic distance from each control point to all sample points
	DenseMatrix geodesic_distance_;


	// Cholesky solver for the linear system
	#ifdef USE_CHOLMOD
    Eigen::CholmodDecomposition<SparseMatrix> solver_heatflow_, solver_geodesics_, solver_inner_init_vals_;
	#else
    Eigen::SimplicialLDLT<SparseMatrix> solver_heatflow_, solver_geodesics_, solver_inner_init_vals_;
	#endif

	double compute_geodesic_time_step();

    void compute_geodesic_distance_and_init_inner_points_coordinates();

	void compute_boundary_values();

	// Fill in the barycentric coordinates field of the CageBoundaryFacetInfo struct.
	void compute_barycentric_coordinates(CageBoundaryFacetInfo &info);

    // Normalize the sample point positions such that their centroid is at the origin, and the maximum distance from a sample point to the centroid is 1.
    void normalize_sample_points();

    // Compute the matrices for face gradient, integrated face gradient, integrated divergence, and the symmetric Laplacian operators
	void construct_operator_matrices();

	// For a given cell with vertex indices stored in point_idx, compute the gradient operator coefficients and
	// integrated gradient operator coefficients w.r.t. each vertex, as well as the measure (area or volume) of the cell.
	void compute_gradient_coefficients(const IndexVector &point_idx, DenseMatrix &grad_coef,
			DenseMatrix &integrated_grad_coef, double &cell_measure);

	void compute_solver_data(WeightingScheme scheme);

	void init(WeightingScheme scheme);
};

}

#endif /* DATASETUP_H_ */
