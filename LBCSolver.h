//	LBCSolver.h
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

#ifndef LBCSOLVER_H_
#define LBCSOLVER_H_


#include "common.h"
#include "DataSetup.h"
#include <omp.h>
#include <algorithm>



namespace LBC{

	struct Param{
		double penalty_weight;

		// The relaxation coefficient for ADMM
		double relaxation_alpha;

		// Relative convergence threshold for primal and dual residuals
		double rel_primal_eps, rel_dual_eps;

		// Absolute threshold for primal and dual residuals
		double abs_eps;

		int max_iterations;

		// How often should we check the convergence of the solver?s
		int convergence_check_frequency;

		// How often should we output the progress of the solver?
		// This variable represents the ratio between the output frequency and the convergence check frequency
		int output_frequency_ratio;


		bool use_timer;



		// Constructor that sets the default values
		Param():
        penalty_weight(10), relaxation_alpha(1.0), rel_primal_eps(1e-6), rel_dual_eps(1e-6),
		abs_eps(1e-8), max_iterations(10000), convergence_check_frequency(10),
		output_frequency_ratio(10),
		use_timer(true)
		{}
	};

	class LBCSolver{
	public:
        LBCSolver(const Param &param, const DataSetup &data_setup):
                param_(param), control_points_(data_setup.get_LBC_solver_control_points()),
				data_points_(data_setup.get_LBC_solver_data_points()),
                grad_weights_(data_setup.get_LBC_solver_grad_weights()),
				G_(data_setup.get_LBC_solver_grad_operator()), E_(data_setup.get_LBC_solver_grad_const()),
                has_init_coord_(false), start_time_(0.0), end_time_(0.0)
		{
			init();
		}

		virtual ~LBCSolver(){}

		void initialize_coordinates(const DenseMatrix &coord);

		const DenseMatrix& get_coordinates() const{
			return Z_;
		}

		void solve();

	protected:

		// Solver parameters
		Param param_;

		// Positions of control points and sample data points.
		// The number of rows is the same as the problem dimension.
		DenseMatrix control_points_, data_points_;

		// Gradient weights on the cells, of size N_f x N_c, where N_f is the number of cells,
		// N_c is the number of control points.
		DenseMatrix grad_weights_;

		// G: Gradient operator for cells, with size (d N_f) x N_v, where N_f is the number of cells,
		// N_v is the number of non-boundary vertices, and d is the problem dimension
		// Gt is the transpose of G
		SparseMatrix G_, Gt_;

		// Constant factor in the face gradient computation
		DenseMatrix E_;

		// Other variables used in the solver (see the LBC paper for their definitions)
		DenseMatrix W_, X_, Y_, M_, Mt_, H_, J_;

		// Relaxation variables for W and X
		DenseMatrix W_h_, X_h_;

		// Scaled dual variables
		DenseMatrix d1_, d2_;

		// The coordinate values computed from Y: Z = Y * M + H
		DenseMatrix Z_, prev_Z_;

		// Temporary storage for face gradient values
		DenseMatrix Grad_, prev_Grad_;

		// Initial values for the coordinate values
		DenseMatrix init_coord_;
		bool has_init_coord_;

		// Temporary storage for computing the rhs of the linear system for Y
		DenseMatrix pre_rhs_;

		// Boolean variables for the progress of the solver
		bool optimization_converge_, optimization_end_;

		// Timer variables
		double start_time_, end_time_;

		// Are the initial data valid?
		bool valid_init_data_;

		// Problem dimension
		int dim_;

		int n_control_points_;
		int n_data_points_;
		int n_cells_;
		int n_total_coef_;
		int n_grad_elements_;
		int n_Y_col_;

		bool check_convergence_;
		bool output_progress_;
		int output_frequency_, convergence_check_frequency_;

		// Variables for primal and dual residuals
		DenseMatrix primal_residual_X_, primal_residual_W_;
		double primal_residual_X_sqr_norm_, primal_residual_W_sqr_norm_;
		double dual_residual_X_sqr_norm_, dual_residual_W_sqr_norm_;
		double primal_residual_sqr_norm_, dual_residual_sqr_norm_;
		double primal_residual_sqr_norm_threshold_, dual_residual_sqr_norm_threshold_;


		// Cholesky solver for the linear system
		#ifdef USE_CHOLMOD
		Eigen::CholmodDecomposition<SparseMatrix> solver_;
		#else
		Eigen::SimplicialLDLT<SparseMatrix> solver_;
		#endif


		// Update steps for variables W, X, Y
		void update_w();
		void update_x();
		void update_y();

		// initialize primal and dual thresholds for solver convergence
		void initialize_thresholds();

		// Pre-process the linear equality constraints, to compute matrices M and H
		void initialize_linear_constraint_data();

		// Pre-factorize the linear system matrix for Y
		void initialize_solver();

		// Initialization for variables used in the solver
		void initialize_variables();

		void init();

		// Update the dual variables and check if the solver converges
		void update_dual_variables(int iter_num);


		// Timer methods
		void start_timer();
		void end_timer();
		void show_elapsed_time();
	};

}

#endif /* LBCSOLVER_H_ */
