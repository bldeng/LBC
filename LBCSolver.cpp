//	LBCSolver.cpp
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

#include "LBCSolver.h"
#include <vector>
#include <iostream>
#include <algorithm>

namespace LBC{

void LBCSolver::update_x()
{
	#pragma omp for
	for(int i = 0; i < n_grad_elements_; ++i){
		int f_index = i % n_cells_;
		int r_index = dim_ * f_index;
		int c_index = i / n_cells_;
        double a = param_.penalty_weight * X_.block(r_index, c_index, dim_, 1).norm();
        double current_weight = grad_weights_(f_index, c_index);

		if (a <= current_weight){
            X_.block(r_index, c_index, dim_, 1).setZero();
		}
		else{
            X_.block(r_index, c_index, dim_, 1) *= (1.0 - current_weight/a);
		}
	}
}

void LBCSolver::update_y()
{
	#pragma omp for
	for (int i = 0; i < n_control_points_; ++i){
		X_h_.col(i) = X_.col(i) * param_.relaxation_alpha - Grad_.col(i) * (param_.relaxation_alpha - 1);
		W_h_.col(i) = W_.col(i) * param_.relaxation_alpha - Z_.col(i) * (param_.relaxation_alpha - 1);
		pre_rhs_.col(i) = Gt_ * (X_h_.col(i) - J_.col(i) + d1_.col(i)) + W_h_.col(i) - H_.col(i) + d2_.col(i);
	}

    #pragma omp for
    for (int i = 0; i < n_Y_col_; ++i){
        Y_.col(i) = solver_.solve(pre_rhs_ * Mt_.col(i));
    }
}


void LBCSolver::update_w()
{
    #pragma omp for
    for (int i = 0; i < n_total_coef_; ++i)
    {
        int v_index = i % n_data_points_;
        int c_index = i / n_data_points_;
        W_(v_index, c_index) = std::min(1.0, std::max(0.0, W_(v_index, c_index)));
    }
}

void LBCSolver::update_dual_variables(int iter_num)
{
	#pragma omp sections
	{
		#pragma omp section
		{
			if(check_convergence_){
				prev_Z_ = Z_;
			}
		}

		#pragma omp section
		{
			if(check_convergence_){
				prev_Grad_ = Grad_;
			}
		}

        #pragma omp section
        {
            primal_residual_X_ = X_;
        }

        #pragma omp section
        {
            primal_residual_W_ = W_;
        }
	}

	#pragma omp for
	for (int i = 0; i < n_control_points_; i++)
	{
		Z_.col(i) = Y_ * M_.col(i) + H_.col(i);
	}

	#pragma omp for
	for (int i = 0; i < n_control_points_; i++)
	{
		Grad_.col(i) = G_ * Z_.col(i) + E_.col(i);
	}


    #pragma omp sections
    {
        #pragma omp section
        {
        	if(check_convergence_){
        		dual_residual_X_sqr_norm_ = (Grad_ - prev_Grad_).squaredNorm();
        	}
        }

		#pragma omp section
        {
        	if(check_convergence_){
        		dual_residual_W_sqr_norm_ = (Z_ - prev_Z_).squaredNorm();
        	}
        }

        #pragma omp section
        {
        	if(check_convergence_){
                primal_residual_X_ -= Grad_;
        		primal_residual_X_sqr_norm_ = primal_residual_X_.squaredNorm();
        	}
        }

        #pragma omp section
        {
        	if(check_convergence_){
                primal_residual_W_ -= Z_;
        		primal_residual_W_sqr_norm_ = primal_residual_W_.squaredNorm();
        	}
        }

		#pragma omp section
		{
			d1_ += ( X_h_ - Grad_ );
			X_ = Grad_ - d1_;
		}

		#pragma omp section
		{
			d2_ += ( W_h_ - Z_ );
			W_ = Z_ - d2_;
		}
    }

    #pragma omp single
    {
    	if(check_convergence_){
			primal_residual_sqr_norm_ = primal_residual_X_sqr_norm_ + primal_residual_W_sqr_norm_;
			dual_residual_sqr_norm_ = (dual_residual_X_sqr_norm_ + dual_residual_W_sqr_norm_)
					* param_.penalty_weight * param_.penalty_weight;

			optimization_converge_ = ( primal_residual_sqr_norm_ <= primal_residual_sqr_norm_threshold_
					&& dual_residual_sqr_norm_ <= dual_residual_sqr_norm_threshold_ );

			optimization_end_ = optimization_converge_ || iter_num >= param_.max_iterations;

			if(optimization_converge_){
				std::cout << "Solver converged." << std::endl;
			}
			else if(optimization_end_){
				std::cout << "Maximum iteration reached." << std::endl;
			}

			if (output_progress_ || optimization_end_)
			{
				std::cout << "Iteration " << iter_num << ":" << std::endl;
				std::cout << "Primal residual squared norm: " << primal_residual_sqr_norm_ << ",  primal threshold:" << primal_residual_sqr_norm_threshold_ << std::endl;
				std::cout << "Dual residual squared norm: " << dual_residual_sqr_norm_ << ",  dual threshold:" << dual_residual_sqr_norm_threshold_ << std::endl;
				std::cout << std::endl;
			}
    	}
    }
}

void LBCSolver::initialize_thresholds()
{
	double primal_threshold = std::max(param_.abs_eps, dim_ * n_cells_ * n_control_points_ * param_.rel_primal_eps);
	double dual_threshold = std::max(param_.abs_eps, n_data_points_ * n_control_points_ * param_.rel_dual_eps);
	primal_residual_sqr_norm_threshold_ = primal_threshold * primal_threshold;
	dual_residual_sqr_norm_threshold_ = dual_threshold * dual_threshold;
}

void LBCSolver::initialize_linear_constraint_data()
{
    DenseMatrix Kt, Bt;
    Kt.resize(dim_ + 1, n_control_points_);
    Bt.resize(dim_ + 1, n_data_points_);

    Kt.block(0, 0, dim_, n_control_points_) = control_points_;
    Kt.row(dim_).fill(1.0);
    Bt.block(0, 0, dim_, n_data_points_) = data_points_;
    Bt.row(dim_).fill(1.0);

    Eigen::JacobiSVD< DenseMatrix, Eigen::FullPivHouseholderQRPreconditioner > jsvd(Kt, Eigen::ComputeFullU | Eigen::ComputeFullV);
    int nrank = jsvd.nonzeroSingularValues();
    Mt_ = jsvd.matrixV().block(0, nrank, n_control_points_, n_control_points_ - nrank);
    M_ = Mt_.transpose();
    n_Y_col_ = M_.rows();
    H_ = jsvd.solve(Bt).transpose();	// Least squares solving

}

void LBCSolver::initialize_solver()
{
    //construct identity matrix for data point size
    SparseMatrix Id(n_data_points_, n_data_points_);

    std::vector<TripletD> triplets_I;
    for (int i = 0; i < n_data_points_; i++)
    {
        triplets_I.push_back(TripletD(i, i, 1.0));
    }
    Id.setFromTriplets(triplets_I.begin(), triplets_I.end());

    SparseMatrix M = Gt_*G_ + Id;
    solver_.compute(M);
    if(solver_.info() != Eigen::Success){
    	std::cerr << "Linear system matrix factorization failed" << std::endl;
    	valid_init_data_ = false;
    }
}


void LBCSolver::initialize_coordinates(const DenseMatrix &coord)
{
	init_coord_ = coord;
	if(init_coord_.rows() != n_data_points_ || init_coord_.cols() != n_control_points_){
		std::cerr << "Invalid size of initial coordinates" << std::endl;
		valid_init_data_ = false;
	}

	has_init_coord_ = true;
}


void LBCSolver::init()
{
	valid_init_data_ = true;

    dim_ = control_points_.rows();
    if(data_points_.rows() != dim_){
    	std::cerr << "Dimension mismatch between control points and data points" << std::endl;
    	valid_init_data_ = false;
    }

    n_control_points_ = control_points_.cols();
    n_data_points_ = data_points_.cols();
    n_total_coef_ = n_control_points_ * n_data_points_;
    if(n_control_points_ == 0 || n_data_points_ == 0){
    	std::cerr << "Invalid number of control points or data poins" << std::endl;
    	valid_init_data_ = false;
    }

    Gt_ = G_.transpose();
    n_cells_ = G_.rows() / dim_;
    n_grad_elements_ = n_cells_ * n_control_points_;
    if (G_.rows() % dim_){
    	std::cerr << "Invalid dimension for the gradient operator" << std::endl;
    	valid_init_data_ = false;
    }
    if(n_cells_ == 0){
    	std::cerr << "Invalid number of cells" << std::endl;
    	valid_init_data_ = false;
    }

    if(grad_weights_.rows() != n_cells_ || grad_weights_.cols() != n_control_points_){
    	std::cerr << "Invalid dimension of gradient weights" << std::endl;
    	valid_init_data_ = false;
    }

    convergence_check_frequency_ = std::max(1, param_.convergence_check_frequency);
    output_frequency_ = std::max(1, param_.output_frequency_ratio * convergence_check_frequency_);

    initialize_solver();
    initialize_linear_constraint_data();
    initialize_thresholds();
}

void LBCSolver::initialize_variables()
{
	if(has_init_coord_){
        Z_ = init_coord_;
	}
	else{
        Z_.resize(n_data_points_, n_control_points_);
        Z_.fill(1.0 / n_control_points_);
	}

    Y_ = (Z_ - H_) * Mt_;

    W_ = Z_;
    Grad_ = G_ * W_ + E_;
    J_ = G_ * H_ + E_;
    X_ = Grad_;
    X_h_ = X_;
    W_h_ = W_;
    primal_residual_W_.setZero(W_.rows(), W_.cols());
    primal_residual_X_.setZero(X_.rows(), X_.cols());
    pre_rhs_.setZero(H_.rows(), H_.cols());
    d1_.setZero(X_.rows(), X_.cols());
    d2_.setZero(W_.rows(), W_.cols());

}

void LBCSolver::start_timer(){
	if(param_.use_timer){
		start_time_ = omp_get_wtime();
	}
}

void LBCSolver::end_timer(){
	if(param_.use_timer){
		end_time_ = omp_get_wtime();
	}
}

void LBCSolver::show_elapsed_time(){
	if(param_.use_timer){
		std::cout << "Solving time: " << std::max(0.0, end_time_ - start_time_) << " seconds." << std::endl;
	}
}

void LBCSolver::solve()
{
	if(!valid_init_data_){
		std::cerr << "Invalid data, unable to solve...." << std::endl;
		return;
	}

	initialize_variables();

    start_timer();
    int iter = 0;
    optimization_end_ = false;

    while (!optimization_end_)
    {
        iter++;
        check_convergence_ = (iter % convergence_check_frequency_ == 0 || iter >= param_.max_iterations);
        output_progress_ = (iter % output_frequency_ == 0 );

        #pragma omp parallel
        {
            this->update_x();

            this->update_w();

            this->update_y();

            this->update_dual_variables(iter);
        }
    }

    end_timer();
    show_elapsed_time();
}

}


