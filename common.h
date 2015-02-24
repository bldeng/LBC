//	common.h
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


#ifndef COMMON_H_
#define COMMON_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace LBC{

typedef Eigen::MatrixXd DenseMatrix;
typedef Eigen::Matrix3Xd DenseMatrix3X;
typedef Eigen::VectorXd DenseVector;
typedef Eigen::MatrixXi DenseIndexMatrix;
typedef Eigen::VectorXi IndexVector;
typedef Eigen::Vector2d Vector2d;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::Triplet<double> TripletD;
typedef Eigen::SparseMatrix<double> SparseMatrix;

}

#endif /* COMMON_H_ */
