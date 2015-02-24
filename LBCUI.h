//	LBCUI.h
//
//  Copyright (C) 2015 Zishun Liu <liuzishun@gmail.com>,
//                     Bailin Deng <bldeng@gmail.com>
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


#ifndef LBCUI_H_
#define LBCUI_H_

#include "common.h"
#include <QImage>
#include <QtOpenGL>


class LBCUI
{
public:
    LBCUI()
	:m_edge_length(1.0), m_scale(1.0), m_avg_cage_edge_length(1.0), m_mouse_picking(false), m_has_cage(false),
    m_has_delaunay(false), m_selected_point(-1), m_texture_initialized(false){}

    ~LBCUI();

public:

    bool read_cage(const std::string& filename);

    void delaunay_Triangulation(int numTriangles);    

    void draw(int pixel_scale);

    void setup_texture();

    void set_mouse_select();

    int select(double x, double y);

    void lbc_solver(int weighting_scheme, int max_iter, double relaxation_alpha, double penalty_weight,
                    double primal_residual_tolerance, double dual_residual_tolerance);

private:
    double domain_area();

    void draw_delaunay();

    void draw_mesh();

    void draw_delaunay_texture();

    void draw_mesh_texture();

    double map_value_to_tex_coord(double value, bool log_scale = true);

private:
    Eigen::Matrix3Xd m_delaunay, m_delaunay_display;
    Eigen::MatrixXi m_delaunay_faces;
    Eigen::Vector3d m_center;
    Eigen::Matrix3Xd m_cage, m_cage_display;
    Eigen::MatrixXd m_w;

    double m_edge_length;
    double m_scale;
    double m_avg_cage_edge_length;

    bool m_mouse_picking;
    bool m_has_cage;
    bool m_has_delaunay;
    int m_selected_point;

    QImage m_texture;    
    GLuint m_TextureID;
    bool m_texture_initialized;
};

#endif
