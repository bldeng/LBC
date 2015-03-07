//	LBCUI.cpp
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

#include "LBCSolver.h"
#include "LBCUI.h"
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>


#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#ifndef VOID
#define VOID void
#endif

extern "C"{
#include "external/triangle/triangle.h"
}


bool LBCUI::read_cage(const std::string& filename)
{
    std::ifstream inCurves;
    inCurves.open(filename.c_str());
    if (!inCurves.is_open())
    {
        std::cerr << "Unable to open file" << std::endl;
        return false;
    }

    // We assume the cage vertices are stored consecutively in the cage obj file
    std::vector<double> coordinates;
    while (inCurves.good())
    {
    	std::string line;
        std::getline(inCurves, line);

        //this line is vertex coordinate
        if ((!inCurves.fail()) && (!line.substr(0, 2).compare("v ")))
        {
            std::stringstream ss(line.substr(2));
            double x, y, z;
            ss >> x >> y >> z;

            if(!ss.fail()){
				coordinates.push_back(x);
				coordinates.push_back(y);
				coordinates.push_back(z);
            }
        }
    }

    m_has_cage = coordinates.size() > 6;
    if(!m_has_cage){
    	std::cerr << "Invalid cage file" << std::endl;
    	return false;
    }

    int n_cage_vtx = coordinates.size()/3;
    m_cage = Eigen::Map<Eigen::Matrix3Xd>(&(coordinates[0]), 3, n_cage_vtx);
    m_cage.row(2).fill(0.0);

    Eigen::Vector3d min_coord = m_cage.rowwise().minCoeff(), max_coord = m_cage.rowwise().maxCoeff();
    m_center = (min_coord + max_coord) * 0.5;
    m_scale = (max_coord - min_coord).norm() * 0.5;

    m_avg_cage_edge_length = 0.0;
    for(int i = 0; i < n_cage_vtx; ++ i){
    	m_avg_cage_edge_length += (m_cage.col(i) - m_cage.col((i+1)% n_cage_vtx)).norm();
    }
    m_avg_cage_edge_length /= n_cage_vtx;

    m_cage_display = (m_cage.colwise() - m_center) / m_scale;

    return true;
}

void LBCUI::setup_texture()
{
    if(m_texture_initialized){
        return;
    }

    m_texture.load("./colorbar_texture.png");
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &m_TextureID);
    glBindTexture(GL_TEXTURE_2D, m_TextureID);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    glTexImage2D(GL_TEXTURE_2D, 0, m_texture.depth()/8, m_texture.width(),
         m_texture.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, m_texture.bits());

    m_texture_initialized = true;
}

void LBCUI::set_mouse_select()
{
    m_mouse_picking = !m_mouse_picking;
}

int LBCUI::select(double x, double y)
{
    if (m_mouse_picking && m_cage.cols() > 0)
    {
    	Eigen::Vector3d pt(x, y, 0.0);
    	pt *= m_scale;
    	pt += m_center;

    	Eigen::Matrix3Xd diff = m_cage;
    	diff.colwise() -= pt;
    	int min_idx = -1;
    	double min_dist = diff.colwise().norm().minCoeff(&min_idx);

    	if(min_dist < m_avg_cage_edge_length * 0.2){
    		m_selected_point = min_idx;
    		return m_selected_point;
    	}
    }

    return -1;
}

void LBCUI::delaunay_Triangulation(int numTriangles)
{
    double domainArea = domain_area();
    double averageArea = domainArea/numTriangles*1.6;
    double edge_length = sqrt(2.0/sqrt(3)*averageArea)*1.6;

    std::vector<double> position_x;
    std::vector<double> position_y;

    m_edge_length = edge_length;
    std::cout << averageArea << " " << edge_length << std::endl;
    std::cout << "Done" << std::endl;
    int n_boundary_vertex = 0;

    int n_control_pts = m_cage.cols();
    for (int i=0; i< n_control_pts; ++i)
    {
    	int next_i = (i + 1) % n_control_pts;
        int n_segment = std::ceil((m_cage.col(i) - m_cage.col(next_i)).norm() / edge_length);
        n_boundary_vertex += n_segment;

        for (int j=0; j<n_segment; j++)
        {
            Eigen::Vector3d boundary_pt = m_cage.col(i) + j * (m_cage.col(next_i) - m_cage.col(i)) / n_segment;
            position_x.push_back(boundary_pt(0));
            position_y.push_back(boundary_pt(1));
        }
    }


    struct triangulateio in, out;

    /* Define input points. */

    in.numberofpoints = n_boundary_vertex;
    in.numberofpointattributes = 0;
    in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
    for (int i=0; i<n_boundary_vertex; ++i)
    {
        in.pointlist[2*i] = position_x[i];
        in.pointlist[2*i+1] = position_y[i];
    }
    in.numberofsegments = n_boundary_vertex;
    in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(REAL));
    for (int i=0; i< n_boundary_vertex-1; i++)
    {
        in.segmentlist[2*i] = i;
        in.segmentlist[2*i+1] = i+1;
    }
    in.segmentlist[2*n_boundary_vertex-2] = n_boundary_vertex-1;
    in.segmentlist[2*n_boundary_vertex-1] = 0;

    in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
    for (int i=0; i<n_boundary_vertex; i++)
    {
        in.pointmarkerlist[i] = 1;
    }
    in.segmentmarkerlist = (int *) malloc(in.numberofsegments * sizeof(int));
    for (int i=0; i<n_boundary_vertex; i++)
    {
        in.segmentmarkerlist[i] = 0;
    }
    in.numberofholes = 0;
    in.numberofregions = 0;

    out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
    /* Not needed if -N switch used or number of point attributes is zero: */
    out.pointattributelist = (REAL *) NULL;
    out.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
    out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
    /* Not needed if -E switch used or number of triangle attributes is zero: */
    out.triangleattributelist = (REAL *) NULL;
    /* Needed only if segments are output (-p or -c) and -P not used: */
    out.segmentlist = (int *) NULL;
    /* Needed only if segments are output (-p or -c) and -P and -B not used: */
    out.segmentmarkerlist = (int *) NULL;

    std::ostringstream ostr_switches;
    ostr_switches << "pza" << averageArea << "q33.0Y";
    std::string switches = ostr_switches.str();
    std::vector<char> switch_chars(switches.length() + 1, '\0');
    switches.copy(&(switch_chars[0]), switches.length());
    triangulate(&(switch_chars[0]), &in, &out, (struct triangulateio *) NULL);

    m_delaunay.resize(3, out.numberofpoints);
    for(int i = 0; i < out.numberofpoints; ++i){
    	m_delaunay.col(i) = Eigen::Vector3d(out.pointlist[2*i], out.pointlist[2*i+1], 0.0);
    }

    m_delaunay_faces = Eigen::Map<Eigen::MatrixXi>(out.trianglelist, 3, out.numberoftriangles);
    m_delaunay_display = (m_delaunay.colwise() - m_center) / m_scale;


    /* Free all allocated arrays, including those allocated by Triangle. */
    free(in.pointlist);
    free(out.pointlist);
    free(out.pointattributelist);
    free(out.trianglelist);
    free(out.triangleattributelist);

    m_has_delaunay = true;
}

double LBCUI::domain_area()
{
	double signed_area = 0.0;

	int n = m_cage.cols();
	for(int i = 0; i < n - 1; ++i){
		signed_area += m_cage.block(0, i, 2, 2).determinant();
	}

	Eigen::Matrix2d M_;
	M_.col(0) = m_cage.block(0, n-1, 2, 1);
	M_.col(1) = m_cage.block(0, 0, 2, 1);
	signed_area += M_.determinant();

	return 0.5 * std::abs(signed_area);
}

void LBCUI::lbc_solver(int weighting_scheme, int max_iter, double relaxation_alpha, double penalty_weight,
                       double primal_residual_tolerance, double dual_residual_tolerance)
{
	if(!m_has_cage){
		std::cerr << "Error: no cage data" << std::endl;
		return;
	}

	if(!m_has_delaunay){
		std::cerr << "Error: no triangulation data" << std::endl;
		return;
	}

	std::cout << "Setting up input data" << std::endl;

    LBC::DenseMatrix sample_points = m_delaunay.block(0, 0, 2, m_delaunay.cols());
    LBC::DenseIndexMatrix &cell_vertices = m_delaunay_faces;

    // control_point_idx
    LBC::IndexVector control_point_idx(m_cage.cols());
    for (int i = 0; i < m_cage.cols(); ++i)
    {
        Eigen::Vector2d ctrl_pt = m_cage.block(0, i, 2, 1);
        LBC::DenseMatrix candidate_points = sample_points;
        int min_idx = -1;
        (candidate_points.colwise() - ctrl_pt).colwise().norm().minCoeff(&min_idx);
        control_point_idx[i] = min_idx;
    }

    // boundary_facet_info
    std::vector< LBC::DataSetup::CageBoundaryFacetInfo > boundary_facet_info;
    int n_control_pt = m_cage.cols();
    for (int i = 0; i < n_control_pt; ++i)
    {
        int next_i = (i+1) % n_control_pt;

        LBC::IndexVector facet_vertices(2);
        facet_vertices(0) = i;
        facet_vertices(1) = next_i;
        int n_segments = std::ceil((m_cage.col(i) - m_cage.col(next_i)).norm() / m_edge_length);

        LBC::IndexVector boundary_points(n_segments-1);
        for (int k = 1; k < n_segments; ++k)
        {
            Eigen::Vector3d pt_3d = m_cage.col(i) + k * (m_cage.col(next_i) - m_cage.col(i)) / n_segments;
            Eigen::Vector2d pt = pt_3d.head(2);
            LBC::DenseMatrix candidate_points = sample_points;
            int min_idx = -1;
            (candidate_points.colwise() - pt).colwise().norm().minCoeff(&min_idx);
            boundary_points[k-1] = min_idx;
        }

        LBC::DataSetup::CageBoundaryFacetInfo edge(facet_vertices, boundary_points);
        boundary_facet_info.push_back(edge);
    }

    LBC::DataSetup::WeightingScheme scheme = static_cast<LBC::DataSetup::WeightingScheme>(weighting_scheme);
    LBC::DataSetup ds(sample_points, control_point_idx, cell_vertices, boundary_facet_info, scheme);


    LBC::Param param;
    param.max_iterations = max_iter;
    param.relaxation_alpha = relaxation_alpha;
    param.convergence_check_frequency = 10;
    param.output_frequency_ratio = 10;
    param.rel_primal_eps = primal_residual_tolerance;
    param.rel_dual_eps = dual_residual_tolerance;
    param.penalty_weight = penalty_weight;
    LBC::LBCSolver solver(param, ds);

    std::cout << "LBC Solver started" << std::endl;

    solver.solve();

    std::cout << "Finished computation" << std::endl;


    /*
    m_w = ds.get_geodesic_distance();
    for(int i = 0; i < m_w.cols(); ++ i){
        std::cout << "Column " << i << std::endl;
        std::cout << "Max distance: " << m_w.col(i).maxCoeff() <<", Min distance: " << m_w.col(i).minCoeff() << std::endl;
        m_w.col(i) /= m_w.col(i).maxCoeff();
    }
    */

    m_w = ds.get_full_coordinate_values(solver.get_coordinates());
    setup_texture();
}


void LBCUI::draw(int pixel_scale)
{
	if(!m_has_cage){
		return;
	}

	double deafult_width = 3.0;
	double triangle_line_width = std::min(deafult_width, pixel_scale * 0.04 * m_edge_length / m_scale);

	glLineWidth(triangle_line_width);

    if (m_has_delaunay && m_mouse_picking && m_w.rows() == m_delaunay.cols() && m_w.cols() == m_cage.cols()
    		&& m_selected_point >= 0 && m_selected_point < static_cast<int>(m_cage.cols()))
    {
    	draw_delaunay_texture();
    }
    else{
    	draw_delaunay();
    }

    int n_control_pts = m_cage.cols();

    glColor3f(0.5, 0.5, 0.5);
    glLineWidth(deafult_width);
    glBegin(GL_LINES);

    for(int i = 0; i < n_control_pts; ++ i){
    	int next_i = (i+1) % n_control_pts;
    	glVertex2d(m_cage_display(0, i), m_cage_display(1, i));
    	glVertex2d(m_cage_display(0, next_i), m_cage_display(1, next_i));
    }

    glEnd();

    glPointSize(8.0);
    glBegin(GL_POINTS);
    for(int i = 0; i < n_control_pts; ++ i)
    {
        if (m_selected_point == i)
        {
            glColor3f(1.0f, 0.0f, 0.0f);
        }
        else
        {
            glColor3f(0.8f, 0.8f, 0.0f);
        }

        glVertex2d(m_cage_display(0, i), m_cage_display(1, i));
    }

    glEnd();
}

void LBCUI::draw_delaunay()
{
	glPushAttrib(GL_ALL_ATTRIB_BITS);
    glColor3f(0.0, 0.0, 1.0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_TRIANGLES);

    int n_cells = m_delaunay_faces.cols();
    for(int i = 0; i < n_cells; ++ i){
    	for(int k = 0; k < 3; ++ k){
    		int col_idx = m_delaunay_faces(k, i);
    		glVertex2d(m_delaunay_display(0, col_idx), m_delaunay_display(1, col_idx));
    	}
    }

    glEnd();
    glPopAttrib();
}

void LBCUI::draw_delaunay_texture()
{
	glPushAttrib(GL_ALL_ATTRIB_BITS);
    glEnable( GL_TEXTURE_2D );
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glBegin(GL_TRIANGLES);


    int n_cells = m_delaunay_faces.cols();
    for(int i = 0; i < n_cells; ++ i){
    	for(int k = 0; k < 3; ++ k){
    		int pt_idx = m_delaunay_faces(k, i);
    		double x = map_value_to_tex_coord( m_w(pt_idx, m_selected_point) );

            GLfloat tex[2];
            tex[0] = x;
            tex[1] = x;
            glTexCoord2fv(tex);
    		glVertex2d(m_delaunay_display(0, pt_idx), m_delaunay_display(1, pt_idx));
    	}
    }

    glEnd();
    glDisable( GL_TEXTURE_2D );
    glPopAttrib();
}

double LBCUI::map_value_to_tex_coord(double x, bool log_scale)
{
	if(log_scale){
        if (x < 1.e-5){
            x = 0.0;
        }
        else{
            x = std::min(1.0, std::max(0.0, 1.0 + log10(x)*0.2));
        }

        return x;
	}
	else{
		return std::min(1.0, std::max(0.0, x));
	}
}
