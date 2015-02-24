//	OpenGL.cpp
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


#include "OpenGL.h"
#include <algorithm>
#include <iostream>

COpenGL::COpenGL(const QGLFormat &format, QWidget *parent)
    : QGLWidget(format, parent)
{ 
    lbc = new LBCUI;
    selected = -1;
}

COpenGL::~COpenGL()
{
}

void COpenGL::initializeGL()
{
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_DOUBLEBUFFER);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
}

void COpenGL::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    float xSpan = 1;
    float ySpan = 1;

    // Make sure the shorter side has length 1
    if(w > h){
    	xSpan = float(w) / h;
    }
    else{
    	ySpan = float(h) / w;
    }

    glOrtho(-xSpan, xSpan, -ySpan, ySpan, -1, 1);

    // Use the entire window for rendering.
    glViewport(0, 0, w, h);
}

void COpenGL::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glPushMatrix();

    Render();
    glPopMatrix();
}

void COpenGL::Render()
{
    lbc->draw(std::min(width(), height()));
}

void COpenGL::mousePressEvent(QMouseEvent *e)
{
    if (e->button() == Qt::LeftButton)
    {
        double x = e->pos().x();
        double y = e->pos().y();

        double base_length = std::min(width(), height());

        selected = lbc->select((2*x-width())/base_length, (-2*y+height())/base_length);

        updateGL();
    }
}
