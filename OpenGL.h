//	OpenGL.h
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


#ifndef OPENGL_H
#define OPENGL_H

#include <QGLFormat>
#include <QGLWidget>
#include <QMouseEvent>
#include <QEvent>
#include "LBCUI.h"


class COpenGL : public QGLWidget
{
public:
    COpenGL(const QGLFormat &format, QWidget *parent = 0);
    ~COpenGL();

protected:

    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
    void mousePressEvent(QMouseEvent *e);

    int selected;

public:

    void Render();

public:    
    LBCUI *lbc;
};

#endif
