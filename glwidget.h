//	glwidget.h
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

#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QtGui>
#include "OpenGL.h"
#include "ui_glwidget.h"

class GLWIDGET : public QMainWindow
{
    Q_OBJECT

public:
    GLWIDGET(QWidget *parent = 0);
    ~GLWIDGET();

private:    
    QCheckBox   *checkBox_DisplayWeights_;
    QPushButton *pushButton_Cage_;
    QLabel		*label_numTriangles_;
    QLabel		*label_weighting_;
    QComboBox	*combobox_weighting_;
    QSpinBox    *spinBox_numTriangles_;
    QPushButton *pushButton_triangulation_;
    QPushButton *pushButton_LBC_solver_;
    QGroupBox   *groupbox_triangulation_;
    QGroupBox   *groupbox_solver_;


    QSpinBox    *spinBox_max_iter_;
    QDoubleSpinBox    *dSpinBox_relaxation_alpha_;
    QDoubleSpinBox    *dSpinBox_penalty_;
    QDoubleSpinBox    *dSpinBox_primal_res_tolerance_;
    QDoubleSpinBox    *dSpinBox_dual_res_tolerance_;

    COpenGL     *gl_Widget_;

private slots:
    void load_cage();
    void triangulation();
    void select();
    void lbc_solver();

private:
    Ui::GLWIDGETClass ui;
    bool triangulation_ready_, solver_ready_, display_ready_;

    void update_controls();
};

#endif // GLWIDGET_H
