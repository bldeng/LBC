QT += core opengl gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = LBC
TEMPLATE = app

SOURCES += external/triangle/triangle.c \
        main.cpp \
        glwidget.cpp \
        OpenGL.cpp \
        DataSetup.cpp \
        LBCSolver.cpp \
        LBCUI.cpp
        

HEADERS += external/triangle/triangle.h \
        glwidget.h \
        OpenGL.h \
        common.h \
        DataSetup.h \
        LBCSolver.h \
        LBCUI.h

FORMS    += glwidget.ui

DEFINES += TRILIBRARY ANSI_DECLARATORS

# uncomment the following to use cholmod
# DEFINES += USE_CHOLMOD

INCLUDEPATH *= $${_PRO_FILE_PWD_}/external/eigen
INCLUDEPATH *= /usr/local/include/eigen3/
INCLUDEPATH *= /usr/include/eigen3/

DESTDIR = $${OUT_PWD}

unix:!macx{
    LIBS += -lgomp
    QMAKE_CXXFLAGS += -fopenmp
    QMAKE_POST_LINK = cp $${_PRO_FILE_PWD_}/colorbar_texture.png $${OUT_PWD}

    # Libraries for Cholmod
    contains(DEFINES, USE_CHOLMOD) {
        INCLUDEPATH *= /usr/include/suitesparse
        LIBS += -lamd -lcamd -lcolamd -lccolamd -lcholmod
    }
}

win32{
    DEFINES += NO_TIMER NOMINMAX
    QMAKE_CXXFLAGS *= /openmp
    QMAKE_CXXFLAGS *= /MP

    # Copy texture file to output
    VIZ_TEXTURE += $${_PRO_FILE_PWD_}/colorbar_texture.png
    VIZ_TEXTURE ~= s,/,\\,g
    QMAKE_POST_LINK += $$quote(cmd /c copy /y $${VIZ_TEXTURE} $${DESTDIR})
}
