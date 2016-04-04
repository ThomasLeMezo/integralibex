#-------------------------------------------------
#
# Project created by QtCreator 2014-11-14T09:49:56
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = integrationIBEX
TEMPLATE = app

INCLUDEPATH += /usr/local/include/ibex /opt/VIBES/client-api/C++/src/

QMAKE_CXXFLAGS +=-I/usr/local/include/ibex -O0 -DNDEBUG -Wno-deprecated -frounding-math -std=c++11 -fopenmp
QMAKE_LFLAGS +=  -fopenmp

LIBS += -L/usr/local/include/ibex -libex -lprim

SOURCES += main.cpp \
    border.cpp \
    pave.cpp \
    scheduler.cpp \
    utils.cpp \
    tests.cpp \
    /opt/VIBES/client-api/C++/src/vibes.cpp \
    graph.cpp \
    inclusion.cpp \
    graphdot.cpp \
    imageintegral.cpp

HEADERS  +=     border.h \
    pave.h \
    scheduler.h \
    utils.h \
    tests.h \
    /opt/VIBES/client-api/C++/src/vibes.h \
    graph.h \
    inclusion.h \
    graphdot.h \
    imageintegral.h

INCLUDEPATH += /usr/local/include/opencv2
LIBS += -L/usr/local/lib
LIBS += -lopencv_core
LIBS += -lopencv_imgproc
LIBS += -lopencv_highgui
LIBS += -lopencv_ml
LIBS += -lopencv_video
LIBS += -lopencv_features2d
LIBS += -lopencv_calib3d
LIBS += -lopencv_objdetect
LIBS += -lopencv_contrib
LIBS += -lopencv_legacy
LIBS += -lopencv_flann




