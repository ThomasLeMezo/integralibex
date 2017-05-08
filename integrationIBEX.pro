#-------------------------------------------------
#
# Project created by QtCreator 2014-11-14T09:49:56
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += qt

TARGET = integrationIBEX
TEMPLATE = app

QMAKE_CXXFLAGS +=-I/usr/local/include/ -O3 -DNDEBUG -Wno-deprecated -frounding-math -std=c++14 -fopenmp
QMAKE_LFLAGS += -fopenmp

INCLUDEPATH += /usr/local/include/
INCLUDEPATH += /opt/VIBES/client-api/C++/src/
INCLUDEPATH += /usr/local/include/ibex/
INCLUDEPATH += /usr/local/include/ibex-geometry/
LIBS += -libex
LIBS += -lprim
LIBS += -libex-geometry

SOURCES += main.cpp \
    border.cpp \
    pave.cpp \
    scheduler.cpp \
    utils.cpp \
    tests.cpp \
    /opt/VIBES/client-api/C++/src/vibes.cpp \
    graph.cpp \
    inclusion.cpp \
    graphdot.cpp

HEADERS  +=     border.h \
    pave.h \
    scheduler.h \
    utils.h \
    tests.h \
    /opt/VIBES/client-api/C++/src/vibes.h \
    graph.h \
    inclusion.h \
    graphdot.h

LIBS += -L/usr/local/lib

DISTFILES += \
    NOTES




