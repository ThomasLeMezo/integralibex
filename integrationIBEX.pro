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

QMAKE_CXXFLAGS +=-I/usr/local/include/ibex -O0 -DNDEBUG -Wno-deprecated -frounding-math -std=c++11

LIBS += -L/usr/local/include/ibex -libex -lprim

SOURCES += main.cpp \
    border.cpp \
    pave.cpp \
    scheduler.cpp \
    utils.cpp \
    tests.cpp \
    /opt/VIBES/client-api/C++/src/vibes.cpp \
    graph.cpp

HEADERS  +=     border.h \
    pave.h \
    scheduler.h \
    utils.h \
    tests.h \
    /opt/VIBES/client-api/C++/src/vibes.h \
    graph.h




