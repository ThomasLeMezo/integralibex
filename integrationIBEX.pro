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

QMAKE_CXXFLAGS +=-I/usr/local/include/ibex -I/usr/local/include/vtk-7.1 -O0 -DNDEBUG -Wno-deprecated -frounding-math -std=c++11 -fopenmp
QMAKE_LFLAGS +=  -fopenmp

LIBS += -L/usr/local/include/ibex -libex -lprim -lppl -lgmpxx -lgmp

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
    conversion.cpp

HEADERS  +=     border.h \
    pave.h \
    scheduler.h \
    utils.h \
    tests.h \
    /opt/VIBES/client-api/C++/src/vibes.h \
    graph.h \
    inclusion.h \
    graphdot.h \
    conversion.h


INCLUDEPATH += "/usr/local/include/vtk-7.1"
LIBS += -L/usr/local/lib \
        -lvtkalglib-7.1 \
        -lvtkChartsCore-7.1 \
        -lvtkCommonColor-7.1 \
        -lvtkCommonComputationalGeometry-7.1 \
        -lvtkCommonCore-7.1 \
        -lvtkCommonDataModel-7.1 \
        -lvtkCommonExecutionModel-7.1 \
        -lvtkCommonMath-7.1 \
        -lvtkCommonMisc-7.1 \
        -lvtkCommonSystem-7.1 \
        -lvtkCommonTransforms-7.1 \
        -lvtkDICOMParser-7.1 \
        -lvtkDomainsChemistry-7.1 \
        -lvtkDomainsChemistryOpenGL2-7.1 \
        -lvtkexoIIc-7.1 \
        -lvtkexpat-7.1 \
        -lvtkFiltersAMR-7.1 \
        -lvtkFiltersCore-7.1 \
        -lvtkFiltersExtraction-7.1 \
        -lvtkFiltersFlowPaths-7.1 \
        -lvtkFiltersGeneral-7.1 \
        -lvtkFiltersGeneric-7.1 \
        -lvtkFiltersGeometry-7.1 \
        -lvtkFiltersHybrid-7.1 \
        -lvtkFiltersHyperTree-7.1 \
        -lvtkFiltersImaging-7.1 \
        -lvtkFiltersModeling-7.1 \
        -lvtkFiltersParallel-7.1 \
        -lvtkFiltersParallelImaging-7.1 \
        -lvtkFiltersProgrammable-7.1 \
        -lvtkFiltersSelection-7.1 \
        -lvtkFiltersSMP-7.1 \
        -lvtkFiltersSources-7.1 \
        -lvtkFiltersStatistics-7.1 \
        -lvtkFiltersTexture-7.1 \
        -lvtkFiltersVerdict-7.1 \
        -lvtkfreetype-7.1 \
        -lvtkGeovisCore-7.1 \
        -lvtkglew-7.1 \
        -lvtkhdf5-7.1 \
        -lvtkhdf5_hl-7.1 \
        -lvtkImagingColor-7.1 \
        -lvtkImagingCore-7.1 \
        -lvtkImagingFourier-7.1 \
        -lvtkImagingGeneral-7.1 \
        -lvtkImagingHybrid-7.1 \
        -lvtkImagingMath-7.1 \
        -lvtkImagingMorphological-7.1 \
        -lvtkImagingSources-7.1 \
        -lvtkImagingStatistics-7.1 \
        -lvtkImagingStencil-7.1 \
        -lvtkInfovisCore-7.1 \
        -lvtkInfovisLayout-7.1 \
        -lvtkInteractionImage-7.1 \
        -lvtkInteractionStyle-7.1 \
        -lvtkInteractionWidgets-7.1 \
        -lvtkIOAMR-7.1 \
        -lvtkIOCore-7.1 \
        -lvtkIOEnSight-7.1 \
        -lvtkIOExodus-7.1 \
        -lvtkIOExport-7.1 \
        -lvtkIOGeometry-7.1 \
        -lvtkIOImage-7.1 \
        -lvtkIOImport-7.1 \
        -lvtkIOInfovis-7.1 \
        -lvtkIOLegacy-7.1 \
        -lvtkIOLSDyna-7.1 \
        -lvtkIOMINC-7.1 \
        -lvtkIOMovie-7.1 \
        -lvtkIONetCDF-7.1 \
        -lvtkIOParallel-7.1 \
        -lvtkIOParallelXML-7.1 \
        -lvtkIOPLY-7.1 \
        -lvtkIOSQL-7.1 \
        -lvtkIOVideo-7.1 \
        -lvtkIOXML-7.1 \
        -lvtkIOXMLParser-7.1 \
        -lvtkjpeg-7.1 \
        -lvtkjsoncpp-7.1 \
        -lvtklibxml2-7.1 \
        -lvtkmetaio-7.1 \
        -lvtkNetCDF-7.1 \
        -lvtkNetCDF_cxx-7.1 \
        -lvtkoggtheora-7.1 \
        -lvtkParallelCore-7.1 \
        -lvtkproj4-7.1 \
        -lvtkRenderingAnnotation-7.1 \
        -lvtkRenderingContext2D-7.1 \
        -lvtkRenderingContextOpenGL2-7.1 \
        -lvtkRenderingCore-7.1 \
        -lvtkRenderingFreeType-7.1 \
        -lvtkRenderingImage-7.1 \
        -lvtkRenderingLabel-7.1 \
        -lvtkRenderingLOD-7.1 \
        -lvtkRenderingOpenGL2-7.1 \
        -lvtkRenderingVolume-7.1 \
        -lvtkRenderingVolumeOpenGL2-7.1 \
        -lvtksqlite-7.1 \
        -lvtktiff-7.1 \
        -lvtkverdict-7.1 \
        -lvtkViewsContext2D-7.1 \
        -lvtkViewsCore-7.1 \
        -lvtkViewsInfovis-7.1 \
        -lvtkChartsCore-7.1 \
        -lvtkCommonColor-7.1 \
        -lvtkCommonComputationalGeometry-7.1 \
        -lvtkCommonCore-7.1 \
        -lvtkCommonDataModel-7.1 \
        -lvtkCommonExecutionModel-7.1 \
        -lvtkCommonMath-7.1 \
        -lvtkCommonMisc-7.1 \
        -lvtkCommonSystem-7.1 \
        -lvtkCommonTransforms-7.1 \
        -lvtkDICOMParser-7.1 \
        -lvtkDomainsChemistry-7.1 \
        -lvtkDomainsChemistryOpenGL2-7.1 \
        -lvtkexoIIc-7.1 \
        -lvtkexpat-7.1 \
        -lvtkFiltersAMR-7.1 \
        -lvtkFiltersCore-7.1 \
        -lvtkFiltersExtraction-7.1 \
        -lvtkFiltersFlowPaths-7.1 \
        -lvtkFiltersGeneral-7.1 \
        -lvtkFiltersGeneric-7.1 \
        -lvtkFiltersGeometry-7.1 \
        -lvtkFiltersHybrid-7.1 \
        -lvtkFiltersHyperTree-7.1 \
        -lvtkFiltersImaging-7.1 \
        -lvtkFiltersModeling-7.1 \
        -lvtkFiltersParallel-7.1 \
        -lvtkFiltersParallelImaging-7.1 \
        -lvtkFiltersProgrammable-7.1 \
        -lvtkFiltersSelection-7.1 \
        -lvtkFiltersSMP-7.1 \
        -lvtkFiltersSources-7.1 \
        -lvtkFiltersStatistics-7.1 \
        -lvtkFiltersTexture-7.1 \
        -lvtkFiltersVerdict-7.1 \
        -lvtkfreetype-7.1 \
        -lvtkGeovisCore-7.1 \
        -lvtkglew-7.1 \
        -lvtkhdf5-7.1 \
        -lvtkhdf5_hl-7.1 \
        -lvtkImagingColor-7.1 \
        -lvtkImagingCore-7.1 \
        -lvtkImagingFourier-7.1 \
        -lvtkImagingGeneral-7.1 \
        -lvtkImagingHybrid-7.1 \
        -lvtkImagingMath-7.1 \
        -lvtkImagingMorphological-7.1 \
        -lvtkImagingSources-7.1 \
        -lvtkImagingStatistics-7.1 \
        -lvtkImagingStencil-7.1 \
        -lvtkInfovisCore-7.1 \
        -lvtkInfovisLayout-7.1 \
        -lvtkInteractionImage-7.1 \
        -lvtkInteractionStyle-7.1 \
        -lvtkInteractionWidgets-7.1 \
        -lvtkIOAMR-7.1 \
        -lvtkIOCore-7.1 \
        -lvtkIOEnSight-7.1 \
        -lvtkIOExodus-7.1 \
        -lvtkIOExport-7.1 \
        -lvtkIOGeometry-7.1 \
        -lvtkIOImage-7.1 \
        -lvtkIOImport-7.1 \
        -lvtkIOInfovis-7.1 \
        -lvtkIOLegacy-7.1 \
        -lvtkIOLSDyna-7.1 \
        -lvtkIOMINC-7.1 \
        -lvtkIOMovie-7.1 \
        -lvtkIONetCDF-7.1 \
        -lvtkIOParallel-7.1 \
        -lvtkIOParallelXML-7.1 \
        -lvtkIOPLY-7.1 \
        -lvtkIOSQL-7.1 \
        -lvtkIOVideo-7.1 \
        -lvtkIOXML-7.1 \
        -lvtkIOXMLParser-7.1 \
        -lvtkjpeg-7.1 \
        -lvtkjsoncpp-7.1 \
        -lvtklibxml2-7.1 \
        -lvtkmetaio-7.1 \
        -lvtkNetCDF-7.1 \
        -lvtkNetCDF_cxx-7.1 \
        -lvtkoggtheora-7.1 \
        -lvtkParallelCore-7.1 \
        -lvtkproj4-7.1 \
        -lvtkRenderingAnnotation-7.1 \
        -lvtkRenderingContext2D-7.1 \
        -lvtkRenderingContextOpenGL2-7.1 \
        -lvtkRenderingCore-7.1 \
        -lvtkRenderingFreeType-7.1 \
        -lvtkRenderingImage-7.1 \
        -lvtkRenderingLabel-7.1 \
        -lvtkRenderingLOD-7.1 \
        -lvtkRenderingOpenGL2-7.1 \
        -lvtkRenderingVolume-7.1 \
        -lvtkRenderingVolumeOpenGL2-7.1 \
        -lvtksqlite-7.1 \
        -lvtktiff-7.1 \
        -lvtkverdict-7.1 \
        -lvtkViewsContext2D-7.1 \
        -lvtkViewsCore-7.1 \
        -lvtkViewsInfovis-7.1 \
        -lvtkzlib-7.1 \




