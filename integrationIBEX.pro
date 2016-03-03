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

QMAKE_CXXFLAGS +=-I/usr/local/include/ibex -I/usr/local/include/vtk-7.0 -O0 -DNDEBUG -Wno-deprecated -frounding-math -std=c++11 -fopenmp
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


INCLUDEPATH += "/usr/local/include/vtk-7.0"
LIBS += -L/usr/local/lib \
        -lvtkalglib-7.0 \
        -lvtkChartsCore-7.0 \
        -lvtkCommonColor-7.0 \
        -lvtkCommonComputationalGeometry-7.0 \
        -lvtkCommonCore-7.0 \
        -lvtkCommonDataModel-7.0 \
        -lvtkCommonExecutionModel-7.0 \
        -lvtkCommonMath-7.0 \
        -lvtkCommonMisc-7.0 \
        -lvtkCommonSystem-7.0 \
        -lvtkCommonTransforms-7.0 \
        -lvtkDICOMParser-7.0 \
        -lvtkDomainsChemistry-7.0 \
        -lvtkDomainsChemistryOpenGL2-7.0 \
        -lvtkexoIIc-7.0 \
        -lvtkexpat-7.0 \
        -lvtkFiltersAMR-7.0 \
        -lvtkFiltersCore-7.0 \
        -lvtkFiltersExtraction-7.0 \
        -lvtkFiltersFlowPaths-7.0 \
        -lvtkFiltersGeneral-7.0 \
        -lvtkFiltersGeneric-7.0 \
        -lvtkFiltersGeometry-7.0 \
        -lvtkFiltersHybrid-7.0 \
        -lvtkFiltersHyperTree-7.0 \
        -lvtkFiltersImaging-7.0 \
        -lvtkFiltersModeling-7.0 \
        -lvtkFiltersParallel-7.0 \
        -lvtkFiltersParallelImaging-7.0 \
        -lvtkFiltersProgrammable-7.0 \
        -lvtkFiltersSelection-7.0 \
        -lvtkFiltersSMP-7.0 \
        -lvtkFiltersSources-7.0 \
        -lvtkFiltersStatistics-7.0 \
        -lvtkFiltersTexture-7.0 \
        -lvtkFiltersVerdict-7.0 \
        -lvtkfreetype-7.0 \
        -lvtkGeovisCore-7.0 \
        -lvtkglew-7.0 \
        -lvtkhdf5-7.0 \
        -lvtkhdf5_hl-7.0 \
        -lvtkImagingColor-7.0 \
        -lvtkImagingCore-7.0 \
        -lvtkImagingFourier-7.0 \
        -lvtkImagingGeneral-7.0 \
        -lvtkImagingHybrid-7.0 \
        -lvtkImagingMath-7.0 \
        -lvtkImagingMorphological-7.0 \
        -lvtkImagingSources-7.0 \
        -lvtkImagingStatistics-7.0 \
        -lvtkImagingStencil-7.0 \
        -lvtkInfovisCore-7.0 \
        -lvtkInfovisLayout-7.0 \
        -lvtkInteractionImage-7.0 \
        -lvtkInteractionStyle-7.0 \
        -lvtkInteractionWidgets-7.0 \
        -lvtkIOAMR-7.0 \
        -lvtkIOCore-7.0 \
        -lvtkIOEnSight-7.0 \
        -lvtkIOExodus-7.0 \
        -lvtkIOExport-7.0 \
        -lvtkIOGeometry-7.0 \
        -lvtkIOImage-7.0 \
        -lvtkIOImport-7.0 \
        -lvtkIOInfovis-7.0 \
        -lvtkIOLegacy-7.0 \
        -lvtkIOLSDyna-7.0 \
        -lvtkIOMINC-7.0 \
        -lvtkIOMovie-7.0 \
        -lvtkIONetCDF-7.0 \
        -lvtkIOParallel-7.0 \
        -lvtkIOParallelXML-7.0 \
        -lvtkIOPLY-7.0 \
        -lvtkIOSQL-7.0 \
        -lvtkIOVideo-7.0 \
        -lvtkIOXML-7.0 \
        -lvtkIOXMLParser-7.0 \
        -lvtkjpeg-7.0 \
        -lvtkjsoncpp-7.0 \
        -lvtklibxml2-7.0 \
        -lvtkmetaio-7.0 \
        -lvtkNetCDF-7.0 \
        -lvtkNetCDF_cxx-7.0 \
        -lvtkoggtheora-7.0 \
        -lvtkParallelCore-7.0 \
        -lvtkproj4-7.0 \
        -lvtkRenderingAnnotation-7.0 \
        -lvtkRenderingContext2D-7.0 \
        -lvtkRenderingContextOpenGL2-7.0 \
        -lvtkRenderingCore-7.0 \
        -lvtkRenderingFreeType-7.0 \
        -lvtkRenderingImage-7.0 \
        -lvtkRenderingLabel-7.0 \
        -lvtkRenderingLOD-7.0 \
        -lvtkRenderingOpenGL2-7.0 \
        -lvtkRenderingVolume-7.0 \
        -lvtkRenderingVolumeOpenGL2-7.0 \
        -lvtksqlite-7.0 \
        -lvtktiff-7.0 \
        -lvtkverdict-7.0 \
        -lvtkViewsContext2D-7.0 \
        -lvtkViewsCore-7.0 \
        -lvtkViewsInfovis-7.0 \
        -lvtkChartsCore-7.0 \
        -lvtkCommonColor-7.0 \
        -lvtkCommonComputationalGeometry-7.0 \
        -lvtkCommonCore-7.0 \
        -lvtkCommonDataModel-7.0 \
        -lvtkCommonExecutionModel-7.0 \
        -lvtkCommonMath-7.0 \
        -lvtkCommonMisc-7.0 \
        -lvtkCommonSystem-7.0 \
        -lvtkCommonTransforms-7.0 \
        -lvtkDICOMParser-7.0 \
        -lvtkDomainsChemistry-7.0 \
        -lvtkDomainsChemistryOpenGL2-7.0 \
        -lvtkexoIIc-7.0 \
        -lvtkexpat-7.0 \
        -lvtkFiltersAMR-7.0 \
        -lvtkFiltersCore-7.0 \
        -lvtkFiltersExtraction-7.0 \
        -lvtkFiltersFlowPaths-7.0 \
        -lvtkFiltersGeneral-7.0 \
        -lvtkFiltersGeneric-7.0 \
        -lvtkFiltersGeometry-7.0 \
        -lvtkFiltersHybrid-7.0 \
        -lvtkFiltersHyperTree-7.0 \
        -lvtkFiltersImaging-7.0 \
        -lvtkFiltersModeling-7.0 \
        -lvtkFiltersParallel-7.0 \
        -lvtkFiltersParallelImaging-7.0 \
        -lvtkFiltersProgrammable-7.0 \
        -lvtkFiltersSelection-7.0 \
        -lvtkFiltersSMP-7.0 \
        -lvtkFiltersSources-7.0 \
        -lvtkFiltersStatistics-7.0 \
        -lvtkFiltersTexture-7.0 \
        -lvtkFiltersVerdict-7.0 \
        -lvtkfreetype-7.0 \
        -lvtkGeovisCore-7.0 \
        -lvtkglew-7.0 \
        -lvtkhdf5-7.0 \
        -lvtkhdf5_hl-7.0 \
        -lvtkImagingColor-7.0 \
        -lvtkImagingCore-7.0 \
        -lvtkImagingFourier-7.0 \
        -lvtkImagingGeneral-7.0 \
        -lvtkImagingHybrid-7.0 \
        -lvtkImagingMath-7.0 \
        -lvtkImagingMorphological-7.0 \
        -lvtkImagingSources-7.0 \
        -lvtkImagingStatistics-7.0 \
        -lvtkImagingStencil-7.0 \
        -lvtkInfovisCore-7.0 \
        -lvtkInfovisLayout-7.0 \
        -lvtkInteractionImage-7.0 \
        -lvtkInteractionStyle-7.0 \
        -lvtkInteractionWidgets-7.0 \
        -lvtkIOAMR-7.0 \
        -lvtkIOCore-7.0 \
        -lvtkIOEnSight-7.0 \
        -lvtkIOExodus-7.0 \
        -lvtkIOExport-7.0 \
        -lvtkIOGeometry-7.0 \
        -lvtkIOImage-7.0 \
        -lvtkIOImport-7.0 \
        -lvtkIOInfovis-7.0 \
        -lvtkIOLegacy-7.0 \
        -lvtkIOLSDyna-7.0 \
        -lvtkIOMINC-7.0 \
        -lvtkIOMovie-7.0 \
        -lvtkIONetCDF-7.0 \
        -lvtkIOParallel-7.0 \
        -lvtkIOParallelXML-7.0 \
        -lvtkIOPLY-7.0 \
        -lvtkIOSQL-7.0 \
        -lvtkIOVideo-7.0 \
        -lvtkIOXML-7.0 \
        -lvtkIOXMLParser-7.0 \
        -lvtkjpeg-7.0 \
        -lvtkjsoncpp-7.0 \
        -lvtklibxml2-7.0 \
        -lvtkmetaio-7.0 \
        -lvtkNetCDF-7.0 \
        -lvtkNetCDF_cxx-7.0 \
        -lvtkoggtheora-7.0 \
        -lvtkParallelCore-7.0 \
        -lvtkproj4-7.0 \
        -lvtkRenderingAnnotation-7.0 \
        -lvtkRenderingContext2D-7.0 \
        -lvtkRenderingContextOpenGL2-7.0 \
        -lvtkRenderingCore-7.0 \
        -lvtkRenderingFreeType-7.0 \
        -lvtkRenderingImage-7.0 \
        -lvtkRenderingLabel-7.0 \
        -lvtkRenderingLOD-7.0 \
        -lvtkRenderingOpenGL2-7.0 \
        -lvtkRenderingVolume-7.0 \
        -lvtkRenderingVolumeOpenGL2-7.0 \
        -lvtksqlite-7.0 \
        -lvtktiff-7.0 \
        -lvtkverdict-7.0 \
        -lvtkViewsContext2D-7.0 \
        -lvtkViewsCore-7.0 \
        -lvtkViewsInfovis-7.0 \
        -lvtkzlib-7.0 \




