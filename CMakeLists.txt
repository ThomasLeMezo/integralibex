project(integrationIBEX)
cmake_minimum_required(VERSION 2.8)

##### IBEX #####
include_directories(${CMAKE_CURRENT_LIST_DIR} "/home/lemezoth/Ibex/soplex-1.7.2/lib /home/lemezoth/Ibex/soplex-1.7.2/src /usr/local/include/ibex")
FILE(GLOB LOCAL_LIBRARY_DIRS /usr/local/include/ibex /home/lemezoth/Ibex/soplex-1.7.2/lib /home/lemezoth/Ibex/soplex-1.7.2/src)
INCLUDE_DIRECTORIES(${LOCAL_LIBRARY_DIRS})
LINK_DIRECTORIES(${LOCAL_LIBRARY_DIRS})

##### VIBES #####
SET(vibes_SOURCES /opt/VIBES/client-api/C++/src/vibes.cpp)
INCLUDE_DIRECTORIES(/opt/VIBES/client-api/C++/src)

SET(SRC
    main.cpp
    border.cpp
    pave.cpp
    ${vibes_SOURCES}
)

SET(HEADERS
    border.h
    pave.h

)

ADD_EXECUTABLE(${PROJECT_NAME} ${SRC})
TARGET_LINK_LIBRARIES(${PROJECT_NAME}
    ${IBEX}
    ${SOPLEX}
    ${CLP}
    ${PRIM})
