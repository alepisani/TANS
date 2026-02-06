cmake_minimum_required(VERSION 3.15)
project(MyROOTProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Trova ROOT
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist)
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Geom RGL)

# Include directories di ROOT
include(${ROOT_USE_FILE})

# Aggiungi la directory include
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# File sorgenti
set(SOURCES
    src/main.cpp
    src/point.cpp
    src/particle.cpp
    src/event.cpp
    )

# Header files per il dizionario
set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/include/point.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/particle.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/event.h)

# Genera il dizionario ROOT
ROOT_GENERATE_DICTIONARY(G__MyDict ${HEADERS} LINKDEF src/LinkDef.h)

# Crea l'eseguibile
add_executable(sim ${SOURCES} G__MyDict.cxx)

# Linka le librerie di ROOT
target_link_libraries(sim ${ROOT_LIBRARIES})