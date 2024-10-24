# collect project headers
file(GLOB MQT_QECC_HEADERS CONFIGURE_DEPENDS ${MQT_QECC_INCLUDE_BUILD_DIR}/*.hpp)
file(GLOB LDPC_HEADERS CONFIGURE_DEPENDS ${MQT_QECC_INCLUDE_BUILD_DIR}/ldpc/*.hpp)

# main project library
add_library(
        ${PROJECT_NAME}
        ${MQT_QECC_HEADERS}
        ${LDPC_HEADERS}
        ldpc/gf2sparse_linalg.cpp
        ldpc/osd.cpp
        ldpc/sort.cpp
        ldpc/util.cpp
        DecoderComparison.cpp
        DecodingSimulator.cpp
        GF2.cpp
        LightsOut.cpp
        MaxSATDecoder.cpp
        UFDecoder.cpp
        UFHeuristic.cpp
        union_find.cpp
        Utils.cpp
)

# set include directories
target_include_directories(${PROJECT_NAME} PUBLIC ${MQT_QECC_INCLUDE_BUILD_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${MQT_QECC_INCLUDE_BUILD_DIR}/ldpc/robin_map)

# add options and warnings to the library
target_link_libraries(
        ${PROJECT_NAME}
        PRIVATE MQT::ProjectWarnings MQT::ProjectOptions
        PUBLIC nlohmann_json::nlohmann_json)

# add z3 package
find_package(Z3 REQUIRED)
if (Z3_FOUND)
    message(STATUS "Found Z3: ${Z3_INCLUDE_DIRS} ${Z3_LIBRARIES}")
    target_include_directories(${PROJECT_NAME} PUBLIC ${Z3_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${Z3_LIBRARIES})
else ()
    message(FATAL_ERROR "Z3 library not found")
endif ()

# connect necessary files with z3
add_executable(lights_out LightsOut.cpp)
target_link_libraries(lights_out ${PROJECT_NAME} ${Z3_LIBRARIES})

# connect necessary files with z3
add_executable(max_sat MaxSATDecoder.cpp)
target_link_libraries(max_sat ${PROJECT_NAME} ${Z3_LIBRARIES})

# add MQT alias
add_library(MQT::QECC ALIAS ${PROJECT_NAME})

# add the ECC framework code
add_subdirectory(ecc)

if (BUILD_MQT_QECC_BINDINGS)
    add_subdirectory(python)
endif ()
