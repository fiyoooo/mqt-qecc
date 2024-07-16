# add MQT::QFR library
set(BUILD_MQT_CORE_TESTS
        OFF
        CACHE BOOL "Build MQT Core tests")
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/mqt-core" "extern/mqt-core" EXCLUDE_FROM_ALL)

# collect source and header files
file(GLOB_RECURSE HEADER_FILES "${PROJECT_SOURCE_DIR}/include/*.hpp")
file(GLOB SOURCE_FILES "${PROJECT_SOURCE_DIR}/src/*.cpp")

# main project library
add_library(${PROJECT_NAME}_lib
        ${SOURCE_FILES}
        ${HEADER_FILES}
        ${PROJECT_SOURCE_DIR}/include/ldpc/bp.hpp
        ${PROJECT_SOURCE_DIR}/include/ldpc/gf2sparse.hpp
        ${PROJECT_SOURCE_DIR}/include/ldpc/gf2sparse_linalg.hpp
        ${PROJECT_SOURCE_DIR}/include/ldpc/osd.hpp
        ${PROJECT_SOURCE_DIR}/include/ldpc/rng.hpp
        ${PROJECT_SOURCE_DIR}/include/ldpc/sort.hpp
        ${PROJECT_SOURCE_DIR}/include/ldpc/sparse_matrix_base.hpp
        ${PROJECT_SOURCE_DIR}/include/ldpc/sparse_matrix_util.hpp
        ${PROJECT_SOURCE_DIR}/include/ldpc/union_find.hpp
        ${PROJECT_SOURCE_DIR}/include/ldpc/util.hpp
)

# set include directories
target_include_directories(${PROJECT_NAME}_lib PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_BINARY_DIR}/include
)

# add options and warnings
target_link_libraries(${PROJECT_NAME}_lib PUBLIC project_options project_warnings)

# add submodules
target_link_libraries(${PROJECT_NAME}_lib PUBLIC nlohmann_json ldpc)

# add flint package
find_package(FLINT REQUIRED)
if (FLINT_FOUND)
    message(STATUS "Found FLINT: ${FLINT_INCLUDE_DIRS} ${FLINT_LIBRARIES}")
    target_include_directories(${PROJECT_NAME}_lib PUBLIC ${FLINT_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME}_lib PUBLIC ${FLINT_LIBRARIES})
else ()
    message(FATAL_ERROR "FLINT library not found")
endif ()


# add z3 package
find_package(Z3 REQUIRED)
if (Z3_FOUND)
    message(STATUS "Found Z3: ${Z3_INCLUDE_DIRS} ${Z3_LIBRARIES}")
    target_include_directories(${PROJECT_NAME}_lib PUBLIC ${Z3_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME}_lib PUBLIC ${Z3_LIBRARIES})
else ()
    message(FATAL_ERROR "Z3 library not found")
endif ()

# add MQT alias
add_library(MQT::${PROJECT_NAME}_lib ALIAS ${PROJECT_NAME}_lib)

if (BUILD_MQT_QECC_BINDINGS)
    add_subdirectory(python)
endif ()
