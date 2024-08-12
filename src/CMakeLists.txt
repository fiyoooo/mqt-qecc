# main project library
add_library(
        ${PROJECT_NAME}
        ${MQT_QECC_INCLUDE_BUILD_DIR}/Code.hpp
        ${MQT_QECC_INCLUDE_BUILD_DIR}/Codes.hpp
        ${MQT_QECC_INCLUDE_BUILD_DIR}/Decoder.hpp
        ${MQT_QECC_INCLUDE_BUILD_DIR}/DecodingRunInformation.hpp
        ${MQT_QECC_INCLUDE_BUILD_DIR}/DecodingSimulator.hpp
        ${MQT_QECC_INCLUDE_BUILD_DIR}/LightsOut.hpp
        ${MQT_QECC_INCLUDE_BUILD_DIR}/MaxSATDecoder.hpp
        ${MQT_QECC_INCLUDE_BUILD_DIR}/QeccException.hpp
        ${MQT_QECC_INCLUDE_BUILD_DIR}/TreeNode.hpp
        ${MQT_QECC_INCLUDE_BUILD_DIR}/UFDecoder.hpp
        ${MQT_QECC_INCLUDE_BUILD_DIR}/UFHeuristic.hpp
        ${MQT_QECC_INCLUDE_BUILD_DIR}/Utils.hpp
        DecodingSimulator.cpp
        LightsOut.cpp
        MaxSATDecoder.cpp
        UFDecoder.cpp
        UFHeuristic.cpp
        Utils.cpp)

# set include directories
target_include_directories(${PROJECT_NAME} PUBLIC ${MQT_QECC_INCLUDE_BUILD_DIR})

# add options and warnings
target_link_libraries(
        ${PROJECT_NAME}
        PRIVATE MQT::ProjectWarnings MQT::ProjectOptions
        PUBLIC nlohmann_json::nlohmann_json ldpc)

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

if (BUILD_MQT_QECC_BINDINGS)
    add_subdirectory(python)
endif ()
