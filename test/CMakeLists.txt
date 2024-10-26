package_add_test_with_working_dir(
        mqt-qecc-test
        MQT::QECC
        ${CMAKE_CURRENT_BINARY_DIR}
        test_decoder_comparison.cpp
        test_maxsat.cpp
        test_ufdheuristic.cpp
        test_ufdmaxsat.cpp
        test_utils.cpp
        test_originalUfd.cpp
        test_simulator.cpp
        test_treenode.cpp
        # TestBPDecoder.cpp
        TestLsd.cpp
        TestOsdDecoder.cpp
)

# files containing codes for tests
file(
        COPY "${PROJECT_SOURCE_DIR}/examples/test/"
        DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/resources/codes/"
        PATTERN "*")

package_add_test(mqt-qecc-framework-test MQT::QECCFramework test_ecc_functionality.cpp)
target_link_libraries(mqt-qecc-framework-test PRIVATE MQT::CoreDD)