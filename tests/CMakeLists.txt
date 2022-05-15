enable_testing()


set(GTEST_ROOT ${PROJECT_SOURCE_DIR}/imported/googletest)
message(STATUS "GTEST_ROOT='${GTEST_ROOT}'")
include(GoogleTest)
include_directories(${GTEST_ROOT}/googletest/include)
file(GLOB TEST_SRC *.cpp)

add_executable(tests
    ${TEST_SRC}
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

target_link_libraries(tests 
    PUBLIC 
        gtest
        gtest_main
)

gtest_discover_tests(tests)
