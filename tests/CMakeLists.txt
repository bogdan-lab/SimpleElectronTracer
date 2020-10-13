
add_executable(tracer_tests
                tracer_tests.cpp)

target_link_libraries(tracer_tests PUBLIC tracer_lib gtest pthread)

add_test(
	NAME AllTracerTests
        COMMAND tracer_tests
)
