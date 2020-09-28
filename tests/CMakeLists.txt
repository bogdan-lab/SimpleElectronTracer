
add_executable(tracer_tests
		utils_tests.cpp)

target_link_libraries(tracer_tests PRIVATE gtest_main tracer_lib)

add_test(
	NAME AllTracerTests
	COMMAND tracer_tests
)
