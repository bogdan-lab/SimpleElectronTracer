
add_executable(tracer_tests
                main.cpp
                matrix_tests.cpp
		particle_tests.cpp
		surface_tests.cpp
		reflector_tests.cpp
		vector_tests.cpp
		)

target_link_libraries(tracer_tests PUBLIC tracer_lib gtest pthread)

add_test(
	NAME AllTracerTests
        COMMAND tracer_tests
)
