

add_executable(make_step_benchmark
		make_step_benchmark.cpp)

target_link_libraries(make_step_benchmark PRIVATE benchmark pthread tracer_lib)


add_executable(rand_vel_benchmark
		rand_vel_benchmark.cpp)
target_link_libraries(rand_vel_benchmark PRIVATE benchmark pthread tracer_lib)


