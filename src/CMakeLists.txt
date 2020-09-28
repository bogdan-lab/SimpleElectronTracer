
add_library(tracer_lib STATIC
	    particle.cpp
	    reflector.cpp
	    surface.cpp
	    utils.cpp		
)

add_executable(pt_tracer
		main.cpp		
)

target_link_libraries(pt_tracer PRIVATE tracer_lib)

