
add_library(tracer_lib STATIC
	    particle.cpp
	    reflector.cpp
	    surface.cpp
	    math.cpp		
)

target_link_libraries(tracer_lib PUBLIC CONAN_PKG::fmt)

set_target_properties(tracer_lib PROPERTIES
                        OUTPUT_NAME_DEBUG tracer_lib_debug
                        OUTPUT_NAME_RELEASE tracer_lib)

add_executable(pt_tracer
		main.cpp		
)

target_link_libraries(pt_tracer PRIVATE tracer_lib
                                        CONAN_PKG::nlohmann_json
                                        CONAN_PKG::lyra)

set_target_properties(pt_tracer PROPERTIES
                        OUTPUT_NAME_DEBUG pt_tracer_debug
                        OUTPUT_NAME_RELEASE pt_tracer_release)
