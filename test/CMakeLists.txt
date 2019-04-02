include_directories("${CMAKE_SOURCE_DIR}/include")

add_executable(sampler_test main.cpp cdf_tester.cpp)

target_link_libraries(sampler_test cdf_sampler)

install(TARGETS sampler_test DESTINATION bin)
