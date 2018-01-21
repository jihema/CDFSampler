add_executable(sampler_test main.cpp CDFSamplerTest.cpp)

target_link_libraries(sampler_test cdf_sampler)

install(TARGETS sampler_test DESTINATION bin)
