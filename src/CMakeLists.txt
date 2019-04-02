include_directories("${CMAKE_SOURCE_DIR}/include")

add_library(cdf_sampler cdf1_sampler.cpp cdf2_sampler.cpp)

install(TARGETS cdf_sampler DESTINATION lib)
