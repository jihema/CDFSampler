add_library(cdf_sampler cdf1_sampler.cpp cdf2_sampler.cpp)

install(TARGETS cdf_sampler DESTINATION lib)
