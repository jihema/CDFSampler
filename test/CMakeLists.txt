set(TARGET test_rectified_sampling)

add_executable(${TARGET} main.cpp RectifiedSamplingTest.cpp)

target_link_libraries(${TARGET}
    ix_shader
    ix_app
    oslexec
    dneg_bxdf_clarisse
    dneg_dnpath
    dneg_clarisse_shaders
)