# CMakeLists.txt for building documentation.

#############################
# Doxygen API documentation #
#############################

find_package(Doxygen)

if(DOXYGEN_FOUND)

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile.in 
                    ${CMAKE_CURRENT_BINARY_DIR}/doxygen/Doxyfile 
                    @ONLY)
    
    add_custom_target(doxygen_doc ALL
                      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxygen/Doxyfile
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxygen
                      COMMENT "Generating API documentation with Doxygen" VERBATIM)
    
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxygen/html DESTINATION doc/doxygen)

endif(DOXYGEN_FOUND)
