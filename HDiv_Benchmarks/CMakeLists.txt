
# Configure the path.h.in file.
# This file can only use the @VARIABLE@ syntax in the file
configure_file(path.h.in ${PROJECT_SOURCE_DIR}/HDiv_Benchmarks/path.h @ONLY)


file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(Hdiv3DBenchmarks ${sources} ${header})

target_link_libraries(Hdiv3DBenchmarks pz ${CBLAS_LIB} ${LAPACK_LIB}) 
