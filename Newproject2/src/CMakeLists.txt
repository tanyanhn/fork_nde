add_subdirectory(Core)
add_subdirectory(RegularGrid)
add_subdirectory(Multigrid)
#add_subdirectory(IrregularGrid)

add_library(bays STATIC ${Core_Files} ${RegularGrid_Files}
  ${Multigrid_Files})
target_link_libraries(bays ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
  ${LAPACKE_LIBRARIES} ${JSONCPP_LIBRARIES} ${MPI_LIBRARIES} gfortran)
