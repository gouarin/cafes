#ADD_EXECUTABLE(two_part dton/2D/two_part.cpp)
#TARGET_LINK_LIBRARIES(two_part ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${VTK_LIBRARIES})

#ADD_EXECUTABLE(sem2d sem/sem.cpp)
#TARGET_LINK_LIBRARIES(sem2d ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${VTK_LIBRARIES})

ADD_EXECUTABLE(convergence dton/2D/convergence/two_part.cpp)
TARGET_LINK_LIBRARIES(convergence ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${VTK_LIBRARIES})

ADD_EXECUTABLE(convergence_same_velocities dton/2D/convergence/same_velocities/two_part_same_velocities.cpp)
TARGET_LINK_LIBRARIES(convergence_same_velocities ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${VTK_LIBRARIES})

ADD_EXECUTABLE(convergence_same_velocities_reference dton/2D/convergence/same_velocities/two_part_same_velocities_reference.cpp)
TARGET_LINK_LIBRARIES(convergence_same_velocities_reference ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${VTK_LIBRARIES})
