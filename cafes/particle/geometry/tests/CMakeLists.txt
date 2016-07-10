set(TEST_SOURCES
  box.cpp
  )

add_executable(geometry_test ${TEST_SOURCES})
target_link_libraries(geometry_test Catch)

add_test(NAME GeometryTest COMMAND geometry_test)