add_executable(TestMG TestMG.cpp)
target_link_libraries(TestMG bays)

add_executable(Test1DMG Test1DMG.cpp)
target_link_libraries(Test1DMG bays)

add_executable(Test2DMG Test2DMG.cpp)
target_link_libraries(Test2DMG bays)

#add_executable(TestIR TestIR.cpp)
#target_link_libraries(TestIR bays)

#add_test(NAME UnitTest COMMAND UnitTest WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)

