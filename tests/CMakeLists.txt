enable_testing()

add_executable(test_math test_math.cpp)
target_link_libraries(test_math gpf)
add_test(NAME test_math COMMAND test_math)

add_executable(test_disjoint_set test_disjoint_set.cpp)
target_link_libraries(test_disjoint_set gpf)
add_test(NAME test_disjoint_set COMMAND test_disjoint_set)

add_executable(test_graphcut test_graphcut.cpp)
target_link_libraries(test_graphcut gpf)
add_test(NAME test_graphcut COMMAND test_graphcut)

add_executable(test_predicates test_predicates.cpp)
target_link_libraries(test_predicates gpf)
add_test(NAME test_predicates COMMAND test_predicates)
