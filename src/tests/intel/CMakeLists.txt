add_executable(test_numerical_limits test_numerical_limits.cpp)
target_link_libraries(test_numerical_limits ${myLibs})

# intel test do not make sense when using gcc
# add_executable(test_intel_compiler test_intel_compiler.cpp)
# target_link_libraries(test_intel_compiler ${myLibs})
# install(TARGETS test_intel_compiler RUNTIME DESTINATION bin)
