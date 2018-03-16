add_executable(test_intel_compiler test_intel_compiler.cpp)
target_link_libraries(test_intel_compiler ${myLibs})
install(TARGETS test_intel_compiler RUNTIME DESTINATION bin)
