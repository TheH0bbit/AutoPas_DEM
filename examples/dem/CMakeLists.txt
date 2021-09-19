file(
    GLOB_RECURSE
    DEM_SRC
    "*.cpp"
    "*.h"
)

add_executable(dem-test ${DEM_SRC})

target_link_libraries(dem-test PUBLIC autopas autopasTools)

# -----------------test-----------------

add_subdirectory(tests)

# add check for current target
# cmake-format: off
add_test(
        NAME dem.test
        COMMAND dem-tests
        CONFIGURATIONS checkExamples
)
# cmake-format: on

# add the executable to checkExamples as dependency
add_dependencies(checkExamples dem-test)