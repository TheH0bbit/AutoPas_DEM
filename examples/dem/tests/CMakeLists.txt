file(
    GLOB_RECURSE
    DEMTest_SRC
    "./*.cpp"
    "./*.h"
)

add_executable(demTests ${DEMTest_SRC})

target_include_directories(
    demTests
    PUBLIC ${PROJECT_SOURCE_DIR}/tests/testAutopas ${PROJECT_SOURCE_DIR}/examples/md-flexible
)

target_link_libraries(
    demTests
    PUBLIC
        autopas
        autopasTools
        gmock # gmock includes the gtest target
        yaml-cpp
)

include(GoogleTest)

gtest_discover_tests(demTests TEST_PREFIX "demTests/")