# Standalone BS.1770-4 loudness meter acceptance test. It needs neither
# googletest nor the OpenAL library (the meter is self-contained), and uses
# no audio device (null backend): it synthesizes its own test tones.
add_executable(bs1770_test
    bs1770_test.cpp
    ${OpenAL_SOURCE_DIR}/core/bs1770.cpp
)
target_include_directories(bs1770_test PRIVATE ${OpenAL_SOURCE_DIR})
set_target_properties(bs1770_test PROPERTIES ${ALSOFT_STD_VERSION_PROPS})
add_test(NAME bs1770_test COMMAND bs1770_test)

add_executable(OpenAL_Tests)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        main
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

target_link_libraries(OpenAL_Tests PRIVATE
	OpenAL
	GTest::gtest_main
)

target_sources(OpenAL_Tests PRIVATE
example.t.cpp
)

# This needs to come last
include(GoogleTest)
gtest_discover_tests(OpenAL_Tests)
