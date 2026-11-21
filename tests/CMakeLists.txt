# Self-test for runtime SOFA HRTF loading. Does not use googletest, so it
# builds without fetching anything.
add_executable(sofa_hrtf_test sofa_hrtf_test.cpp)
set_target_properties(sofa_hrtf_test PROPERTIES ${ALSOFT_STD_VERSION_PROPS}
    OUTPUT_NAME "sofa_hrtf_test.exe")
target_include_directories(sofa_hrtf_test PRIVATE
    ${OpenAL_BINARY_DIR}
    ${OpenAL_SOURCE_DIR}
    ${OpenAL_SOURCE_DIR}/include
    ${OpenAL_SOURCE_DIR}/common
    ${OpenAL_SOURCE_DIR}/gsl/include
    ${OpenAL_SOURCE_DIR}/fmt-11.2.0/include)
target_compile_definitions(sofa_hrtf_test PRIVATE
    AL_ALEXT_PROTOTYPES
    "SOFA_TEST_DIR=\"${OpenAL_SOURCE_DIR}/tests/data\""
    "MHR_TEST_FILE=\"${OpenAL_SOURCE_DIR}/hrtf/Default HRTF.mhr\"")
# The loader entry point is internal (hidden visibility in the shared
# library), so compile the loader sources into the test rather than exporting
# them from the public library.
target_link_libraries(sofa_hrtf_test PRIVATE OpenAL alsoft::fmt)
target_sources(sofa_hrtf_test PRIVATE
    ${OpenAL_SOURCE_DIR}/core/hrtf_loader.cpp
    ${OpenAL_SOURCE_DIR}/core/logging.cpp
    ${OpenAL_SOURCE_DIR}/common/envutils.cpp)
if(HAVE_MYSOFA)
    target_link_libraries(sofa_hrtf_test PRIVATE alsoft::mysofa)
    target_sources(sofa_hrtf_test PRIVATE ${OpenAL_SOURCE_DIR}/core/hrtf_loader_sofa.cpp)
    target_compile_definitions(sofa_hrtf_test PRIVATE HAVE_MYSOFA=1)
else()
    target_compile_definitions(sofa_hrtf_test PRIVATE HAVE_MYSOFA=0)
endif()

add_test(NAME sofa_hrtf_test COMMAND sofa_hrtf_test)

# The googletest-based tests fetch googletest at configure time, which fails
# on machines without network access. They can be disabled to build just the
# self-contained SOFA test.
option(ALSOFT_GTEST_TESTS "Build the googletest-based test programs" ON)
if(ALSOFT_GTEST_TESTS)
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        main
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(OpenAL_Tests)

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
endif()
