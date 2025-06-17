set(LIBRARY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/lbfgs.cc")
set(LIBRARY_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/lbfgs.h"
"${CMAKE_CURRENT_SOURCE_DIR}/arithmetic_ansi.h")
set(LIBRARY_INCLUDES "./" "${CMAKE_BINARY_DIR}/configured_files/include")

add_library("LibLBFGS" STATIC ${LIBRARY_SOURCES} ${LIBRARY_HEADERS})
target_include_directories("LibLBFGS" PUBLIC ${LIBRARY_INCLUDES})

find_package(OpenMP REQUIRED)

target_link_libraries(
    "LibLBFGS"
    PRIVATE nlohmann_json::nlohmann_json
            fmt::fmt
            spdlog::spdlog
            cxxopts::cxxopts
    PUBLIC  OpenMP::OpenMP_CXX)

if(${ENABLE_WARNINGS})
    target_set_warnings(
        TARGET
        "LibLBFGS"
        ENABLE
        ${ENABLE_WARNINGS}
        AS_ERRORS
        ${ENABLE_WARNINGS_AS_ERRORS})
endif()

if(${ENABLE_LTO})
    target_enable_lto(
        TARGET
        "LibLBFGS"
        ENABLE
        ON)
endif()

if(${ENABLE_CLANG_TIDY})
    add_clang_tidy_to_target("LibLBFGS")
endif()
