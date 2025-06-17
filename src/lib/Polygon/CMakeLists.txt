set(LIBRARY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/polygon.cc")
set(LIBRARY_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/polygon.h")
set(LIBRARY_INCLUDES "./" "${CMAKE_BINARY_DIR}/configured_files/include")

add_library("LibPolygon" STATIC ${LIBRARY_SOURCES} ${LIBRARY_HEADERS})
target_include_directories("LibPolygon" PUBLIC ${LIBRARY_INCLUDES})

find_package(OpenMP REQUIRED)

target_link_libraries(
    "LibPolygon"
    PRIVATE "LibVector" 
            nlohmann_json::nlohmann_json
            fmt::fmt
            spdlog::spdlog
            cxxopts::cxxopts
    PUBLIC  OpenMP::OpenMP_CXX)

if(${ENABLE_WARNINGS})
    target_set_warnings(
        TARGET
        "LibPolygon"
        ENABLE
        ${ENABLE_WARNINGS}
        AS_ERRORS
        ${ENABLE_WARNINGS_AS_ERRORS})
endif()

if(${ENABLE_LTO})
    target_enable_lto(
        TARGET
        "LibPolygon"
        ENABLE
        ON)
endif()

if(${ENABLE_CLANG_TIDY})
    add_clang_tidy_to_target("LibPolygon")
endif()
