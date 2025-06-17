set(LIBRARY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/optimal_transport.cc")
set(LIBRARY_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/optimal_transport.h")
set(LIBRARY_INCLUDES "./" "${CMAKE_BINARY_DIR}/configured_files/include")

add_library("LibOT" STATIC ${LIBRARY_SOURCES} ${LIBRARY_HEADERS})
target_include_directories("LibOT" PUBLIC ${LIBRARY_INCLUDES})

find_package(OpenMP REQUIRED)

target_link_libraries(
    "LibOT"
    PRIVATE "LibVector" 
            "LibLBFGS"
            "LibSTB"
            "LibPolygon"
            "LibPS"
            nlohmann_json::nlohmann_json
            fmt::fmt
            spdlog::spdlog
            cxxopts::cxxopts
    PUBLIC  OpenMP::OpenMP_CXX)

if(${ENABLE_WARNINGS})
    target_set_warnings(
        TARGET
        "LibOT"
        ENABLE
        ${ENABLE_WARNINGS}
        AS_ERRORS
        ${ENABLE_WARNINGS_AS_ERRORS})
endif()

if(${ENABLE_LTO})
    target_enable_lto(
        TARGET
        "LibOT"
        ENABLE
        ON)
endif()

if(${ENABLE_CLANG_TIDY})
    add_clang_tidy_to_target("LibOT")
endif()
