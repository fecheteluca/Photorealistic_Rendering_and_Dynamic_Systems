set(LIBRARY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/renderer.cc")
set(LIBRARY_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/renderer.h")
set(LIBRARY_INCLUDES "./" "${CMAKE_BINARY_DIR}/configured_files/include")

add_library("LibRenderer" STATIC ${LIBRARY_SOURCES} ${LIBRARY_HEADERS})
target_include_directories("LibRenderer" PUBLIC ${LIBRARY_INCLUDES})

find_package(OpenMP REQUIRED)

target_link_libraries(
    "LibRenderer"
    PRIVATE "LibCamera"
            "LibScene"
            "LibVector"
            "LibSphere"
            "LibRay"
            "LibSTB"
            "LibTriangleMesh"
            "LibGeometry"
            "LibPolygon"
            "LibLBFGS"
            "LibOT"
            "LibPS"
            nlohmann_json::nlohmann_json
            fmt::fmt
            spdlog::spdlog
            cxxopts::cxxopts
    PUBLIC OpenMP::OpenMP_CXX)

if(${ENABLE_WARNINGS})
    target_set_warnings(
        TARGET
        "LibRenderer"
        ENABLE
        ${ENABLE_WARNINGS}
        AS_ERRORS
        ${ENABLE_WARNINGS_AS_ERRORS})
endif()

if(${ENABLE_LTO})
    target_enable_lto(
        TARGET
        "LibRenderer"
        ENABLE
        ON)
endif()

if(${ENABLE_CLANG_TIDY})
    add_clang_tidy_to_target("LibRenderer")
endif()
