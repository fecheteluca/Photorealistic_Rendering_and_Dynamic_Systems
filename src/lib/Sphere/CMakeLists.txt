set(LIBRARY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/sphere.cc")
set(LIBRARY_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/sphere.h")
set(LIBRARY_INCLUDES "./" "${CMAKE_BINARY_DIR}/configured_files/include")

add_library("LibSphere" STATIC ${LIBRARY_SOURCES} ${LIBRARY_HEADERS})
target_include_directories("LibSphere" PUBLIC ${LIBRARY_INCLUDES})
target_link_libraries(
    "LibSphere"
    PRIVATE "LibVector"
            "LibRay"
            "LibGeometry"
            nlohmann_json::nlohmann_json
            fmt::fmt
            spdlog::spdlog
            cxxopts::cxxopts)

if(${ENABLE_WARNINGS})
    target_set_warnings(
        TARGET
        "LibSphere"
        ENABLE
        ${ENABLE_WARNINGS}
        AS_ERRORS
        ${ENABLE_WARNINGS_AS_ERRORS})
endif()

if(${ENABLE_LTO})
    target_enable_lto(
        TARGET
        "LibSphere"
        ENABLE
        ON)
endif()

if(${ENABLE_CLANG_TIDY})
    add_clang_tidy_to_target("LibSphere")
endif()
