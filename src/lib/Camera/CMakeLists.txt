set(LIBRARY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/camera.cc")
set(LIBRARY_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/camera.h")
set(LIBRARY_INCLUDES "./" "${CMAKE_BINARY_DIR}/configured_files/include")

add_library("LibCamera" STATIC ${LIBRARY_SOURCES} ${LIBRARY_HEADERS})
target_include_directories("LibCamera" PUBLIC ${LIBRARY_INCLUDES})
target_link_libraries(
    "LibCamera"
    PRIVATE "LibScene"
            "LibSphere"
            "LibVector"
            "LibRay"
            "LibSTB" 
            nlohmann_json::nlohmann_json
            fmt::fmt
            spdlog::spdlog
            cxxopts::cxxopts)

if(${ENABLE_WARNINGS})
    target_set_warnings(
        TARGET
        "LibCamera"
        ENABLE
        ${ENABLE_WARNINGS}
        AS_ERRORS
        ${ENABLE_WARNINGS_AS_ERRORS})
endif()

if(${ENABLE_LTO})
    target_enable_lto(
        TARGET
        "LibCamera"
        ENABLE
        ON)
endif()

if(${ENABLE_CLANG_TIDY})
    add_clang_tidy_to_target("LibCamera")
endif()
