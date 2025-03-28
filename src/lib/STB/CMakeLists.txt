set(LIBRARY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/stb_image.cc")
set(LIBRARY_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/stb_image.h" 
"${CMAKE_CURRENT_SOURCE_DIR}/stb_image_write.h")
set(LIBRARY_INCLUDES "./" "${CMAKE_BINARY_DIR}/configured_files/include")

add_library("LibSTB" STATIC ${LIBRARY_SOURCES} ${LIBRARY_HEADERS})
target_include_directories("LibSTB" PUBLIC ${LIBRARY_INCLUDES})
target_link_libraries(
    "LibSTB"
    PRIVATE nlohmann_json::nlohmann_json
            fmt::fmt
            spdlog::spdlog
            cxxopts::cxxopts)

if(${ENABLE_WARNINGS})
    target_set_warnings(
        TARGET
        "LibSTB"
        ENABLE
        ${ENABLE_WARNINGS}
        AS_ERRORS
        ${ENABLE_WARNINGS_AS_ERRORS})
endif()

if(${ENABLE_LTO})
    target_enable_lto(
        TARGET
        "LibSTB"
        ENABLE
        ON)
endif()

if(${ENABLE_CLANG_TIDY})
    add_clang_tidy_to_target("LibSTB")
endif()
