cmake_minimum_required(VERSION 3.0)

project(sbgECom VERSION 2.0.0 DESCRIPTION "SBG-ECOM Serial Communcation Interface for SBG AHRS")

# GLOB all source files
file(GLOB_RECURSE SRC ${PROJECT_SOURCE_DIR}/src/*.c)
file(GLOB_RECURSE COMMON_SRC ${PROJECT_SOURCE_DIR}/common/*.c)
if (NOT MSVC)
	list(REMOVE_ITEM COMMON_SRC "${CMAKE_CURRENT_SOURCE_DIR}/common/interfaces/sbgInterfaceSerialWin.c")
else ()
	list(REMOVE_ITEM COMMON_SRC "${CMAKE_CURRENT_SOURCE_DIR}/common/interfaces/sbgInterfaceSerialUnix.c")
endif()

# Create library
message(LIBRARIES ${CMAKE_CURRENT_SOURCE_DIR})

px4_add_library(sbgECom ${SRC} ${COMMON_SRC})
target_include_directories(sbgECom PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/ ${CMAKE_CURRENT_SOURCE_DIR}/common/)
target_compile_options(sbgECom PUBLIC -Wno-bad-function-cast -Wno-double-promotion )

if (MSVC)
	target_link_libraries(sbgECom Ws2_32)
	target_compile_definitions(sbgECom PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

set_target_properties(sbgECom PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    PUBLIC_HEADER src/sbgECom.h)

