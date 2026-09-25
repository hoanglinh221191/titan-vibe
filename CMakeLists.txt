cmake_minimum_required(VERSION 3.20)

project(titan_vibe VERSION 1.0.3 LANGUAGES CXX)

file(GLOB TITAN_VIBE_CORE_SOURCES CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/src/pbctopo/*.cpp")

add_library(titan_vibe STATIC
  ${TITAN_VIBE_CORE_SOURCES}
  src/mol_topology.cpp
  src/solv_pbc_pbctopo.cpp
)

target_compile_features(titan_vibe PUBLIC cxx_std_23)
target_include_directories(titan_vibe PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>
)
target_compile_definitions(titan_vibe PRIVATE
  TITAN_BUILD_PROFILE_TOKEN=portable
  TITAN_CXX23_FLOOR_CLANG20_GXX13=1
  TITAN_RUNTIME_DISPATCH_HOTSPOTS=0
)

if(MSVC)
  target_compile_options(titan_vibe PRIVATE /O2 /W4 /arch:AVX2)
else()
  target_compile_options(titan_vibe PRIVATE -O2 -Wall -Wextra -mavx2 -mfma)
endif()

include(GNUInstallDirs)
install(TARGETS titan_vibe
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY src/pbctopo
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.h")
install(FILES
  src/MOL.h
  src/MOL_info.h
  src/MOL_para.h
  src/Multi_MOL.h
  src/arguments.h
  src/box_type.h
  src/BoxMatrix.h
  src/index_v2.h
  src/mol_topology.h
  src/pbctopo_hull.h
  src/titan_build_info.h
  src/titan_executor.h
  src/titan_gyration.h
  src/titan_sha256.h
  src/titan_simd_dispatch.h
  src/titan_thread_compat.h
  src/titan_ui.h
  src/titan_vec3.h
  src/titan_version_info.h
  src/topology_include_utils.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
