cmake_minimum_required(VERSION 3.20)

# Setup Raven Project
PROJECT(raven LANGUAGES CXX)

# Fetch remote source files
include(FetchContent)
FetchContent_Populate(
    RavenSource
    URL ${RAVEN_URL}
    URL_HASH SHA256=${RAVEN_SHA256}
    SOURCE_DIR RavenSource
)

# Find NetCDF library
list(APPEND CMAKE_MODULE_PATH ${BUILD_HELPERS})
find_package(NetCDF REQUIRED)

# Find header & source
file(GLOB HEADER "${ravensource_SOURCE_DIR}/*.h")
file(GLOB SOURCE "${ravensource_SOURCE_DIR}/*.cpp")

# Add source files to build
add_executable(raven
  ${SOURCE}
  ${HEADER}
)
source_group("Header Files" FILES ${HEADER})
source_group("Source Files" FILES ${SOURCE})

# Link relevant libraries to build
target_link_libraries(raven PRIVATE NetCDF::NetCDF)
target_compile_features(raven PUBLIC cxx_std_11)
set_target_properties(raven PROPERTIES LINKER_LANGUAGE CXX)

# Install binary to environment $PATH
install(TARGETS raven DESTINATION ${SKBUILD_SCRIPTS_DIR})