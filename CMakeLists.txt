cmake_minimum_required(VERSION 3.14...3.30)

project(primes
        LANGUAGES CXX
        DESCRIPTION "Fast small primes generator.")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Begin of gitache configuration.
set(GITACHE_PACKAGES libcwd_r)
include(cwm4/cmake/StableGitache)
# End of gitache configuration.
 
# This project uses aicxx modules.
include(cwm4/cmake/AICxxProject)
 
#Exit if someone tries to contaminates the source directory with an in-source build.
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "Please do out-of-source builds.\nCleanup: \"rm -rf CMake CMakeCache.txt CMakeFiles/\"")
endif()
 
include(AICxxSubmodules)

add_executable(fastprimes6 fastprimes6.cxx)
target_link_libraries(fastprimes6 PRIVATE ${AICXX_OBJECTS_LIST})

add_executable(Nth-prime Nth-prime.cxx)
target_link_libraries(Nth-prime PRIVATE ${AICXX_OBJECTS_LIST})

add_executable(series series.cxx)
target_link_libraries(series PRIVATE ${AICXX_OBJECTS_LIST})

add_executable(identity identity.cxx)
target_link_libraries(identity PRIVATE ${AICXX_OBJECTS_LIST})

add_executable(atkin atkin.cxx)

add_executable(atkin3 atkin3.cxx)
target_link_libraries(atkin3 PRIVATE ${AICXX_OBJECTS_LIST})
target_compile_options(atkin3 PRIVATE -Wall -Wextra)

add_executable(sum_first_n_primes sum_first_n_primes.cxx)
target_link_libraries(sum_first_n_primes PRIVATE ${AICXX_OBJECTS_LIST})

# Add executable target from table.cpp.
add_executable(table_executable table.cpp)

# Custom target to generate table.
add_custom_target(table
  COMMAND $<TARGET_FILE:table_executable> > table
  DEPENDS table_executable
  COMMENT "Generating table using table_executable."
)
