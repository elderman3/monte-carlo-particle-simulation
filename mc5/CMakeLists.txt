cmake_minimum_required(VERSION 3.24)

# WSL: prefer the /usr/local CUDA toolkit over distro-packaged nvcc
set(CMAKE_CUDA_COMPILER "/usr/local/cuda-13.0/bin/nvcc" CACHE FILEPATH "")

# RTX 4070 Ti
set(CMAKE_CUDA_ARCHITECTURES 89 CACHE STRING "")

project(monte-carlo-particle-simulation LANGUAGES CXX CUDA)

add_executable(mc
  mc.cpp
  io.cpp
  geominput.cpp
  mc.h
  stats.cpp
  volume.cpp
  geomops.cu
  simulation.cu
)
target_include_directories(mc PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Standards
set_target_properties(mc PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED YES
  CUDA_STANDARD 20
  CUDA_STANDARD_REQUIRED YES
)

# Pick GPU arch (set to yours; or list multiple)
# Examples: 75(Turing), 86(Ampere), 89(Ada), 90(Hopper)
set_property(TARGET mc PROPERTY CUDA_ARCHITECTURES 89)

# If you have CUDA code split across multiple .cu files that call each other:
set_target_properties(mc PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Useful CUDA compile behaviors
set_target_properties(mc PROPERTIES
  CUDA_RUNTIME_LIBRARY Shared          # cudart DLL (default-ish on Windows)
  CUDA_EXTENSIONS OFF
)

# Warnings / optimizations / debug
if (MSVC)
  target_compile_options(mc PRIVATE
    /W4 /permissive- /EHsc
    $<$<CONFIG:Release>:/O2 /DNDEBUG>
    $<$<CONFIG:Debug>:/Od /Zi>
  )

  # Flags for NVCC when using MSVC as host compiler
  target_compile_options(mc PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math>
    $<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/W4>
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/EHsc>
    $<$<COMPILE_LANGUAGE:CUDA>:$<$<CONFIG:Release>:-O3>>
    $<$<COMPILE_LANGUAGE:CUDA>:$<$<CONFIG:Debug>:-G>>
  )
else()
  # Non-MSVC example (optional)
  target_compile_options(mc PRIVATE
    -Wall -Wextra

    # suppress unused
    -Wno-unused-variable
    -Wno-unused-parameter
    $<$<CONFIG:Release>:-O3 -DNDEBUG>
    $<$<CONFIG:Debug>:-O0 -g>
  )
endif()

# If you use __device__ symbols across translation units and hit unresolved device symbols,
# also consider enabling RDC explicitly:
set_target_properties(mc PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)


# Build as follows:
# rm -rf build
# cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
# clear; cmake --build build -v
# ./build/mc