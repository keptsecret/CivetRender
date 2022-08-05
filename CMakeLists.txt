# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.17)

project (CIVET LANGUAGES CXX C)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

list (INSERT CMAKE_MODULE_PATH 0 "${CMAKE_SOURCE_DIR}/cmake_helpers")

set (GPU_SHADER_MODEL "" CACHE STRING "")
#add_definitions(-DCMAKE_TOOLCHAIN_FILE=C:/Users/Nut/vcpkg/scripts/buildsystems/vcpkg.cmake)

find_package( Sanitizers )
find_package ( Threads )
find_package( OpenGL REQUIRED )

set_property (GLOBAL PROPERTY USE_FOLDERS ON)

if (MSVC)
	list (APPEND CIVET_DEFINITIONS "CIVET_IS_MSVC" "_CRT_SECURE_NO_WARNINGS")
	list (APPEND CIVET_DEFINITIONS "CIVET_IS_MSVC" "_ENABLE_EXTENDED_ALIGNED_STORAGE")
endif ()

include ( CheckCXXSourceCompiles )
include ( CheckCXXSourceRuns )

########################################
# thirdparty libraries

set (BUILD_SHARED_LIBS OFF)

add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/src/thirdparty)

########################################
# os-compiler stuff

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	list (APPEND CIVET_CXX_FLAGS " -std=c++17")

	# ICC will default to -fp-model fast=1, which performs value-unsafe optimizations which will
	# cause pbrt_test to fail. For safety, -fp-model precise is explicitly set here by default.
	set (FP_MODEL "precise" CACHE STRING "The floating point model to compile with.")
	set_property (CACHE FP_MODEL PROPERTY STRINGS "precise" "fast=1" "fast=2")

	list (APPEND CIVET_CXX_FLAGS "-fp-model" "${FP_MODEL}")
endif ()

add_library (civet_warnings INTERFACE)
target_compile_options (
		civet_warnings
		INTERFACE
		"$<$<CXX_COMPILER_ID:MSVC>:$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcompiler >/wd4244>" # int -> float conversion
		"$<$<CXX_COMPILER_ID:MSVC>:$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcompiler >/wd4267>" # size_t -> int conversion
		"$<$<CXX_COMPILER_ID:MSVC>:$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcompiler >/wd4305>" # double constant assigned to float
		"$<$<CXX_COMPILER_ID:MSVC>:$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcompiler >/wd4552>" # result of expression not used
		"$<$<CXX_COMPILER_ID:MSVC>:$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcompiler >/wd4838>" # double -> int conversion
		"$<$<CXX_COMPILER_ID:MSVC>:$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcompiler >/wd4843>" # double -> float conversion
		"$<$<CXX_COMPILER_ID:MSVC>:$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcompiler >/wd26451>" # arithmetic on 4-byte value, then cast to 8-byte
		"$<$<CXX_COMPILER_ID:MSVC>:$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcompiler >/wd26495>" # uninitialized member variable
		"$<$<CXX_COMPILER_ID:MSVC>:$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcompiler >/wd4334>" # 32 to 64 bit displacement
)

CHECK_CXX_SOURCE_COMPILES (
		"int main() { constexpr int x = 0; }"
		HAVE_CONSTEXPR )
IF ( HAVE_CONSTEXPR )
	add_definitions ( -D CIVET_HAVE_CONSTEXPR )
	add_definitions ( -D CIVET_CONSTEXPR=constexpr )
ELSE ()
	add_definitions ( -D CIVET_CONSTEXPR=const )
ENDIF ()

CHECK_CXX_SOURCE_COMPILES (
		"struct alignas(32) Foo { char x; }; int main() { }"
		HAVE_ALIGNAS )
IF ( HAVE_ALIGNAS )
	add_definitions ( -D CIVET_HAVE_ALIGNAS )
ENDIF ()

CHECK_CXX_SOURCE_COMPILES (
		"int main() { int x = alignof(double); }"
		HAVE_ALIGNOF )
IF ( HAVE_ALIGNOF )
	add_definitions ( -D CIVET_HAVE_ALIGNOF )
ENDIF ()

########################################
# Aligned memory allocation

CHECK_CXX_SOURCE_COMPILES ( "
#include <malloc.h>
int main() { void * ptr = _aligned_malloc(1024, 32); }
" HAVE__ALIGNED_MALLOC )

CHECK_CXX_SOURCE_COMPILES ( "
#include <stdlib.h>
int main() {
  void *ptr;
  posix_memalign(&ptr, 32, 1024);
} " HAVE_POSIX_MEMALIGN )

CHECK_CXX_SOURCE_COMPILES ( "
#include <malloc.h>
int main() {
    void *ptr = memalign(32, 1024);
} " HAVE_MEMALIGN )

IF ( HAVE__ALIGNED_MALLOC )
	add_definitions ( -D CIVET_HAVE__ALIGNED_MALLOC )
ELSEIF ( HAVE_POSIX_MEMALIGN )
	add_definitions ( -D CIVET_HAVE_POSIX_MEMALIGN )
ELSEIF ( HAVE_MEMALIGN )
	add_definitions ( -D CIVET_HAVE_MEMALIGN )
ELSE ()
	message ( SEND_ERROR "Unable to find a way to allocate aligned memory" )
ENDIF ()

########################################
# thread-local variables

CHECK_CXX_SOURCE_COMPILES ( "
#ifdef __CYGWIN__
// Hack to work around https://gcc.gnu.org/bugzilla/show_bug.cgi?id=64697
#error \"No thread_local on cygwin\"
#endif  // __CYGWIN__
thread_local int x; int main() { }
" HAVE_THREAD_LOCAL )

CHECK_CXX_SOURCE_COMPILES ( "
__declspec(thread) int x; int main() { }
" HAVE_DECLSPEC_THREAD )

CHECK_CXX_SOURCE_COMPILES ( "
__thread int x; int main() { }
" HAVE___THREAD )

IF ( HAVE_THREAD_LOCAL )
	add_definitions ( -D CIVET_THREAD_LOCAL=thread_local )
ELSEIF ( HAVE___THREAD )
	add_definitions ( -D CIVET_THREAD_LOCAL=__thread )
ELSEIF ( HAVE_DECLSPEC_THREAD )
	add_definitions ( -D "CIVET_THREAD_LOCAL=__declspec(thread)" )
ELSE ()
	message ( SEND_ERROR "Unable to find a way to declare a thread-local variable")
ENDIF ()

###########################################################################
# CUDA

add_library (cuda_build_configuration INTERFACE)

include(CheckLanguage)
check_language(CUDA)

if (CMAKE_CUDA_COMPILER)
	find_package (CUDA REQUIRED)

	if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.17.0")
		set (CMAKE_CUDA_STANDARD 17)
	endif ()

	# https://github.com/VIAME/VIAME/blob/aa6c5f56a898b08e4da102c400b453e23952199c/CMakeLists.txt#L291
	if (NOT CUDA_VERSION_PATCH)
		if (CUDA_NVCC_EXECUTABLE AND
				CUDA_NVCC_EXECUTABLE STREQUAL CMAKE_CUDA_COMPILER AND
				CMAKE_CUDA_COMPILER_VERSION MATCHES [=[([0-9]+)\.([0-9]+)\.([0-9]+)]=])
			set (CUDA_VERSION_PATCH "${CMAKE_MATCH_3}")
		elseif (CUDA_NVCC_EXECUTABLE)
			execute_process (COMMAND ${CUDA_NVCC_EXECUTABLE} "--version" OUTPUT_VARIABLE NOUT)
			if (NOUT MATCHES [=[ V([0-9]+)\.([0-9]+)\.([0-9]+)]=])
				set (CUDA_VERSION_PATCH "${CMAKE_MATCH_3}")
			endif ()
		endif ()
	endif ()

	message (STATUS "Found CUDA: ${CUDA_VERSION_MAJOR}.${CUDA_VERSION_MINOR}.${CUDA_VERSION_PATCH}")
	if (CUDA_VERSION_MAJOR EQUAL 11 AND CUDA_VERSION_MINOR EQUAL 3 AND CUDA_VERSION_PATCH LESS 109)
		message (SEND_ERROR "Unfortunately, civet triggers an internal compiler error in CUDA 11.3.0. Please either use CUDA 11.0-11.2 or 11.3.1 or later.")
	endif ()
	enable_language (CUDA)
	list (APPEND CIVET_DEFINITIONS "CIVET_BUILD_GPU")

	set (CUDA_ENABLED ON)
	include_directories (${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

	target_compile_options (
			civet_warnings
			INTERFACE
			#"$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xptxas --warn-on-double-precision-use>"
			"$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcudafe --diag_suppress=partial_override>"
			"$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcudafe --diag_suppress=virtual_function_decl_hidden>"
			"$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcudafe --diag_suppress=integer_sign_change>"
			"$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcudafe --diag_suppress=declared_but_not_referenced>"
			# WAR invalid warnings about this with "if constexpr"
			"$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcudafe --diag_suppress=implicit_return_from_non_void_function>"
	)

	string (APPEND CMAKE_CUDA_FLAGS " -Xnvlink -suppress-stack-size-warning")

	target_compile_options (
			cuda_build_configuration
			INTERFACE
			"$<$<COMPILE_LANGUAGE:CUDA>:--std=c++17;--use_fast_math;--expt-relaxed-constexpr;--extended-lambda;--forward-unknown-to-host-compiler>"
			# The "$<NOT:$<BOOL:$<TARGET_PROPERTY:CUDA_PTX_COMPILATION>>>" part is to not add debugging symbols when generating PTX files for OptiX; see https://github.com/mmp/pbrt-v4/issues/69#issuecomment-715499748.
			"$<$<COMPILE_LANGUAGE:CUDA>:$<IF:$<AND:$<CONFIG:Debug>,$<NOT:$<BOOL:$<TARGET_PROPERTY:CUDA_PTX_COMPILATION>>>>,-G;-g,-lineinfo;-maxrregcount;128>>"
	)

	if (GPU_SHADER_MODEL STREQUAL "")
		#Get CUDA compute capability
		set (CHECK_CUDA_OUTPUT_EXE ${CMAKE_BINARY_DIR}/checkcuda)
		execute_process (COMMAND ${CMAKE_CUDA_COMPILER} -lcuda ${CMAKE_SOURCE_DIR}/cmake_helpers/checkcuda.cu -o ${CHECK_CUDA_OUTPUT_EXE}
				RESULT_VARIABLE BUILD_CHECK_CUDA_RETURN_CODE)

		if (NOT ${BUILD_CHECK_CUDA_RETURN_CODE} EQUAL 0)
			message (SEND_ERROR "Was unable to build checkcuda, consider manually setting GPU_SHADER_MODEL")
		endif ()

		execute_process (COMMAND ${CHECK_CUDA_OUTPUT_EXE}
				RESULT_VARIABLE CUDA_RETURN_CODE
				OUTPUT_VARIABLE CHECK_CUDA_OUTPUT)

		if (NOT ${CUDA_RETURN_CODE} EQUAL 0)
			message (STATUS ${CHECK_CUDA_OUTPUT})
			set (CUDA_ENABLED OFF)
		else ()
			set (ARCH "${CHECK_CUDA_OUTPUT}")
			message (STATUS "Detected CUDA Architecture: ${ARCH}")
			string (APPEND CMAKE_CUDA_FLAGS " --gpu-architecture=${ARCH}")
		endif ()
	else ()
		set (ARCH "${GPU_SHADER_MODEL}")
		message (STATUS "Specified CUDA Architecture: ${ARCH}")
		string (APPEND CMAKE_CUDA_FLAGS " --gpu-architecture=${ARCH}")
	endif ()

	# find CUDA's bin2c executable
	get_filename_component (cuda_compiler_bin "${CMAKE_CUDA_COMPILER}" DIRECTORY)
	find_program (BIN2C
			NAMES bin2c
			PATHS ${cuda_compiler_bin}
			DOC "Path to the CUDA SDK bin2c executable."
			NO_DEFAULT_PATH)
	if (NOT BIN2C)
		message (FATAL_ERROR
				"bin2c not found:\n"
				"  CMAKE_CUDA_COMPILER='${CMAKE_CUDA_COMPILER}'\n"
				"  cuda_compiler_bin='${cuda_compiler_bin}'\n"
				)
	endif ()
else()
	message(STATUS "CUDA not found")
endif()

###########################################################################
# Main civet stuff

set ( CIVET_CORE_HEADERS
		src/core/civet.h
		src/core/engine.h
		src/core/input_manager.h
		src/core/geometry/vecmath.h
		src/core/geometry/ray.h
		src/core/geometry/transform.h
		src/core/geometry/quaternion.h
		src/core/interaction.h
		src/core/shape.h
		src/core/primitive.h
		src/core/mesh.h
		src/core/material.h
		src/core/spectrum.h
		src/core/camera.h
		src/core/sampler.h
		src/core/filter.h
		src/core/film.h
		src/core/shader.h
		)

set ( CIVET_CORE_SOURCE
		src/core/engine.cpp
		src/core/input_manager.cpp
		src/core/geometry/transform.cpp
		src/core/geometry/quaternion.cpp
		src/core/interaction.cpp
		src/core/shape.cpp
		src/core/primitive.cpp
		src/core/mesh.cpp
		src/core/material.cpp
		src/core/spectrum.cpp
		src/core/camera.cpp
		src/core/sampler.cpp
		src/core/filter.cpp
		src/core/film.cpp
		src/core/shader.cpp
		)

file ( GLOB CIVET_SOURCE
		src/shapes/*
		src/utils/*
		src/cameras/*
		)

if (CUDA_ENABLED)
	set ( CIVET_GPU_SOURCE
			src/gpu/testCuda.cu
			src/gpu/testCuda.h
			)

	set_source_files_properties(
			${CIVET_CORE_HEADERS}
			${CIVET_CORE_SOURCE}
			${CIVET_SOURCE}

			PROPERTIES LANGUAGE CUDA
	)
endif()

add_library ( civet_lib STATIC
		${CIVET_CORE_HEADERS}
		${CIVET_CORE_SOURCE}
		${CIVET_GPU_SOURCE}
		${CIVET_SOURCE}
		)
add_library ( civet::civet_lib ALIAS civet_lib )

target_compile_definitions ( civet_lib PRIVATE ${CIVET_DEFINITIONS} )

set_property ( TARGET civet_lib PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON )
set_property ( TARGET civet_lib PROPERTY CUDA_ARCHITECTURES all-major )

target_include_directories( civet_lib PUBLIC
		src
		src/thirdparty
		${ASSIMP_INCLUDE}
		${STB_INCLUDE}
		${GLFW_INCLUDE}
		${GLAD_INCLUDE}
		)

target_compile_options (civet_lib PUBLIC ${CIVET_CXX_FLAGS})
target_link_libraries (civet_lib PRIVATE civet_warnings $<$<BOOL:CUDA_ENABLED>:cuda_build_configuration>)

add_sanitizers ( civet_lib )

if (WIN32)
	# Avoid a name clash when building on Visual Studio
	set_target_properties ( civet_lib PROPERTIES OUTPUT_NAME libcivet )
endif()

set ( ALL_CIVET_LIBS
		civet_lib
		${CMAKE_THREAD_LIBS_INIT}
		assimp
		glfw
		glad
		OpenGL::GL
		)

if (CUDA_ENABLED)
	set_property (TARGET civet_lib PROPERTY CUDA_SEPARABLE_COMPILATION ON)
endif()

if (WIN32)
	list (APPEND ALL_CIVET_LIBS "dbghelp" "wsock32" "ws2_32")
endif ()

###########################################################################
# Main executable

add_executable ( civet_exe src/main/civet.cpp )
add_executable ( civet::civet_exe ALIAS civet_exe )

target_compile_definitions ( civet_exe PRIVATE ${CIVET_DEFINITIONS} )
target_compile_options ( civet_exe PRIVATE ${CIVET_CXX_FLAGS} )
target_include_directories ( civet_exe PRIVATE src src/thirdparty )
target_link_libraries ( civet_exe ${ALL_CIVET_LIBS} civet_warnings )

set_target_properties ( civet_exe PROPERTIES OUTPUT_NAME civet )

add_sanitizers ( civet_exe )

set_property (TARGET civet_exe PROPERTY FOLDER "main")

# Installation

install ( TARGETS
		civet_exe
		DESTINATION
		bin
		)

install ( TARGETS
		civet_lib
		DESTINATION
		lib
		)
