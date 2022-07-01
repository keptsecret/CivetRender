# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.17)

project (CIVET LANGUAGES CXX C)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

find_package ( Threads )
find_package( OpenGL REQUIRED )

include ( CheckCXXSourceCompiles )

########################################
# os-compiler stuff

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	list (APPEND CIVET_CXX_FLAGS "-std=c++17")

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

	enable_language (CUDA)
	message(STATUS "CUDA Support enabled.")

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

	target_compile_options (
			cuda_build_configuration
			INTERFACE
			"$<$<COMPILE_LANGUAGE:CUDA>:--std=c++17;--use_fast_math;--expt-relaxed-constexpr;--extended-lambda;--forward-unknown-to-host-compiler>"
			# The "$<NOT:$<BOOL:$<TARGET_PROPERTY:CUDA_PTX_COMPILATION>>>" part is to not add debugging symbols when generating PTX files for OptiX; see https://github.com/mmp/pbrt-v4/issues/69#issuecomment-715499748.
			"$<$<COMPILE_LANGUAGE:CUDA>:$<IF:$<AND:$<CONFIG:Debug>,$<NOT:$<BOOL:$<TARGET_PROPERTY:CUDA_PTX_COMPILATION>>>>,-G;-g,-lineinfo;-maxrregcount;128>>"
	)
else()
	message(STATUS "CUDA not found")
endif()

###########################################################################
# Main civet stuff

set ( CIVET_SOURCE
		src/core/civet.h
		src/core/geometry.h
		src/core/ray.h)

if (CUDA_ENABLED)
	set ( CIVET_GPU_SOURCE
			src/gpu/testCuda.cu
			src/gpu/testCuda.h
			)

	set_source_files_properties(
			${CIVET_SOURCE}

			PROPERTIES LANGUAGE CUDA
	)
endif()

include_directories( src )
include_directories( src/core )

add_library ( civet_lib STATIC
		${CIVET_SOURCE}
		${CIVET_GPU_SOURCE}
		)
if (CUDA_ENABLED)
	set_target_properties(civet_lib PROPERTIES LINKER_LANGUAGE CUDA)
else()
	set_target_properties(civet_lib PROPERTIES LINKER_LANGUAGE CXX)
endif()

set ( CIVET_CXX_FEATURES
		cxx_auto_type
		cxx_explicit_conversions
		cxx_lambdas
		cxx_nullptr
		cxx_range_for
		cxx_static_assert
		)
target_compile_features ( civet_lib PRIVATE ${CIVET_CXX_FEATURES} )

IF (WIN32)
	# Avoid a name clash when building on Visual Studio
	SET_TARGET_PROPERTIES ( civet_lib PROPERTIES OUTPUT_NAME libcivet )
ENDIF()

set_property(TARGET civet_lib PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
target_compile_options(civet_lib PUBLIC ${CIVET_CXX_FLAGS})
target_link_libraries (civet_lib PRIVATE civet_warnings $<$<BOOL:CUDA_ENABLED>:cuda_build_configuration>)

SET ( ALL_CIVET_LIBS
		civet_lib
		)

add_executable ( civet_exe src/civet.cpp )

set_target_properties ( civet_exe PROPERTIES OUTPUT_NAME civet )
target_compile_features ( civet_exe PRIVATE ${CIVET_CXX11_FEATURES} )
target_link_libraries ( civet_exe ${ALL_CIVET_LIBS} )

# Installation

INSTALL ( TARGETS
		civet_exe
		DESTINATION
		bin
		)

INSTALL ( TARGETS
		civet_lib
		DESTINATION
		lib
		)
