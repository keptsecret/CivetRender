cmake_minimum_required (VERSION 3.17)

if (MSVC)
    add_definitions(/D _CRT_SECURE_NO_WARNINGS /Dstrdup=_strdup
            /wd4018 /wd4100 /wd4101 /wd4127 /wd4146 /wd4232 /wd4242 /wd4244 /wd4245 /wd4267 /wd4305 /wd4309
            /wd4310 /wd4334 /wd4456 /wd4464 /wd4668 /wd4701 /wd4703 /wd4711 /wd4756 /wd4820 /wd5045 /wd5250)
endif ()

###########################################################################
# assimp

set (ASSIMP_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/assimp/include PARENT_SCOPE)

set (BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set (ASSIMP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set (ASSIMP_BUILD_ASSIMP_TOOLS OFF CACHE BOOL "" FORCE)
set (ASSIMP_INSTALL OFF CACHE BOOL "" FORCE)
set (ASSIMP_BUILD_ASSIMP_VIEW OFF CACHE BOOL "" FORCE)

add_subdirectory (assimp)

set_property (TARGET assimp PROPERTY FOLDER "thirdparty")

###########################################################################
# glfw / glad

set (GLFW_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/glfw/include PARENT_SCOPE)
set (GLAD_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/glad/include PARENT_SCOPE)

set (GLFW_LIBRARY_TYPE STATIC CACHE STRING "" FORCE)
set (GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set (GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set (GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

# in the case of hybrid graphics on Windows
if (USE_HYBRID_GRAPHICS)
    message (STATUS "Setting GLFW to compile for hybrid graphics")
    set (GLFW_USE_HYBRID_HPG ON CACHE BOOL "" FORCE)
endif ()

add_subdirectory (glfw)
add_subdirectory (glad)

set_property (TARGET glfw PROPERTY FOLDER "thirdparty")
set_property (TARGET glad PROPERTY FOLDER "thirdparty")

###########################################################################
# imgui

set (IMGUI_PATH ${CMAKE_CURRENT_SOURCE_DIR}/imgui)
file (GLOB IMGUI_SOURCE ${IMGUI_PATH}/*.cpp
        ${IMGUI_PATH}/backends/imgui_impl_glfw.cpp
        ${IMGUI_PATH}/backends/imgui_impl_opengl3.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ImGuiFileDialog/*.cpp)

add_library (imgui STATIC ${IMGUI_SOURCE})
target_include_directories (imgui PUBLIC ${IMGUI_PATH} ${IMGUI_PATH}/backends)
target_link_libraries (imgui PUBLIC glfw glad)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIMGUI_IMPL_OPENGL_LOADER_GLAD2")

set (IMGUI_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/imgui PARENT_SCOPE)

# imgui filedialog extension
### uses dirent so we include that too

set (DIRENT_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/ImGuiFileDialog/dirent PARENT_SCOPE)
set (IMGUI_EXT_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/ImGuiFileDialog PARENT_SCOPE)

###########################################################################
# stb

set (STB_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/stb PARENT_SCOPE)