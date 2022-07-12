add_library (assimp STATIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/assimpd.lib)

target_include_directories(assimp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

set_property (TARGET assimp PROPERTY FOLDER "thirdparty")

set_target_properties(assimp PROPERTIES LINKER_LANGUAGE CXX)