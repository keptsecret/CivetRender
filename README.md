# CivetRender

Civet Rendering Engine is a small rendering engine, planning to fully support in-depth rasterization and raytracing methods.  
A lot of the engine is as much a learning process for the author as a (hopefully) useful tool.  

## Features

By the end, Civet will be a fully featured rendering engine, written in OpenGL C++ and CUDA.  
The raytracing implementation is derived from [PBRT](https://github.com/mmp/pbrt-v3), and will support importing external scene files like `.obj`.
GPU accelerated rendering features is baked in but not currently fully implemented.

It is yet to be seen as to how far the author gets with this project.

### Main engine

* Vector math and matrix transforms library, fully compatible with OpenGL (with CUDA support underway)
* Load model and scene objects to triangle meshes with [assimp](https://github.com/assimp/assimp/)
* Load image textures with [stb](https://github.com/nothings/stb) 
* _(UD) GUI and improved mouse input_

### Realtime OpenGL renderer

* Diffuse and specular texture mapping using UVs
* Shadow mapping for directional lights and omnidirectional shadow mapping for point lights
* Normal mapping for surface micro details
* Forward rendering method supporting up to 4 directional lights and 32 point lights in each scene
* _(UD) forward+ rendering to support \[theoretically\] infinite lights and alpha blending_
* _(UD) Disney PBR based on Epic's Unreal Engine implementation_

![Early screenshot](./resources/screenshots/sht_76ceb66.png)

### _Ray-tracer (UD)_

* Acceleration structure with BVH
* Subdivision surfaces
* Animated transforms for rendered motion blur
* _(UD) Disney BSDF materials: diffuse, specular, glossy, transmission_
* _(UD) Path integrator and build scenes from realtime version_

### _CUDA support (UD)_

The CUDA library does not have alternatives for structures used in the ray-tracer such as `std::vector`
and smart-pointers (e.g. `std::shared_ptr`, `std::unique_ptr`).
As such, the author will have to implement custom alternatives or use different solutions altogether.
Since the current goal is to get the entire project up and running as a useable product first,
GPU compute will probably be saved for a rainy day, or until those have been completed. 

_UD: under development_
