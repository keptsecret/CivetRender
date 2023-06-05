#ifndef CIVET_SHADER_H
#define CIVET_SHADER_H

#include <core/civet.h>
#include <core/geometry/transform.h>
#include <core/geometry/vecmath.h>
#include <glad/glad.h>

namespace civet {

class Shader {
public:
	Shader() {}
	Shader(const char* vertex_path, const char* fragment_path, const char* geometry_path = nullptr);

	void use();

	void setBool(const std::string& name, bool value) const;
	void setInt(const std::string& name, int value) const;
	void setFloat(const std::string& name, float value) const;

	void setVec2(const std::string& name, float x, float y);
	void setVec3(const std::string& name, float x, float y, float z);
	void setVec3(const std::string& name, const Vector3f& value);
	void setIVec3(const std::string& name, int x, int y, int z);
	void setIVec3(const std::string& name, const Point3i& value);
	void setVec4(const std::string& name, float x, float y, float z, float w);
	void setMat4(const std::string& name, const Matrix4& value);

	unsigned int ID;

private:
	void compileShader(const std::string& shader_path, const GLuint id);
};

} // namespace civet

#endif // CIVET_SHADER_H
