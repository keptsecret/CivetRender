#ifndef CIVET_SHADER_H
#define CIVET_SHADER_H

#include <core/civet.h>
#include <core/geometry/transform.h>
#include <core/geometry/vecmath.h>
#include <glad/glad.h>

namespace civet {

class Shader {
public:
	Shader(const char* vertex_path, const char* fragment_path);

	virtual void use(std::vector<Point3f>&);

	void setBool(const std::string& name, bool value) const;
	void setInt(const std::string& name, int value) const;
	void setFloat(const std::string& name, float value) const;

	void setVec3(const std::string& name, float x, float y, float z);
	void setVec3(const std::string& name, const Vector3f& value);
	void setVec4(const std::string& name, float x, float y, float z, float w);
	void setMat4(const std::string& name, const Matrix4& value);

	unsigned int ID;

private:
	void compileShader(const std::string& shader_path, const GLuint id);
};

} // namespace civet

#endif // CIVET_SHADER_H
