#include <rendering/gbuffer.h>

namespace civet {

GBuffer::GBuffer() {
	FBO = 0;
	depth_map = 0;
}

GBuffer::~GBuffer() {
	if (!FBO) {
		glDeleteFramebuffers(1, &FBO);
	}

	if (!textures[0]) {
		glDeleteTextures(num_textures, textures);
	}

	if (!depth_map) {
		glDeleteTextures(1, &depth_map);
	}
}

bool GBuffer::init(unsigned int width, unsigned int height) {
	glGenFramebuffers(1, &FBO);
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, FBO);

	glGenTextures(num_textures, textures);
	glGenTextures(1, &depth_map);
	glGenTextures(1, &raw_texture);
	glGenTextures(1, &final_texture);

	for (int i = 0; i < num_textures; i++) {
		glBindTexture(GL_TEXTURE_2D, textures[i]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, width, height, 0, GL_RGB, GL_FLOAT, nullptr);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0 + i, GL_TEXTURE_2D, textures[i], 0);
	}

	glBindTexture(GL_TEXTURE_2D, depth_map);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH32F_STENCIL8, width, height, 0, GL_DEPTH_STENCIL, GL_FLOAT_32_UNSIGNED_INT_24_8_REV, nullptr);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_TEXTURE_2D, depth_map, 0);

	glBindTexture(GL_TEXTURE_2D, raw_texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, width, height, 0, GL_RGB, GL_FLOAT, nullptr);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT4, GL_TEXTURE_2D, raw_texture, 0);

	glBindTexture(GL_TEXTURE_2D, final_texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGB, GL_FLOAT, nullptr);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT5, GL_TEXTURE_2D, final_texture, 0);

	GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);

	if (status != GL_FRAMEBUFFER_COMPLETE) {
		printf("FB error, status: 0x%x\n", status);
		return false;
	}

	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	glCheckError("ERROR::GBuffer::init: OpenGL error code");
	return true;
}

void GBuffer::start() {
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, FBO);
	glDrawBuffer(GL_COLOR_ATTACHMENT4);
	glClear(GL_COLOR_BUFFER_BIT);
}

void GBuffer::bindGeomPass() {
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, FBO);
	GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3 };
	glDrawBuffers(num_textures, draw_buffers);
}

void GBuffer::bindStencilPass() {
	glDrawBuffer(GL_NONE);
}

void GBuffer::bindLightPass() {
	glDrawBuffer(GL_COLOR_ATTACHMENT4);

	for (unsigned int i = 0 ; i < num_textures; i++) {
		glActiveTexture(GL_TEXTURE0 + i);
		glBindTexture(GL_TEXTURE_2D, textures[GBUFFER_TEXTURE_POSITION + i]);
	}
}

void GBuffer::bindPostProcessPass() {
	glDrawBuffer(GL_COLOR_ATTACHMENT5);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, raw_texture);
}

void GBuffer::bindFinalPass() {
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	glBindFramebuffer(GL_READ_FRAMEBUFFER, FBO);
	glReadBuffer(GL_COLOR_ATTACHMENT5);
}

} // namespace civet