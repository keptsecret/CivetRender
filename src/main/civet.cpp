#include <core/civet.h>

#include <core/engine.h>

using namespace civet;

int main(int argc, char* argv[]) {
	Engine* engine = Engine::getSingleton();
	if (engine->init()) {
		return -1;
	}
	engine->start();

	return 0;
}
