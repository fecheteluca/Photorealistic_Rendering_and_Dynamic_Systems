#include "renderer.h"

int main() {
    Scene scene = get_standard_scene();
    Camera camera = Camera(Vector(0, 0, 55), 60, 2048, 2048);
    Renderer renderer = Renderer(scene, camera);
    renderer.render_shading_and_shadows();
    return 0;
}
