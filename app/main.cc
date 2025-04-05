#include "renderer.h"

int main() {
    Scene scene = get_standard_scene();
    Camera camera = Camera(Vector(0, 0, 55), 60, 2048, 2048);
    Vector vec_light_source = Vector(-10, 20, 40);
    Renderer renderer = Renderer(scene, camera, vec_light_source, 40000);
    renderer.render_shading_and_shadows();
    return 0;
}
