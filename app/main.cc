#include "renderer.h"

int main() {
    Scene scene = get_intermediate_scene_refraction_reflection();
    Camera camera = Camera(Vector(0, 0, 55), 60, 2048, 2048);
    Renderer renderer = Renderer(scene, camera);
    renderer.render_refractions_reflections_shadows();
    return 0;
}
