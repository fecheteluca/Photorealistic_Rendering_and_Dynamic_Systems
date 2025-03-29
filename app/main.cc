#include "ray_tracer.h"

int main() {
    int W = 512;
    int H = 512;

    // Create an image buffer with all pixels initialized to black
    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            // Set each pixel to red
            image[(i * W + j) * 3 + 0] = 255; // Red
            image[(i * W + j) * 3 + 1] = 0;   // Green
            image[(i * W + j) * 3 + 2] = 0;   // Blue
        }
    }
    stbi_write_png("../../images/image.png", W, H, 3, image.data(), 0);
    return 0;
}
