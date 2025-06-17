#include "renderer.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// Renderer class definitions
Renderer::Renderer(
    Scene aux_scene, 
    Camera aux_camera,
    int aux_nr_rays,
    int aux_intensity_depth
) {
    scene = aux_scene;
    camera = aux_camera;
    W = camera.get_width();
    H = camera.get_height();

    // Create an image buffer with all pixels initialized to black
    image.resize(W * H * 3, 0);
    std::fill(image.begin(), image.end(), 0);

    nr_rays = aux_nr_rays;
    intensity_depth = aux_intensity_depth;
}

void Renderer::set_image_pixel_color(const int& i, const int& j, Vector& vec_albedo) {
    double gamma = 2.2;
    image[(i * W + j) * 3 + 0] = std::min(pow(vec_albedo.get_x(), 1.0 / gamma) * 254, 254.0); // Red
    image[(i * W + j) * 3 + 1] = std::min(pow(vec_albedo.get_y(), 1.0 / gamma) * 254, 254.0); // Green
    image[(i * W + j) * 3 + 2] = std::min(pow(vec_albedo.get_z(), 1.0 / gamma) * 254, 254.0); // Blue
}

Vector Renderer::get_image_pixel_color(const int& i, const int& j) {
    return Vector(image[(i * W + j) * 3 + 0], image[(i * W + j) * 3 + 1], image[(i * W + j) * 3 + 2]);
}

void Renderer::display_image(const std::string& filename) {
    stbi_write_png(filename.c_str(), W, H, 3, image.data(), 0);
}


void Renderer::render(const std::string& filename) {
    const int total_pixels = W * H;
    static std::atomic<int> pixel_count{0};
    pixel_count = 0; 

    auto start_time = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector vec_albedo_average = Vector();
            for (int k = 0; k < nr_rays; k++) {
                Ray ray = camera.get_ray(i, j);
                Intersection intersection = scene.get_closest_hit(ray);
                if (intersection.flag) {
                    Vector vec_albedo = scene.get_intensity(ray, intensity_depth, false);
                    vec_albedo_average = vec_albedo_average + vec_albedo;
                }
            }
            vec_albedo_average = (1.0 / (static_cast<double>(nr_rays))) * vec_albedo_average;
            set_image_pixel_color(i, j, vec_albedo_average);

            int local_count = ++pixel_count;
            if (local_count % 5000 == 0) {
                double fraction = double(local_count) / double(total_pixels);

                const int barWidth = 50;
                int fill = static_cast<int>(fraction * barWidth);

                #pragma omp critical
                {
                    std::cout << "\r[";
                    for (int b = 0; b < barWidth; ++b) {
                        if (b < fill)
                            std::cout << "=";  
                        else
                            std::cout << " ";  
                    }
                    std::cout << "] " 
                              << int(fraction * 100.0) << "%  " 
                              << std::flush;  
                }
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << std::endl << "Render loop took " << elapsed.count() << " seconds.\n";
    
    display_image(filename);
}

// SVG function definitions
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol) {
        FILE* f = fopen(filename.c_str(), "w+"); 
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        for (int i=0; i<polygons.size(); i++) {
            fprintf(f, "<g>\n");
            fprintf(f, "<polygon points = \""); 
            for (int j = 0; j < polygons[i].vertices.size(); j++) {
                fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
            }
            fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
            fprintf(f, "</g>\n");
        }
        fprintf(f, "</svg>\n");
        fclose(f);
}

void save_svg_animated(const std::vector<Polygon> &polygons,
                       const std::vector<std::string>& colors,
                       std::string filename, int frameid, int nbframes) {
    FILE* f;
    if (frameid == 0) {
        f = fopen(filename.c_str(), "w+");
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        fprintf(f, "<g>\n");
    } else {
        f = fopen(filename.c_str(), "a+");
    }
    fprintf(f, "<g>\n");
    for (int i = 0; i < polygons.size(); i++) {
        fprintf(f, "<polygon points = \""); 
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000-polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\" fill = \"%s\" stroke = \"black\"/>\n", colors[i].c_str());
    }
    fprintf(f, "<animate\n");
    fprintf(f, "    id = \"frame%u\"\n", frameid);
    fprintf(f, "    attributeName = \"display\"\n");
    fprintf(f, "    values = \"");
    for (int j = 0; j < nbframes; j++) {
        if (frameid == j) {
            fprintf(f, "inline");
        } else {
            fprintf(f, "none");
        }
        fprintf(f, ";");
    }
    fprintf(f, "none\"\n    keyTimes = \"");
    for (int j = 0; j < nbframes; j++) {
        fprintf(f, "%2.3f", j / (double)(nbframes));
        fprintf(f, ";");
    }
    fprintf(f, "1\"\n   dur = \"5s\"\n");
    fprintf(f, "    begin = \"0s\"\n");
    fprintf(f, "    repeatCount = \"indefinite\"/>\n");
    fprintf(f, "</g>\n");
    if (frameid == nbframes - 1) {
        fprintf(f, "</g>\n");
        fprintf(f, "</svg>\n");
    }
    fclose(f);
}

// Test case for Polygon Clipping (Sutherland–Hodgman)
void test_polygon_clipping() {
    Polygon subject({Vector(0.2, 0.2), Vector(0.8, 0.2), Vector(0.8, 0.8), Vector(0.2, 0.8)});
    Polygon clip({Vector(0.5, 0.1), Vector(0.9, 0.5), Vector(0.5, 0.9), Vector(0.1, 0.5)});
    Polygon clipped = polygon_clipping(subject, clip);
    save_svg({subject}, "../../images/svg_files/test_subject.svg", "#cccccc");
    save_svg({clip}, "../../images/svg_files/test_clip.svg", "#ffcccc");
    save_svg({clipped}, "../../images/svg_files/test_clipped.svg", "#ccffcc");
}

// Test case for Voronoi Parallel Linear Enumeration
void test_and_visualize_voronoi() {
    int num_sites = 200;
    std::vector<Vector> sites;
    sites.reserve(num_sites);
    for (int i = 0; i < num_sites; ++i) {
        Vector point = Vector((double)rand() / (double)RAND_MAX, (double)rand() / (double)RAND_MAX);
        sites.push_back(point);
    }
    auto cells = voronoi_ple(sites);
    save_svg(cells, "../../images/svg_files/voronoi.svg", "#ffcccc");
}

// Test case for Power Diagram
void run_and_visualize_power_diagram() {
    int num_sites = 100;
    std::vector<Vector> sites;
    sites.reserve(num_sites);
    for (int i = 0; i < num_sites; ++i) {
        Vector point = Vector((double)rand() / (double)RAND_MAX, (double)rand() / (double)RAND_MAX);
        sites.push_back(point);
    }
    std::vector<double> weights;
    weights.reserve(sites.size());
    weights[0] = 0.01;
    for (int  i = 1; i < sites.size(); i++) {
        weights[i] = 0.0001 * i;
    }
    auto cells = weighted_voronoi_ple(sites, weights);

    save_svg(cells, "../../images/svg_files/voronoi_pd.svg", "#ffcccc");    
}

// Test case for the Semi-Discrete Optimal Transport
void run_and_optimize_power_diagram() {
    const int N = 100;

    std::vector<Vector> local_sites;
    std::vector<double> local_lambdas;
    local_sites.reserve(N);
    local_lambdas.reserve(N);
    for (int i = 0; i < N; ++i) {
        Vector point = Vector((double)rand() / (double)RAND_MAX, (double)rand() / (double)RAND_MAX);
        local_sites.push_back(point);
        local_lambdas.push_back(1.0 / N);
    }

    ::sites      = local_sites;
    ::lambdas    = local_lambdas;

    std::cout << "Global sites size = " << ::sites.size()
              << ", lambdas size = " << ::lambdas.size() << "\n";

    std::vector<lbfgsfloatval_t> x(N, 0.0);

    lbfgs_parameter_t parameters;
    lbfgs_parameter_init(&parameters);
    parameters.max_iterations    = 500;
    parameters.epsilon           = 1e-6;
    parameters.past              = 3;
    parameters.delta             = 1e-6;
    parameters.linesearch        = LBFGS_LINESEARCH_BACKTRACKING;

    lbfgsfloatval_t final_F = 0.0;
    int ret = lbfgs(N, x.data(), &final_F, evaluate, lbfgs_progress, nullptr, &parameters);

    std::cout
      << "L-BFGS finished with code " << ret
      << ", final F = " << final_F
      << std::endl;

    std::vector<double> w(N);
    for (int i = 0; i < N; ++i) w[i] = x[i];
    auto cells = weighted_voronoi_ple(::sites, w);

    save_svg(cells, "../../images/svg_files/optimized_power_diagram.svg", "#ffcccc");
}

// Test case for the Particle Simulation using Semi-Discrete Optimal Transport
void run_and_visualize_particle_simulation() {
    const int N_air = 500;
    const int N_liquid = 150;
    const double dt = 0.03;
    const double eps = 0.05;
    const int nbframes = 50;
    ParticleSimulation sim(N_liquid, N_air, dt, eps);
    for (int frame = 0; frame < nbframes; ++frame) {
        std::vector<Polygon> polygons = sim.get_current_polygons();

        std::vector<std::string> colors;
        colors.reserve(N_liquid + N_air);
        for (int i = 0; i < N_liquid; ++i) colors.emplace_back("#0000FF"); 
        for (int i = 0; i < N_air; ++i)    colors.emplace_back("#FFFFFF"); 

        save_svg_animated(polygons, colors, "../../images/svg_files/my_fluid_animation.svg", frame, nbframes);
        sim.Gallouet_Merigot_step();
    }
}
