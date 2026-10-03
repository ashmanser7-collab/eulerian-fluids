#include <GLFW/glfw3.h>
#include <cmath>
#include <iostream>
#include <algorithm>

// Test if divergence is set to 0,  add way of editing velocities

const double PI = 3.14159265358979323846;
const int size_x = 128;
const int size_y = 128;
const float cell_size = 8;
int width = cell_size * size_x;
int height = cell_size * size_y;
float g = 3;
float density = 1;
float dt = 0.1;


class Simulation {
    public:
    float v[size_x][size_y+1] = {{0}}; // Vertical
    float u[size_x+1][size_y] = {{0}}; // Horizontal
    float vbuffer[size_x][size_y+1] = {{0}}; // Modify buffer not direct values
    float ubuffer[size_x+1][size_y] = {{0}}; // update with alignVelocitiesToBuffers()
    float pressure[size_x][size_y] = {{0}};
    bool solidity[size_x][size_y] = {{1}};
    float fluidMap[size_x][size_y] = {{0}};
    float fluidMapBuffer[size_x][size_y] = {{0}};
    float o = 1; // over-relaxation coeficient

    Simulation() {
        resetSolidityMap();
    }
    Simulation(float value) {
        
        for (int y = 0; y < size_y; y++) {
            for (int x = 0; x < size_x; x++) {
                fluidMap[x][y] = value;
            }
        }

        resetSolidityMap();
    }
    Simulation(bool random_fluid, bool random_velocities) {
        srand(time(0));
        if (random_fluid) {
            for (int y = 0; y < size_y; y++) {
                for (int x = 0; x < size_x; x++) {
                    float value = (rand() % 1000)/1000.0f;
                    fluidMap[x][y] = value;
                }
            }
        }
        
        if (random_velocities) {
            for (int y = 0; y < size_y+1; y++) {
                for (int x = 0; x < size_x; x++) {
                    v[x][y] = (rand() % 1000)/1000.0f;
                }
            }
            for (int y = 0; y < size_y; y++) {
                for (int x = 0; x < size_x+1; x++) {
                    u[x][y] = (rand() % 1000)/1000.0f;
                }
            }
        }

        resetSolidityMap();
    }

    void alignVelocityBuffers() {
        for (int x = 0; x < size_x; x++) {
            for (int y = 0; y < size_y+1; y++) {
                v[x][y] = vbuffer[x][y];
            }
        }
        for (int x = 0; x < size_x+1; x++) {
            for (int y = 0; y < size_y; y++) {
                u[x][y] = ubuffer[x][y];
            }
        }
    }
    void resetSolidityMap() {
        for (int x = 0; x < size_x; x++) {
            for (int y = 0; y < size_y; y++) { // 0 is solid, 1 is fluid
                if (x == 0 || x == size_x-1 || y == 0 || y == size_y-1) {
                    solidity[x][y] = 0;
                } else {
                    solidity[x][y] = 1;
                }
            }
        }
    }

    bool isFluid(int x, int y) {
        bool outOfBounds = x < 0 || x >= size_x || y < 0 || y >= size_y;
        return outOfBounds ? 0 : solidity[x][y];
    }
    float getPressure(int x, int y) {
        if (x < 0 || x >= size_x || y < 0 || y >= size_y) return 0.0f;
        if (solidity[x][y] == 0) return 0.0f;
        return pressure[x][y];
    }

    float divergence(int x, int y) {
        return ((u[x+1][y] - u[x][y]) + (v[x][y+1] - v[x][y]))/cell_size;
    }

    void forceIncompressibility() {
        for (int x = 0; x < size_x; x++) {
            for (int y = 0; y < size_y; y++) {
                int stop = isFluid(x, y-1);
                int sbottom = isFluid(x, y+1);
                int sleft = isFluid(x-1, y);
                int sright = isFluid(x+1, y);

                int edges = stop + sbottom + sright + sleft;

                if (edges == 0 || isFluid(x, y) == 0) {
                    pressure[x][y] = 0;
                    continue;
                }

                float div = divergence(x, y);
                float rhs = (density / dt) * div * cell_size * cell_size;

                float psum = getPressure(x-1, y) + getPressure(x+1, y) + getPressure(x, y-1) + getPressure(x, y+1);

                pressure[x][y] = (1.0f - o) * pressure[x][y] + o * ((psum - rhs) / float(edges));
            }
        }
    }

    void updateVelocities() {
        for (int x = 0; x < size_x; x++) {
            for (int y = 0; y < size_y+1; y++) {
                if (y == 0 || y == size_y) {
                    vbuffer[x][y] = 0.0f;
                    continue;
                }
                bool topFluid  = isFluid(x, y-1);
                bool bottomFluid = isFluid(x, y);
                if (!topFluid || !bottomFluid) {
                    vbuffer[x][y] = 0.0f;
                    continue;
                }

                vbuffer[x][y] = v[x][y] - (getPressure(x, y) - getPressure(x, y-1)) * dt / cell_size / density;
            }
        }
        for (int x = 0; x < size_x+1; x++) {
            for (int y = 0; y < size_y; y++) {
                if (x == 0 || x == size_x) {
                    ubuffer[x][y] = 0.0f;
                    continue;
                }
                bool leftFluid  = isFluid(x-1, y);
                bool rightFluid = isFluid(x, y);
                if (!leftFluid || !rightFluid) {
                    ubuffer[x][y] = 0.0f;
                    continue;
                }

                ubuffer[x][y] = u[x][y] - (getPressure(x, y) - getPressure(x-1, y)) * dt / cell_size / density;
            }
        }
        alignVelocityBuffers();
    }

    void resolvePressures(int n) {
        for (int i = 0; i < n; i++) {
            forceIncompressibility();
        }
        updateVelocities();
        
    }

    void applyVerticalForce(float f) {
        for (int x = 1; x < size_x; x++) {
            for (int y = 1; y < size_y+1; y++) {
                if (!isFluid(x, y)) continue;
                v[x][y] += f*dt;
            }
        }
    }

    void advection() {
        float xvel; float yvel;
        for (int x = 1; x < size_x-1; x++) {
            for (int y = 1; y < size_y; y++) {
                if (!isFluid(x, y-1) || !isFluid(x, y)) {vbuffer[x][y] = 0.0f; continue;}
                float xpos = x * cell_size + 0.5 * cell_size;
                float ypos = y * cell_size;
                std::tie(xvel, yvel) = getVelocityAtWorldPos(xpos, ypos);
                // velocities are returned in cell-units, convert to pixels before backtracing
                xpos -= xvel * dt * cell_size;      ypos -= yvel * dt * cell_size;
                vbuffer[x][y] = getVelocityAtWorldPos(xpos, ypos).second;
            }
        }
        for (int x = 1; x < size_x; x++) {
            for (int y = 1; y < size_y-1; y++) {
                if (!isFluid(x-1, y) || !isFluid(x, y)) {ubuffer[x][y] = 0.0f; continue;}
                float xpos = x * cell_size;
                float ypos = y * cell_size + 0.5 * cell_size;
                std::tie(xvel, yvel) = getVelocityAtWorldPos(xpos, ypos);
                // convert velocity (cells) -> pixels when backtracing
                xpos -= xvel * dt * cell_size;      ypos -= yvel * dt * cell_size;
                ubuffer[x][y] = getVelocityAtWorldPos(xpos, ypos).first;
            }
        }
        alignVelocityBuffers();
    }
    std::pair<float, float> getVelocityAtWorldPos(float x, float y) {
        return {horizontalLerp(x, y), verticalLerp(x, y)};
    }
    float verticalLerp(float x, float y) {
        float xpos = (x/cell_size - 0.5); // position of x reduced to size of grid indicies
        float ypos = (y/cell_size);
        int xi = (int)floorf(xpos); // The index of the top-left corner of the box (remember, it is flipped upside down when drawing)
        int yi = (int)floorf(ypos); // Is this meant to be floor or round???
        float difx = xpos-xi; // value between 0 and 1
        float dify = ypos-yi;
        // clamp indices so we never read out-of-bounds
        if (xi < 0) xi = 0;
        if (yi < 0) yi = 0;
        if (xi > size_x-2) xi = size_x-2;
        if (yi > size_y-1) yi = size_y-1;

        float val1 = v[xi][yi] * (1-difx) + v[xi+1][yi] * (difx);
        float val2 = v[xi][yi+1] * (1-difx) + v[xi+1][yi+1] * (difx);
        return val1 * (1-dify) + val2 * (dify);
    }
    float horizontalLerp(float x, float y) {
        float xpos = (x/cell_size); // position of x reduced to size of grid indicies
        float ypos = (y/cell_size - 0.5);
        int xi = (int)floorf(xpos); // The index of the top-left corner of the box (remember, it is flipped upside down when drawing)
        int yi = (int)floorf(ypos); // Is this meant to be floor or round???
        float difx = xpos-xi; // value between 0 and 1
        float dify = ypos-yi;
        // clamp indices to valid range for u (u has size_x+1 by size_y)
        if (xi < 0) xi = 0;
        if (yi < 0) yi = 0;
        if (xi > size_x-1) xi = size_x-1;
        if (yi > size_y-2) yi = size_y-2;

        float val1 = u[xi][yi] * (1-difx) + u[xi+1][yi] * (difx);
        float val2 = u[xi][yi+1] * (1-difx) + u[xi+1][yi+1] * (difx);
        return val1 * (1-dify) + val2 * (dify);
    }

    // Fluid Stuff
    void alignFluidBuffer() {
        for (int x = 0; x < size_x; x++) {
            for (int y = 0; y < size_y; y++) {
                fluidMap[x][y] = fluidMapBuffer[x][y];
            }
        }
    }
    void advectFluid() {
        for (int x = 1; x < size_x-1; x++) {
            for (int y = 1; y < size_y-1; y++) {
                if (!isFluid(x, y)) {fluidMapBuffer[x][y] = 0; continue;}

                // Get velocity at cell center
                float xpos = (x + 0.5) * cell_size;
                float ypos = (y + 0.5) * cell_size;
                auto [xvel, yvel] = getVelocityAtWorldPos(xpos, ypos);

                // Backtrace
                xpos -= xvel * dt * cell_size;
                ypos -= yvel * dt * cell_size;

                // Interpolate fluid quantity from backtraced position
                fluidMapBuffer[x][y] = interpolateFluid(xpos, ypos);
            }
        }
        alignFluidBuffer();
    }
    // Add this helper function for bilinear interpolation of fluid quantities
    float interpolateFluid(float x, float y) {
        // Convert world position to grid coordinates
        float xpos = (x / cell_size - 0.5);
        float ypos = (y / cell_size - 0.5);

        // Get the indices of the bottom-left cell
        int xi = (int)floorf(xpos);
        int yi = (int)floorf(ypos);

        // Get interpolation weights
        float difx = xpos - xi;
        float dify = ypos - yi;

        // Clamp indices to valid range
        if (xi < 0) xi = 0;
        if (yi < 0) yi = 0;
        if (xi > size_x - 2) xi = size_x - 2;
        if (yi > size_y - 2) yi = size_y - 2;

        // Bilinear interpolation
        float val1 = fluidMap[xi][yi] * (1 - difx) + fluidMap[xi + 1][yi] * difx;
        float val2 = fluidMap[xi][yi + 1] * (1 - difx) + fluidMap[xi + 1][yi + 1] * difx;

        return val1 * (1 - dify) + val2 * dify;
    }

    // Rendering
    void renderPressures() {
        glBegin(GL_POINTS);
        for (int y = 0; y < size_y; y += 1) {
            for (int x = 0; x < size_x; x += 1) {
                float value = getPressure(x, y);
                glColor3d(value, 0.5, 1-value);
                glVertex2f((x+0.5)*cell_size, (y+0.5)*cell_size);
            }
        }
        glEnd();
    }

    void renderDivergence() {
        glBegin(GL_POINTS);
        for (int y = 0; y < size_y; y += 1) {
            for (int x = 0; x < size_x; x += 1) {
                float value = divergence(x, y);
                glColor3d(0.5+value, 0.5, 0.5-value);
                glVertex2f((x+0.5)*cell_size, (y+0.5)*cell_size);
            }
        }
        glEnd();
    }

    void renderVelocities(int coefficient) {
        glColor3d(0.8, 0.8, 0.8);
        glBegin(GL_LINES);
        for (int y = 0; y < size_y; y += 1) {
            for (int x = 0; x < size_x; x += 1) {
                float xvel; float yvel;
                float xpos = (x+0.5)*cell_size; float ypos = (y+0.5)*cell_size;
                std::tie(xvel, yvel) = getVelocityAtWorldPos(xpos, ypos);
                glVertex2f(xpos, ypos);
                glVertex2f(xpos + xvel*coefficient, ypos + yvel*coefficient);
            }
        }
        glEnd();
    }
    void renderVelocities() {
        glBegin(GL_POINTS);
        float xvel; float yvel;
        for (int y = 0; y < size_y; y += 1) {
            for (int x = 0; x < size_x; x += 1) {
                float xpos = (x+0.5)*cell_size; float ypos = (y+0.5)*cell_size;
                std::tie(xvel, yvel) = getVelocityAtWorldPos(xpos, ypos);
                float value = hypotf(xvel, yvel);
                glColor3d(value, 0.3, 1-value);
                glVertex2f(xpos, ypos);
            }
        }
        glEnd();
    }

    void renderVerticalVelocities() {
        glColor3f(0.8f, 0.4f, 0.2f);
        glBegin(GL_LINES);
        for (int y = 0; y < size_y+1; y++) {
            for (int x = 0; x < size_x; x++) {
                float xpos = (x+0.5) * cell_size; float ypos = y * cell_size;
                glVertex2f(xpos, ypos);
                glVertex2f(xpos, ypos + v[x][y] * cell_size);
            }
        }
        glEnd();
    }
    void renderHorizontalVelocities() {
        glColor3f(0.2f, 0.4f, 0.8f);
        glBegin(GL_LINES);
        for (int y = 0; y < size_y; y++) {
            for (int x = 0; x < size_x+1; x++) {
                float xpos = x * cell_size; float ypos = (y+0.5) * cell_size;
                glVertex2f(xpos, ypos);
                glVertex2f(xpos + u[x][y] * cell_size, ypos);
            }
        }
        glEnd();
    }
    void renderFluid() {
        glBegin(GL_POINTS);
        for (int y = 0; y < size_y; y += 1) {
            for (int x = 0; x < size_x; x += 1) {
                float value = fluidMap[x][y]; // Needs to be clamped between 1 and 0
                value = 0.8*value;
                glColor3d(value, value, value);
                glVertex2f((x+0.5)*cell_size, (y+0.5)*cell_size);
            }
        }
        glEnd();
    }
    void renderSolidity() {
        glBegin(GL_POINTS);
        for (int y = 0; y < size_y; y += 1) {
            for (int x = 0; x < size_x; x += 1) {
                float value = fluidMap[x][y]; // Needs to be clamped between 1 and 0
                value = isFluid(x, y)*0.5;
                glColor3d(value, value, value);
                glVertex2f((x+0.5)*cell_size, (y+0.5)*cell_size);
            }
        }
        glEnd();
    }
};

bool mouse_down = false;
void mouse_button(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
        mouse_down = true;
    }
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
        mouse_down = false;
    }
}


int main() {
    srand(time(0));

    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return -1;
    }

    GLFWwindow* window = glfwCreateWindow(width, height, "Eulerian Fluid Simulation", NULL, NULL);

    glfwMakeContextCurrent(window);

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, width, 0, height, 1, -1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPointSize(cell_size);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glfwSetMouseButtonCallback(window, NULL);
    glfwSetWindowPos(window, 10, 75);

    glfwSwapBuffers(window);
    glfwPollEvents();

    Simulation world = Simulation();
    int rad = 5;
    int xc = size_x/2;
    int yc = size_y/2;
    for (int x = -rad; x <= rad; x++) {
        for (int y = -rad; y <= rad; y++) {
            float dist = hypot(x, y);
            if (dist <= rad) {
                world.solidity[xc+x][yc+y] = 0;
                world.solidity[xc+x][yc+y] = 0;
            }
        }
    }
    float mouse_radius = 3;
    float mouse_strength = 2;
    int mouse_mode = 1;
    
    double oldxpos;
    double oldypos;
    double xpos;
    double ypos;

    bool wind_tunnel_velocities = false;
    bool wind_tunnel_smoke = true;
    int display = 0;
    float emission_size = 10;
    float smoke_emission = 0.05;
    float velocity_emission = 5;

    float *affected = &mouse_radius;
    float increment = 1;
    std::string affecting = "mouse radius";
    
    glfwSetMouseButtonCallback(window, mouse_button);

    while (!glfwWindowShouldClose(window))
    {
        glfwGetCursorPos(window, &xpos, &ypos);
        ypos = height - ypos;
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
            world = Simulation();
        }
        if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS) {
            mouse_mode = 1;
            std::cout << "Velocity brush" << std::endl;
        }
        if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS) {
            mouse_mode = 2;
            std::cout << "Smoke brush" << std::endl;
        }
        if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS) {
            mouse_mode = 3;
            std::cout << "Solidity brush" << std::endl;
        }


        if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS) {
            affected = &emission_size;
            increment = 1.0f;
            affecting = "emission size";
            std::cout << "Now affecting " << affecting << std::endl;
        }
        if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS) {
            affected = &mouse_radius;
            increment = 1;
            affecting = "mouse radius";
            std::cout << "Now affecting " << affecting << std::endl;
        }
        if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS) {
            affected = &mouse_strength;
            increment = 0.01f;
            affecting = "mouse strength";
            std::cout << "Now affecting " << affecting << std::endl;
        }
        if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS) {
            affected = &smoke_emission;
            increment = 0.005f;
            affecting = "smoke emission";
            std::cout << "Now affecting " << affecting << std::endl;
        }
        if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS) {
            affected = &velocity_emission;
            increment = 0.1f;
            affecting = "velocity emission";
            std::cout << "Now affecting " << affecting << std::endl;
        }
        if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
            *affected += increment;
            std::cout << affecting << ": " << *affected << std::endl;
        }
        if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS) {
            *affected -= increment;
            std::cout << affecting << ": " << *affected << std::endl;
        }


        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
            wind_tunnel_smoke = !wind_tunnel_smoke;
            std::cout << "Wind tunnel smoke " << (wind_tunnel_smoke ? "on" : "off") << std::endl;
            _sleep(100);
        }
        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
            wind_tunnel_velocities = !wind_tunnel_velocities;
            std::cout << "Wind tunnel velocities " << (wind_tunnel_velocities ? "on" : "off") << std::endl;
            _sleep(100);
        }
        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
            display = (display+1)%6;
            std::cout << "Display changed" << std::endl;
            _sleep(100);
        }

        if (mouse_down) {
            float dx = xpos - oldxpos;
            float dy = ypos - oldypos;
            if (mouse_mode == 1) {
                for (int x = -mouse_radius; x <= mouse_radius; x++) {
                    for (int y = -mouse_radius; y <= mouse_radius; y++) {
                        float dist = hypot(x, y);
                        if (dist <= mouse_radius) {
                            world.v[(int)roundf(xpos/cell_size)+x][(int)roundf(ypos/cell_size)+y] += dy/cell_size/mouse_radius/mouse_radius*mouse_strength;
                            world.u[(int)roundf(xpos/cell_size)+x][(int)roundf(ypos/cell_size)+y] += dx/cell_size/mouse_radius/mouse_radius*mouse_strength;
                        }
                    }
                }
            }
            if (mouse_mode == 2) {
                for (int x = -mouse_radius; x <= mouse_radius; x++) {
                    for (int y = -mouse_radius; y <= mouse_radius; y++) {
                        float dist = hypot(x, y);
                        if (dist <= mouse_radius) {
                            world.fluidMap[(int)roundf(xpos/cell_size)+x][(int)roundf(ypos/cell_size)+y] += dt*mouse_strength;
                        }
                    }
                }
            }
            if (mouse_mode == 3) {
                for (int x = -mouse_radius; x <= mouse_radius; x++) {
                    for (int y = -mouse_radius; y <= mouse_radius; y++) {
                        float dist = hypot(x, y);
                        if (dist <= mouse_radius) {
                            world.solidity[(int)roundf(xpos/cell_size)+x][(int)roundf(ypos/cell_size)+y] = 0;
                            world.fluidMap[(int)roundf(xpos/cell_size)+x][(int)roundf(ypos/cell_size)+y] = 0;
                        }
                    }
                }
            }
        }
        oldxpos = xpos;
        oldypos = ypos;
        if (wind_tunnel_velocities) {
            for (int y = size_y/2 - emission_size/2; y < size_y/2 + emission_size/2; y++) {
                world.u[2][y] += velocity_emission;
            }
        }
        if (wind_tunnel_smoke) {
            for (int y = size_y/2 - emission_size/2; y < size_y/2 + emission_size/2; y++) {
                world.fluidMap[4][y] += smoke_emission;
            }
        }

        if (display == 0) {world.renderFluid();}
        if (display == 1) {world.renderVelocities();}
        if (display == 2) {world.renderVelocities(cell_size);}
        if (display == 3) {world.renderPressures();}
        if (display == 4) {world.renderDivergence();}
        if (display == 5) {world.renderSolidity();}
        
        
        world.advection();
        world.resolvePressures(2);
        world.advectFluid();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}