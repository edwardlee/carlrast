/* On Ubuntu, compile with...
    g++ 350mainClipping.cc -std=c++20 040pixel.o gl3w.o -lglfw -Ofast -march=native
   On Mac, compile with...
    g++ 350mainClipping.cc -std=c++20 040pixel.o gl3w.o -lglfw3 -Ofast -march=native -framework Cocoa -framework IOKit */

#include <GLFW/glfw3.h>

#include "150texture.cc"
#include "260shading.cc"
#include "260depth.cc"
#include "270triangle.cc"
#include "280matrix.cc"
#include "350mesh.cc"
#include "250mesh3D.cc"
#include "300isometry.cc"
#include "300camera.cc"
#include "340landscape.cc"
constexpr int LANDSIZE = 40;
enum Attr {X, Y, Z, SA, TA, NA, OA, PA};
enum Vary {W=3, SV, TV, NV, OV, PV};

/* The first four entries of vary are assumed to be X, Y, Z, W. */
void shadeVertex(
        const double unif[], const double attr[], double (&vary)[]) {
	double attrHomog[4] = {attr[X], attr[Y], attr[Z], 1.};
	mat441Multiply((double(*&&)[4])(unif), attrHomog, (double(&)[4])vary);
	vary[TV] = attr[TA];
	vary[PV] = attr[PA];
}

void shadeFragment(
        const double unif[], Texture &tex, const double vary[], double (&rgbd)[4]) {
	tex.Sample(0., vary[TV]*vary[W], rgbd);
	double intensity = vary[PV] * vary[W];
    for(int i = 0; i < 2; ++i) rgbd[i] *= intensity;
	rgbd[3] = vary[Z];
}

Depth buf(512, 512);
constexpr Shading sha{16, 8, 1, 9, shadeFragment, shadeVertex};
Texture tex(1, 2, 2, {0.8, 0.5});
Land<LANDSIZE> land;
double unif[16];
double viewport[4][4];
Camera cam;
double angle = M_PI_4;

void render() {
	pixClearRGB({0.6, 0.2, 0.1});
	buf.Clear(1.);
	cam.GetProjectionInverseIsometry((double(&)[4][4])unif);
	land.Render<sha>(buf, viewport, unif, tex);
}

void handleKeyUp(int key) {
	if (key == GLFW_KEY_ENTER) {
		tex.filtering ^= 1;
	} else if (key == GLFW_KEY_P) {
	    cam.projectionType ^= 1;
        for(int i = 0; i < 4; ++i)
            cam.projection[i] *= cam.projectionType ? 0.1 : 10.;
	}
}

void handleKeyDownAndRepeat(int key) {
    double *position = cam.iso.translation;
    if (key == GLFW_KEY_W) {
        position[0] += cos(angle);
	    position[1] += sin(angle);
    } else if (key == GLFW_KEY_S) {
        position[0] -= cos(angle);
        position[1] -= sin(angle);
    } else if (key == GLFW_KEY_A)
        angle += M_PI / 120.;
    else if (key == GLFW_KEY_D)
        angle -= M_PI / 120.;
    else if (key == GLFW_KEY_Q)
        position[2] -= 1.;
    else if (key == GLFW_KEY_E)
        position[2] += 1.;
    cam.LookFrom(position, M_PI * 0.6, angle);
}

void handleTimeStep(double oldTime, double newTime) {
	if (floor(newTime) - floor(oldTime) >= 1.)
		glfwSetWindowTitle(pixWindow, (to_string(1./(newTime-oldTime))+" frames/sec").c_str());
	render();
}

int main() {
	pixInitialize(512, 512, "Landscape");
    /* Randomly generate a grid of elevation data. */
    for (double m = 0.56; m <= 1.; m += 0.04)
		land.FaultRandomly(m);
	land.Blur();
	uniform_int_distribution uid(0, LANDSIZE);
	for (int i = 0; i < 4; ++i)
		land.Bump(uid(land.gen), uid(land.gen), 5., 1.);
    /* Marshal resources. */
    land.Build();
	/* Manually re-assign texture coordinates. */
	for (auto &&v : land.vert)
	    v[TA] = v[Z];
    tex.SetTexel(0, 0, {1., 0.7});
    /* Configure viewport and camera. */
    mat44Viewport(512, 512, viewport);
    cam.SetFrustum(M_PI / 6., 10., 10., 512, 512);
    double position[3] = {-5., -5., 20.};
    cam.LookFrom(position, M_PI * 0.6, angle);
	/* Run user interface. */
    pixSetKeyDownHandler(handleKeyDownAndRepeat);
    pixSetKeyRepeatHandler(handleKeyDownAndRepeat);
    pixSetKeyUpHandler(handleKeyUp);
    pixSetTimeStepHandler(handleTimeStep);
    pixRun();
    /* Clean up. */
    pixFinalize();
}
