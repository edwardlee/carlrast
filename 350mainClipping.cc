/* On Ubuntu, compile with...
    g++ 350mainClipping.cc -std=c++20 040pixel.o gl3w.o -lglfw -Ofast -march=native */

#include <GLFW/glfw3.h>

#include "280matrix.cc"
#include "150texture.cc"
#include "260shading.cc"
#include "260depth.cc"
#include "270triangle.cc"
#include "350mesh.cc"
#include "190mesh2D.cc"
#include "250mesh3D.cc"
#include "300isometry.cc"
#include "300camera.cc"
#include "340landscape.cc"
enum Attr {X, Y, Z, SA, TA, NA, OA, PA};
enum Vary {W=3, SV, TV, NV, OV, PV, Q};
enum Unif {MODELING, PROJINVISOM=16};

/* The first four entries of vary are assumed to be X, Y, Z, W. */
void shadeVertex(
        int unifDim, const double unif[], int attrDim, const double attr[], 
        int varyDim, double (&vary)[]) {
	double attrHomog[4] = {attr[X], attr[Y], attr[Z], 1.};
	double modHomog[4];
	mat441Multiply((double(*)[4])(&unif[MODELING]), attrHomog, modHomog);
	mat441Multiply((double(*)[4])(&unif[PROJINVISOM]), modHomog, (double(&)[4])vary);
	copy_n(&attr[SA], 5, &vary[SV]);
    vary[Q] = 1.;
}

void shadeFragment(
        int unifDim, const double unif[], int texNum, const texTexture *tex[], 
        int varyDim, const double vary[], double (&rgbd)[4]) {
	double sample[tex[0]->texelDim];
	texSample(tex[0], 0., vary[TV]/vary[Q], sample);
	sample[0] = sample[1] * 0.2 + 0.8;
	sample[1] = sample[1] * 0.2 + 0.6;
	double intensity = vary[PV] / sqrt(vary[NV]*vary[NV] + vary[OV]*vary[OV] + vary[PV]*vary[PV]);
    for(int i = 0; i < 3; ++i) rgbd[i] = intensity * sample[i];
	rgbd[3] = vary[Z];
}

depthBuffer buf;
shaShading sha;
texTexture texture;
const texTexture *textures[1] = {&texture};
const texTexture **tex = textures;
Landscape<LANDSIZE> landMesh;
double unif[16 + 16] = {
	1., 0., 0., 0., 
	0., 1., 0., 0., 
	0., 0., 1., 0., 
	0., 0., 0., 1., 
	1., 0., 0., 0., 
	0., 1., 0., 0., 
	0., 0., 1., 0., 
	0., 0., 0., 1.};
double viewport[4][4];
Camera cam;
double angle = M_PI_4;

void render() {
	pixClearRGB(0.8, 0.8, 1.);
	depthClearDepths(&buf, 1000000000.);
	double projInvIsom[4][4];
	cam.GetProjectionInverseIsometry(projInvIsom);
    copy_n((double *)projInvIsom, 16, &unif[PROJINVISOM]);
	landMesh.Render(&buf, viewport, &sha, unif, tex);
}

void handleKeyUp(
        int key, int shiftIsDown, int controlIsDown, int altOptionIsDown, 
        int superCommandIsDown) {
	if (key == GLFW_KEY_ENTER) {
		if (texture.filtering == texLINEAR)
			texSetFiltering(&texture, texNEAREST);
		else
			texSetFiltering(&texture, texLINEAR);
	} else if (key == GLFW_KEY_P) {
	    cam.projectionType ^= 1;
        cam.SetFrustum(M_PI / 6., 10., 10., 512, 512);
	}
}

void handleKeyDownAndRepeat(
        int key, int shiftIsDown, int controlIsDown, int altOptionIsDown, 
        int superCommandIsDown) {
    double position[3];
    copy_n(cam.isometry.translation, 3, position);
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
    /* Randomly generate a grid of elevation data. */
    Land land;
    for (double m = 0.56; m <= 1; m += 0.04)
		land.FaultRandomly(m);
	land.Blur();
	uniform_int_distribution uid(0, LANDSIZE);
	for (int i = 0; i < 4; ++i)
		land.Bump(uid(land.gen), uid(land.gen), 5., 1.);
    /* Marshal resources. */
	if (pixInitialize(512, 512, "Landscape") != 0)
		return 1;
	if (depthInitialize(&buf, 512, 512) != 0) {
	    pixFinalize();
		return 5;
	}
	if (texInitializeFile(&texture, "awesome.png") != 0) {
	    depthFinalize(&buf);
	    pixFinalize();
		return 2;
	}
    landMesh.Build(1., land.data);
	/* Manually re-assign texture coordinates. */
	for (double *v : landMesh.vert) {
	    v[TA] = v[Z];
	}
	/* Configure texture. */
    texSetFiltering(&texture, texNEAREST);
    texSetLeftRight(&texture, texREPEAT);
    texSetTopBottom(&texture, texREPEAT);
    /* Configure shader program. */
    sha.unifDim = 16 + 16;
    sha.attrDim = 3 + 2 + 3;
    sha.varyDim = 4 + 2 + 3 + 1;
    sha.shadeVertex = shadeVertex;
    sha.shadeFragment = shadeFragment;
    sha.texNum = 1;
    /* Configure viewport and camera. */
    mat44Viewport(512, 512, viewport);
    cam.projectionType = PERSPECTIVE;
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
    texFinalize(&texture);
    depthFinalize(&buf);
    pixFinalize();
}