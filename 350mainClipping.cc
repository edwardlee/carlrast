/* On Ubuntu, compile with...
    g++ 350mainClipping.cc -std=c++20 040pixel.o gl3w.o -lglfw -Ofast -march=native */

#include <cstdio>
#include <ctime>
#include <GLFW/glfw3.h>

#include "040pixel.h"

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

#define LANDSIZE 40

#define ATTRX 0
#define ATTRY 1
#define ATTRZ 2
#define ATTRS 3
#define ATTRT 4
#define ATTRN 5
#define ATTRO 6
#define ATTRP 7
#define VARYX 0
#define VARYY 1
#define VARYZ 2
#define VARYW 3
#define VARYS 4
#define VARYT 5
#define VARYN 6
#define VARYO 7
#define VARYP 8
#define VARY1 9
#define UNIFMODELING 0
#define UNIFPROJINVISOM 16
#define TEXR 0
#define TEXG 1
#define TEXB 2

/* The first four entries of vary are assumed to be X, Y, Z, W. */
void shadeVertex(
        int unifDim, const double unif[], int attrDim, const double attr[], 
        int varyDim, double (&vary)[]) {
	double attrHomog[4] = {attr[ATTRX], attr[ATTRY], attr[ATTRZ], 1.};
	double modHomog[4];
	mat441Multiply((double(*)[4])(&unif[UNIFMODELING]), attrHomog, modHomog);
	mat441Multiply((double(*)[4])(&unif[UNIFPROJINVISOM]), modHomog, (double(&)[4])vary);
	std::copy_n(&attr[ATTRS], 5, &vary[VARYS]);
    vary[VARY1] = 1.;
}

void shadeFragment(
        int unifDim, const double unif[], int texNum, const texTexture *tex[], 
        int varyDim, const double vary[], double (&rgbd)[4]) {
    double vary1[varyDim];
    std::copy_n(vary, varyDim, vary1);
    for(int i = VARYS; i < varyDim; ++i) vary1[i] /= vary[VARY1];
	double sample[tex[0]->texelDim];
	texSample(tex[0], vary1[VARYS], vary1[VARYT], sample);
	sample[0] = sample[1] * 0.2 + 0.8;
	sample[1] = sample[1] * 0.2 + 0.6;
	sample[2] = 0.3;
	double intensity = vary1[VARYP];
	double l = vary1[VARYN]*vary1[VARYN] + vary1[VARYO]*vary1[VARYO] + vary1[VARYP]*vary1[VARYP];
	intensity /= sqrt(l);
    for(int i = 0; i < 3; ++i) rgbd[i] = intensity * sample[i];
	rgbd[3] = vary[VARYZ];
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
camCamera cam;
double angle = M_PI * 0.25;

void render(void) {
	pixClearRGB(0.8, 0.8, 1.);
	depthClearDepths(&buf, 1000000000.);
	double projInvIsom[4][4];
	camGetProjectionInverseIsometry(&cam, projInvIsom);
    std::copy_n((double *)projInvIsom, 16, &unif[UNIFPROJINVISOM]);
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
	    if (cam.projectionType == camORTHOGRAPHIC)
		    camSetProjectionType(&cam, camPERSPECTIVE);
		else
		    camSetProjectionType(&cam, camORTHOGRAPHIC);
        camSetFrustum(&cam, M_PI / 6., 10., 10., 512, 512);
	}
}

void handleKeyDownAndRepeat(
        int key, int shiftIsDown, int controlIsDown, int altOptionIsDown, 
        int superCommandIsDown) {
    double position[3];
    std::copy_n(cam.isometry.translation, 3, position);
    if (key == GLFW_KEY_W) {
        double delta[3] = {cos(angle), sin(angle), 0.};
        for(int i = 0; i < 3; ++i) position[i] += delta[i];
    } else if (key == GLFW_KEY_S) {
        double delta[3] = {cos(angle), sin(angle), 0.};
        for(int i = 0; i < 3; ++i) position[i] -= delta[i];
    } else if (key == GLFW_KEY_A)
        angle += M_PI / 120.;
    else if (key == GLFW_KEY_D)
        angle -= M_PI / 120.;
    else if (key == GLFW_KEY_Q)
        position[2] -= 1.;
    else if (key == GLFW_KEY_E)
        position[2] += 1.;
    camLookFrom(&cam, position, M_PI * 0.6, angle);
}

void handleTimeStep(double oldTime, double newTime) {
	if (floor(newTime) - floor(oldTime) >= 1.)
		printf("handleTimeStep: %f frames/sec\n", 1. / (newTime - oldTime));
	render();
}

int main() {
    /* Randomly generate a grid of elevation data. */
    double landData[LANDSIZE * LANDSIZE];
    landFlat(LANDSIZE, landData, 0.);
    time_t t;
	srand((unsigned)time(&t));
    for (int i = 0; i < 12; i += 1)
		landFaultRandomly(LANDSIZE, (double *)landData, 1. - i * 0.04);
	for (int i = 0; i < 4; i += 1)
		landBlur(LANDSIZE, (double *)landData);
	for (int i = 0; i < 4; i += 1)
		landBump(LANDSIZE, (double *)landData, landInt(0, LANDSIZE - 1), 
		    landInt(0, LANDSIZE - 1), 5., 1.);
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
    landMesh.Build(1., landData);
	/* Manually re-assign texture coordinates. */
	for (int i = 0; i < landMesh.vert.size(); i += 1) {
	    landMesh.vert[i][ATTRS] = 0.;
	    landMesh.vert[i][ATTRT] = landMesh.vert[i][ATTRZ];
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
    camSetProjectionType(&cam, camPERSPECTIVE);
    camSetFrustum(&cam, M_PI / 6., 10., 10., 512, 512);
    double position[3] = {-5., -5., 20.};
    camLookFrom(&cam, position, M_PI * 0.6, angle);
	/* Run user interface. */
    render();
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