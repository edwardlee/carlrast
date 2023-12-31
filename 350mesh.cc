#include <array>
using std::array;


/*** Creating and destroying ***/

/* Feel free to read the struct's members, but don't write them, except through 
the accessors below such as SetTriangle, SetVertex. */
template <size_t t, size_t v, size_t a=8>
struct Mesh {
	array<int[3], t> tri;		/* t * 3 ints */
	array<double[a], v> vert;	/* v * a doubles */


/* Sets the trith triangle to have vertex indices i, j, k. */
void SetTriangle(int n, int i, int j, int k) {
	if (0 <= n && n < t) {
		tri[n][0] = i;
		tri[n][1] = j;
		tri[n][2] = k;
	}
}

/* Sets the vertth vertex to have attributes attr. */
void SetVertex(int n, double (&&attr)[a]) {
	ranges::copy(attr, vert[n]);
}

/*** Rendering ***/
void lerp
	(int varyDim, const double va[], const double vb[], double (&vn)[]) {
	double T = (va[2] + va[3])/(va[2] + va[3] - vb[2] - vb[3]);
	for(int i = 0; i < varyDim; ++i)
		vn[i] = std::lerp(va[i], vb[i], T);
}

/* Renders the . If the  and the shading have differing values for 
a, then prints an error message and does not render anything. */
void Render(
        Depth buf, const double viewport[4][4], const shaShading *sha, 
	const double (&unif)[], Texture &tex) {
    if(a != sha->attrDim) {
        fprintf(stderr, "expected as to match\n");
        return;
    }
    for(int i = 0; i < t; ++i) {
        double va[sha->varyDim], vb[sha->varyDim], vc[sha->varyDim];
        sha->shadeVertex(sha->unifDim, unif, sha->attrDim, vert[tri[i][0]], sha->varyDim, va);
        sha->shadeVertex(sha->unifDim, unif, sha->attrDim, vert[tri[i][1]], sha->varyDim, vb);
        sha->shadeVertex(sha->unifDim, unif, sha->attrDim, vert[tri[i][2]], sha->varyDim, vc);
		bool ca = va[3] <= 0. || va[3] < -va[2],
			 cb = vb[3] <= 0. || vb[3] < -vb[2],
			 cc = vc[3] <= 0. || vc[3] < -vc[2];
		double x[sha->varyDim], y[sha->varyDim];
		auto render = [&](const double (&va)[], const double (&vb)[], const double (&vc)[]) {
			double vv[4], vva[sha->varyDim], vvb[sha->varyDim], vvc[sha->varyDim];
			mat441Multiply(viewport, va, vv);
			for(int i = 0; i < 3; ++i) vva[i] = vv[i] / vv[3];
			for(int i = 4; i < sha->varyDim; ++i) vva[i] = va[i]/vv[3];
			mat441Multiply(viewport, vb, vv);
			for(int i = 0; i < 3; ++i) vvb[i] = vv[i] / vv[3];
			for(int i = 4; i < sha->varyDim; ++i) vvb[i] = vb[i]/vv[3];
			mat441Multiply(viewport, vc, vv);
			for(int i = 0; i < 3; ++i) vvc[i] = vv[i] / vv[3];
			for(int i = 4; i < sha->varyDim; ++i) vvc[i] = vc[i]/vv[3];
			triRender(sha, buf, unif, tex, vva, vvb, vvc);
		};
		if(ca) {
			if(cb) {
				if(!cc) {
					lerp(sha->varyDim, va, vc, x);
					lerp(sha->varyDim, vb, vc, y);
					render(y, vc, x);
				}
			} else if(cc) {
				lerp(sha->varyDim, va, vb, x);
				lerp(sha->varyDim, vb, vc, y);
				render(x, vb, y);
			} else {
				lerp(sha->varyDim, va, vb, x);
				lerp(sha->varyDim, va, vc, y);
				render(x, vb, y);
				render(y, vb, vc);
			}
		} else if(cb) {
			if(cc) {
				lerp(sha->varyDim, va, vb, x);
				lerp(sha->varyDim, va, vc, y);
				render(y, va, x);
			} else {
				lerp(sha->varyDim, va, vb, x);
				lerp(sha->varyDim, vb, vc, y);
				render(y, vc, x);
				render(x, vc, va);
			}
		} else if(cc) {
			lerp(sha->varyDim, va, vc, x);
			lerp(sha->varyDim, vb, vc, y);
			render(x, va, y);
			render(y, va, vb);
		} else {
			render(va, vb, vc);
		}
	}
}

int InitializeRectangle(double, double, double, double);
int InitializeEllipse(double, double, double, double, int);
void InitializeBox(double, double, double, double, double, double);
void InitializeRevolution(int, const double[], const double[], const double[], int);
void InitializeSphere(double, int, int);
void InitializeCylinder(double, double, int);
void InitializeCapsule(double, double, int, int);
};
