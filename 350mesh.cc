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
void SetVertex(int n, const double (&attr)[a]) {
	std::ranges::copy(attr, vert[n]);
}

/*** Rendering ***/
void lerp
	(int varyDim, const double va[], const double vb[], double (&vn)[]) {
	double T = (va[2] + va[3])/(va[2] + va[3] - vb[2] - vb[3]);
	for(int i = 0; i < varyDim; ++i)
		vn[i] = std::lerp(va[i], vb[i], T);
}

template<Shading const &sha>
void Render(
        Depth &buf, const double viewport[4][4], 
	const double (&unif)[], Texture &tex) {
    for(int (&T)[3] : tri) {
        double va[sha.varyDim], vb[sha.varyDim], vc[sha.varyDim];
        sha.shadeVertex(unif, vert[T[0]], va);
        sha.shadeVertex(unif, vert[T[1]], vb);
        sha.shadeVertex(unif, vert[T[2]], vc);
		bool ca = va[3] <= 0. || va[3] < -va[2],
			 cb = vb[3] <= 0. || vb[3] < -vb[2],
			 cc = vc[3] <= 0. || vc[3] < -vc[2];
		double x[sha.varyDim], y[sha.varyDim];
		auto render = [&](const double (&va)[], const double (&vb)[], const double (&vc)[]) {
			double vva[sha.varyDim], vvb[sha.varyDim], vvc[sha.varyDim];
			mat441Multiply(viewport, va, vva);
			for(int i = 0; i < 3; ++i) vva[i] /= vva[3];
			for(int i = 4; i < sha.varyDim; ++i) vva[i] = va[i]/vva[3];
			mat441Multiply(viewport, vb, vvb);
			for(int i = 0; i < 3; ++i) vvb[i] /= vvb[3];
			for(int i = 4; i < sha.varyDim; ++i) vvb[i] = vb[i]/vvb[3];
			mat441Multiply(viewport, vc, vvc);
			for(int i = 0; i < 3; ++i) vvc[i] /= vvc[3];
			for(int i = 4; i < sha.varyDim; ++i) vvc[i] = vc[i]/vvc[3];
			triRender<sha>(buf, unif, tex, vva, vvb, vvc);
		};
		if(ca) {
			if(cb) {
				if(!cc) {
					lerp(sha.varyDim, va, vc, x);
					lerp(sha.varyDim, vb, vc, y);
					render(y, vc, x);
				}
			} else if(cc) {
				lerp(sha.varyDim, va, vb, x);
				lerp(sha.varyDim, vb, vc, y);
				render(x, vb, y);
			} else {
				lerp(sha.varyDim, va, vb, x);
				lerp(sha.varyDim, va, vc, y);
				render(x, vb, y);
				render(y, vb, vc);
			}
		} else if(cb) {
			if(cc) {
				lerp(sha.varyDim, va, vb, x);
				lerp(sha.varyDim, va, vc, y);
				render(y, va, x);
			} else {
				lerp(sha.varyDim, va, vb, x);
				lerp(sha.varyDim, vb, vc, y);
				render(y, vc, x);
				render(x, vc, va);
			}
		} else if(cc) {
			lerp(sha.varyDim, va, vc, x);
			lerp(sha.varyDim, vb, vc, y);
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
