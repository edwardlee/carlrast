#include <algorithm>
#include "040pixel.h"
using std::min, std::max;

template<Shading const &sha>
void triRender(
        Depth &buf, const double (&unif)[], 
        Texture &tex, const double (&a)[], const double (&b)[], 
        const double (&c)[]) {
    // 28.4 fixed-point coordinates
    const int Y1 = (16. * a[1]);
    const int Y2 = (16. * b[1]);
    const int Y3 = (16 * c[1]);

    // Bounding rectangle
    int miny = max((min({Y1, Y2, Y3}) + 0xF) >> 4, 0);
    int maxy = min((max({Y1, Y2, Y3}) + 0xF) >> 4, buf.height);

    const int X1 = (16 * a[0]);
    const int X2 = (16. * b[0]);
    const int X3 = (16 * c[0]);

    int minx = max((min({X1, X2, X3}) + 0xF) >> 4, 0);
    int maxx = min((max({X1, X2, X3}) + 0xF) >> 4, buf.width);

    // Deltas
    int DX12 = X2 - X1;
    int DX23 = X3 - X2;
    int DX31 = X1 - X3;

    int DY12 = Y2 - Y1;
    int DY23 = Y3 - Y2;
    int DY31 = Y1 - Y3;

    // Half-edge constants
    int C1 = DY12 * X1 - DX12 * Y1;
    int C2 = DY23 * X2 - DX23 * Y2;
    int C3 = DY31 * X3 - DX31 * Y3;

    // Correct for fill convention
    if(DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
    if(DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
    if(DY31 < 0 || (DY31 == 0 && DX31 > 0)) C3++;

    int CY1 = C1 + DX12 * (miny << 4) - DY12 * (minx << 4);
    int CY2 = C2 + DX23 * (miny << 4) - DY23 * (minx << 4);
    int CY3 = C3 + DX31 * (miny << 4) - DY31 * (minx << 4);

	double den = 1. / (DX23 * DY31 - DY23 * DX31);

    // Fixed-point deltas
    DX12 <<= 4;
    DX23 <<= 4;
    DX31 <<= 4;

    DY12 <<= 4;
    DY23 <<= 4;
    DY31 <<= 4;
    int index = 3 * (miny * buf.width + minx);
    for(int y = miny; y < maxy; ++y) {
        int CX1 = CY1;
        int CX2 = CY2;
        int CX3 = CY3;
        int ind = index;
        for(int x = minx; x < maxx; ++x) {
            if(CX1 > 0 && CX2 > 0 && CX3 > 0) {
				double attr[sha.varyDim];
				for(int i = 0; i < sha.varyDim; ++i)
					attr[i] = (a[i]*CX2 + b[i]*CX3 + c[i]*CX1)*den;
				double rgbd[4];
				sha.shadeFragment(unif, tex, attr, rgbd);
                if(rgbd[3] < buf[y][x]) {
                    pixSetRGB(ind, rgbd);
                    buf[y][x] = rgbd[3];
                }
            }
            CX1 -= DY12;
            CX2 -= DY23;
            CX3 -= DY31;
            ind += 3;
        }
        CY1 += DX12;
        CY2 += DX23;
        CY3 += DX31;
        index += 3 * buf.width;
    }
}
