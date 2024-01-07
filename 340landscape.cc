#include <random>
using namespace std;

/* A landscape is simply a square array of doubles, with each one giving an 
elevation. This file contains functions for generating landscapes. */
template<size_t size>
struct Land : public Mesh<2 * (size - 1) * (size - 1), size * size, 8> {
	array<array<double, size>, size> data{};
	mt19937 gen;
	Land() : gen(random_device{}()) {}

/* Given a line y = m x + b across the landscape (with the x-axis pointing east 
and the y-axis pointing north), raises points north of the line and lowers 
points south of it (or vice-versa). */
void FaultEastWest(
        double m, double b, double raisingNorth) {
    for (int j = 0; j < size; ++j)
        for (int i = 0; i < size; ++i)
            data[i][j] += raisingNorth * bit_cast<char>(j <=> m * i + b);
}

/* Given a line x = m y + b across the landscape (with the x-axis pointing east 
and the y-axis pointing north), raises points east of the line and lower points 
west of it (or vice-versa). */
void FaultNorthSouth(double m, double b, double raisingEast) {
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
			data[i][j] += raisingEast * bit_cast<char>(i <=> m * j + b);
}

/* Randomly chooses a vertical fault and slips the landscape up and down on the 
two sides of that fault. */
void FaultRandomly(double magnitude) {
	double b, m = uniform_real_distribution(-1., 1.)(gen);
	int sign = 2 * bernoulli_distribution(0.5)(gen) - 1;
	if (m > 0)
		b = uniform_real_distribution(-m, 1.)(gen) * size;
	else
		b = uniform_real_distribution(0., 1. - m)(gen) * size;
	double raising = magnitude * uniform_real_distribution(0.5, 1.5)(gen) * sign;
	if (bernoulli_distribution(0.5)(gen))
		FaultEastWest(m, b, raising);
	else
		FaultNorthSouth(m, b, raising);
}

/* Blurs each non-border elevation with the eight elevations around it. */
void Blur() {
	array<array<double, size>, size> copy = data;
	for (int i = 1; i < size - 1; ++i)
		for (int j = 1; j < size - 1; ++j) {
			copy[i][j] = 
				(data[i][j] + 
				data[i + 1][j] + 
				data[i - 1][j] + 
				data[i][j + 1] + 
				data[i][j - 1] + 
				data[i + 1][j + 1] + 
				data[i + 1][j - 1] + 
				data[i - 1][j + 1] + 
				data[i - 1][j - 1]) / 9.;
		}
	data = move(copy);
}

/* Forms a Gaussian hill or valley at (x, y), with width controlled by stddev 
and height/depth controlled by raising. */
void Bump(
        int x, int y, double stddev, double raising) {
    double scalar, distSq;
    scalar = -0.5 / (stddev * stddev);
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j) {
            distSq = (i - x) * (i - x) + (j - y) * (j - y);
            data[i][j] += raising * exp(scalar * distSq);
        }
}

/* Assumes that attributes 0, 1, 2 are XYZ. Assumes that the vertices of the 
triangle are in counter-clockwise order when viewed from 'outside' the 
triangle. Computes the outward-pointing unit normal vector for the triangle. 
The output CANNOT safely alias the input. */
void mesh3DTrueNormal(
        const double a[], const double b[], const double c[], 
        double (&normal)[3]) {
    double bMinusA[3] = {b[0]-a[0],b[1]-a[1],b[2]-a[2]}, cMinusA[3] = {c[0]-a[0],c[1]-a[1],c[2]-a[2]};
	normal[0] = bMinusA[1]*cMinusA[2]-bMinusA[2]*cMinusA[1];
    normal[1] = bMinusA[2]*cMinusA[0]-bMinusA[0]*cMinusA[2];
    normal[2] = bMinusA[0]*cMinusA[1]-bMinusA[1]*cMinusA[0];
	double l = hypot(normal[0], normal[1], normal[2]);
	normal[0] /= l;
	normal[1] /= l;
	normal[2] /= l;
}

/* Assumes that attributes 0, 1, 2 are XYZ. Sets attributes n, n + 1, n + 2 to 
flat-shaded normals. If a vertex belongs to more than triangle, then some 
unspecified triangle's normal wins. */
void mesh3DFlatNormals(int n) {
    double *a_, *b, *c, normal[3];
    for (int *p : this->tri) {
        a_ = this->vert[p[0]];
        b = this->vert[p[1]];
        c = this->vert[p[2]];
        mesh3DTrueNormal(a_, b, c, normal);
		ranges::copy(normal, &a_[n]);
		ranges::copy(normal, &b[n]);
		ranges::copy(normal, &c[n]);
    }
}

/* Assumes that attributes 0, 1, 2 are XYZ. Sets attributes n, n + 1, n + 2 to 
smooth-shaded normals. Does not do anything special to handle multiple vertices 
with the same coordinates. */
void mesh3DSmoothNormals(int n) {
    double *a, *b, *c, normal[3];
    /* For each triangle, add onto the normal at each of its vertices. */
    for (int *p : this->tri) {
        a = this->vert[p[0]];
        b = this->vert[p[1]];
        c = this->vert[p[2]];
        mesh3DTrueNormal(a, b, c, normal);
        for(int i = 0; i < 3; ++i) a[n+i] += normal[i];
        for(int i = 0; i < 3; ++i) b[n+i] += normal[i];
        for(int i = 0; i < 3; ++i) c[n+i] += normal[i];
    }
    /* Normalize the normals. */
    for (double *a : this->vert) {
        double l = hypot(a[n], a[n+1], a[n+2]);
        for(int i = 0; i < 3; ++i) a[n+i] /= l;
    }
}

/* Builds a non-closed 'landscape' mesh based on a grid of Z-values. There are 
size * size Z-values, which arrive in the data parameter. The mesh is made of 
(size - 1) * (size - 1) squares, each made of two triangles. The spacing 
parameter controls the spacing of the X- and Y-coordinates of the vertices. The 
attributes are XYZ position, ST texture, and NOP unit normal vector. To understand the exact 
layout of the data, try this example code:
double zs[3][3] = {
    {10., 9., 7.}, 
    {6., 5., 3.}, 
    {4., 3., -1.}}; */
void Build(double spacing=1.) {
    int i, j;
    int a, b, c, d;
    double diffSWNE, diffSENW;
        /* Build the vertices with normals set to 0. */
        for (i = 0; i < size; ++i)
            for (j = 0; j < size; ++j) {
                ranges::copy((double(&&)[8]){i * spacing, j * spacing, data[i][j], 
                    (double)i, (double)j}, this->vert[i * size + j]);
            }
        /* Build the triangles. */
        for (i = 0; i < size - 1; ++i)
            for (j = 0; j < size - 1; ++j) {
                int index = 2 * (i * (size - 1) + j);
                a = i * size + j;
                b = (i + 1) * size + j;
                c = (i + 1) * size + (j + 1);
                d = i * size + (j + 1);
                diffSWNE = abs(this->vert[a][2] - this->vert[c][2]);
                diffSENW = abs(this->vert[b][2] - this->vert[d][2]);
                if (diffSENW < diffSWNE) {
                    this->SetTriangle(index, d, a, b);
                    this->SetTriangle(index + 1, b, c, d);
                } else {
                    this->SetTriangle(index, a, b, c);
                    this->SetTriangle(index + 1, a, c, d);
                }
            }
        /* Set the normals. */
        mesh3DFlatNormals(5);
}

};