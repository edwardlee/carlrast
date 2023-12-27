#include <random>
constexpr int LANDSIZE = 40;
using namespace std;

/* A landscape is simply a square array of doubles, with each one giving an 
elevation. This file contains functions for generating landscapes.
To turn the landscape into a mesh, use the appropriate 3D mesh initializer 
functions. */
struct Land {
	array<array<double, LANDSIZE>, LANDSIZE> data{};
	mt19937 gen;
	Land() : gen(random_device{}()) {}

/* Given a line y = m x + b across the landscape (with the x-axis pointing east 
and the y-axis pointing north), raises points north of the line and lowers 
points south of it (or vice-versa). */
void FaultEastWest(
        double m, double b, double raisingNorth) {
    for (int j = 0; j < LANDSIZE; ++j)
        for (int i = 0; i < LANDSIZE; ++i)
            data[i][j] += raisingNorth * bit_cast<char>(j <=> m * i + b);
}

/* Given a line x = m y + b across the landscape (with the x-axis pointing east 
and the y-axis pointing north), raises points east of the line and lower points 
west of it (or vice-versa). */
void FaultNorthSouth(double m, double b, double raisingEast) {
    for (int i = 0; i < LANDSIZE; ++i)
        for (int j = 0; j < LANDSIZE; ++j)
			data[i][j] += raisingEast * bit_cast<char>(i <=> m * j + b);
}

/* Randomly chooses a vertical fault and slips the landscape up and down on the 
two sides of that fault. */
void FaultRandomly(double magnitude) {
	double b, m = uniform_real_distribution(-1., 1.)(gen);
	int sign = 2 * bernoulli_distribution(0.5)(gen) - 1;
	if (m > 0)
		b = uniform_real_distribution(-m, 1.)(gen) * LANDSIZE;
	else
		b = uniform_real_distribution(0., 1. - m)(gen) * LANDSIZE;
	double raising = magnitude * uniform_real_distribution(0.5, 1.5)(gen) * sign;
	if (bernoulli_distribution(0.5)(gen))
		FaultEastWest(m, b, raising);
	else
		FaultNorthSouth(m, b, raising);
}

/* Blurs each non-border elevation with the eight elevations around it. */
void Blur() {
	array<array<double, LANDSIZE>, LANDSIZE> copy = data;
	for (int i = 1; i < LANDSIZE - 1; ++i)
		for (int j = 1; j < LANDSIZE - 1; ++j) {
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
    for (int i = 0; i < LANDSIZE; ++i)
        for (int j = 0; j < LANDSIZE; ++j) {
            distSq = (i - x) * (i - x) + (j - y) * (j - y);
            data[i][j] += raising * exp(scalar * distSq);
        }
}
};