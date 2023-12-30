#include <algorithm>
#include <span>
using namespace std;

struct Depth : span<double> {
	size_t width, height;

	Depth(size_t width, size_t height)
		: span(new double[width * height], width* height), width(width), height(height) {}

	/* Sets every depth-value to the given depth. Typically you use this function
	at the start of each frame, passing a large positive value for depth. */
	void Clear(double depth) { ranges::fill(*this, depth); }

	span<double> operator[](size_t j) { return this->subspan(j * width); }
};