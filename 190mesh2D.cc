


/*** 2D  builders ***/

/* Initializes a  to two triangles forming a rectangle of the given sides. 
The four attributes are X, Y, S, T. Do not call Initialize separately; it 
is called inside this function. Don't forget to call Finalize when done. */
template <size_t t, size_t v, size_t a>
int Mesh<t, v, a>::InitializeRectangle(
        double left, double right, double bottom, double top) {
    SetTriangle(0, 0, 1, 2);
    SetTriangle(1, 0, 2, 3);
    SetVertex(0, {left, bottom, 0., 0.});
    SetVertex(1, {right, bottom, 1., 0.});
    SetVertex(2, {right, top, 1., 1.});
    SetVertex(3, {left, top, 0., 1.});
    return 0;
}

/* Initializes a  to sideNum triangles forming an ellipse of the given 
center (x, y) and radii rx, ry. The four attributes are X, Y, S, T. Do not call 
Initialize separately; it is called inside this function. Don't forget to 
call Finalize when done. */
template <size_t t, size_t v, size_t a>
int Mesh<t, v, a>::InitializeEllipse(
        double x, double y, double rx, double ry, int sideNum) {
    int i, error;
    double theta, cosTheta, sinTheta;
    //error = Initialize(sideNum, sideNum + 1, 2 + 2);
    SetVertex(0, {x, y, 0.5, 0.5});
    for (i = 0; i < sideNum; i += 1) {
        SetTriangle(i, i+1, 0, (i + 1) % sideNum + 1);
        theta = i * 2.0 * M_PI / sideNum;
        cosTheta = cos(theta);
        sinTheta = sin(theta);
        SetVertex(i + 1, {x + rx * cosTheta, y + ry * sinTheta, 
            0.5 * cosTheta + 0.5, 0.5 * sinTheta + 0.5});
    }
    return 0;
}

/* Other ideas: capsule, annulus, arbitrary simple polygon, ... */


