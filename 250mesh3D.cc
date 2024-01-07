/*** 3D mesh builders ***/

/* Builds a mesh for a parallelepiped (box) of the given size. The attributes 
are XYZ position, ST texture, and NOP unit normal vector. The normals are 
discontinuous at the edges (flat shading, not smooth). To facilitate this, some 
vertices have equal XYZ but different NOP, for 24 vertices in all. Don't forget 
to meshFinalize when finished. */
struct Box : public Mesh<12, 24, 8> {
    Box(double left, double right, double bottom, double top, double base, 
        double lid) {
        /* Make the triangles. */
        SetTriangle(0, 0, 2, 1);
        SetTriangle(1, 0, 3, 2);
        SetTriangle(2, 4, 5, 6);
        SetTriangle(3, 4, 6, 7);
        SetTriangle(4, 8, 10, 9);
        SetTriangle(5, 8, 11, 10);
        SetTriangle(6, 12, 13, 14);
        SetTriangle(7, 12, 14, 15);
        SetTriangle(8, 16, 18, 17);
        SetTriangle(9, 16, 19, 18);
        SetTriangle(10, 20, 21, 22);
        SetTriangle(11, 20, 22, 23);
        /* Make the vertices after 0, using vertex 0 as temporary storage. */
        vert = {{{left, bottom, base, 0., 0., 0., 0., -1.},
            {right, bottom, base, 1., 0., 0., 0., -1.},
            {right, top, base, 1., 1., 0., 0., -1.},
            {left, top, base, 0., 1., 0., 0., -1.},
            {left, bottom, lid, 0., 0., 0., 0., 1.},
            {right, bottom, lid, 1., 0., 0., 0., 1.},
            {right, top, lid, 1., 1., 0., 0., 1.},
            {left, top, lid, 0., 1., 0., 0., 1.},
            {left, top, base, 0., 1., 0., 1., 0.},
            {right, top, base, 1., 1., 0., 1., 0.},
            {right, top, lid, 1., 1., 0., 1., 0.},
            {left, top, lid, 0., 1., 0., 1., 0.},
            {left, bottom, base, 0., 0., 0., -1., 0.},
            {right, bottom, base, 1., 0., 0., -1., 0.},
            {right, bottom, lid, 1., 0., 0., -1., 0.},
            {left, bottom, lid, 0., 0., 0., -1., 0.},
            {right, top, base, 1., 1., 1., 0., 0.},
            {right, bottom, base, 1., 0., 1., 0., 0.},
            {right, bottom, lid, 1., 0., 1., 0., 0.},
            {right, top, lid, 1., 1., 1., 0., 0.},
            {left, top, base, 0., 1., -1., 0., 0.},
            {left, bottom, base, 0., 0., -1., 0., 0.},
            {left, bottom, lid, 0., 0., -1., 0., 0.},
            {left, top, lid, 0., 1., -1., 0., 0.}}};
    }
};


/* Rotates a 2-dimensional vector through an angle. The output can safely alias 
the input. */
void mesh3DRotateVector(double theta, const double v[2], double *vRot) {
    double cosTheta = cos(theta);
    double sinTheta = sin(theta);
    double vRot0 = cosTheta * v[0] - sinTheta * v[1];
    vRot[1] = sinTheta * v[0] + cosTheta * v[1];
    vRot[0] = vRot0;
}

/* Rotate a curve about the Z-axis. Can be used to make a sphere, spheroid, 
capsule, circular cone, circular cylinder, box, etc. The z-values should be in 
ascending order --- or at least the first z should be less than the last. The 
first and last r-values should be 0., and no others. Probably the t-values 
should be in ascending or descending order. The sideNum parameter controls the 
fineness of the mesh. The attributes are XYZ position, ST texture, and NOP unit 
normal vector. The normals are smooth. Don't forget to meshFinalize when 
finished. */
template<size_t t, size_t v, size_t a>
void Mesh<t, v, a>::InitializeRevolution(int zNum, const double z[], const double r[],
    const double th[], int sideNum) {
    int i, j;
    /* Make the bottom triangles. */
    for (i = 0; i < sideNum; i += 1)
        SetTriangle(i, 0, i + 2, i + 1);
    /* Make the top triangles. */
    for (i = 0; i < sideNum; i += 1)
        SetTriangle(sideNum + i, v - 1, 
            v - 1 - (sideNum + 1) + i, 
            v - 1 - (sideNum + 1) + i + 1);
    /* Make the middle triangles. */
    for (j = 1; j <= zNum - 3; j += 1)
        for (i = 0; i < sideNum; i += 1) {
            SetTriangle(2 * sideNum * j + 2 * i,
                (j - 1) * (sideNum + 1) + 1 + i, 
                j * (sideNum + 1) + 1 + i + 1, 
                j * (sideNum + 1) + 1 + i);
            SetTriangle(2 * sideNum * j + 2 * i + 1,
                (j - 1) * (sideNum + 1) + 1 + i, 
                (j - 1) * (sideNum + 1) + 1 + i + 1, 
                j * (sideNum + 1) + 1 + i + 1);
        }
    /* Make the vertices, using vertex 0 as temporary storage. */
    // double *v = vert;
    for (j = 1; j <= zNum - 2; j += 1) {
        // Form the sideNum + 1 vertices in the jth layer.
        double p[3] = {z[j + 1] - z[j], 0., r[j] - r[j + 1]};
        for(int i = 0; i < 3; ++i) p[i] /= z[j + 1] - z[j] + r[j] - r[j + 1];
        double q[3] = {z[j] - z[j - 1], 0., r[j - 1] - r[j]};
        for(int i = 0; i < 3; ++i) q[i] /= z[j] - z[j - 1] + r[j - 1] - r[j];
        for(int i = 0; i < 3; ++i) p[i] += q[i];
        for(int i = 0; i < 3; ++i) p[i] *= 0.5;
        double ve[8] = {r[j], 0., z[j], 1., th[j], p[0], p[1], p[2]};
        SetVertex(j * (sideNum + 1), ve);
        ve[3] = 0.;
        SetVertex((j - 1) * (sideNum + 1) + 1, ve);
        for (i = 1; i < sideNum; i += 1) {
            mesh3DRotateVector(2 * M_PI / sideNum, ve, (double*)ve);
            ve[3] += 1. / sideNum;
            mesh3DRotateVector(2 * M_PI / sideNum, &ve[5], (double*)&ve[5]);
            SetVertex((j - 1) * (sideNum + 1) + 1 + i, ve);
        }
    }
    /* Form the top vertex. */
    SetVertex(v - 1, {0., 0., z[zNum - 1], 0., 1., 0., 0., 1.});
    /* Finally form the bottom vertex, which is set implicitly. */
    SetVertex(0, {0., 0., z[0], 0., 0., 0., 0., -1.});
}

/* Builds a mesh for a sphere, centered at the origin, of radius r. The sideNum 
and layerNum parameters control the fineness of the mesh. The attributes are 
XYZ position, ST texture, and NOP unit normal vector. The normals are smooth. 
Don't forget to meshFinalize when finished. */
template<size_t layerNum, size_t sideNum>
struct Sphere : Mesh<(layerNum-1) * sideNum * 2, (layerNum-1) * (sideNum+1) 
    + 2, 8> {
    Sphere(double r) {
        double *ts = (double *)malloc((layerNum + 1) * 3 * sizeof(double));
        if (ts == nullptr)
            return;
        else {
            double *zs = &ts[layerNum + 1];
            double *rs = &ts[2 * layerNum + 2];
            for (int i = 0; i <= layerNum; i += 1) {
                ts[i] = (double)i / layerNum;
                zs[i] = -r * cos(ts[i] * M_PI);
                rs[i] = r * sin(ts[i] * M_PI);
            }
            this->InitializeRevolution(layerNum + 1, zs, rs, ts, sideNum);
            free(ts);
        }
    }
};

/* Builds a mesh for a circular cylinder with spherical caps, centered at the 
origin, of radius r and length l > 2 * r. The sideNum and layerNum parameters 
control the fineness of the mesh. The attributes are XYZ position, ST texture, 
and NOP unit normal vector. The normals are smooth. Don't forget to meshFinalize 
when finished. */
template<size_t layerNum, size_t sideNum>
struct Capsule : public Mesh<layerNum * sideNum * 4, layerNum * (sideNum+1) * 2 
    + 2, 8> {
    Capsule(double r, double l) {
    int error, i;
    double theta;
    double *ts = (double *)malloc((2 * layerNum + 2) * 3 * sizeof(double));
    if (ts) {
        double *zs = &ts[2 * layerNum + 2];
        double *rs = &ts[4 * layerNum + 4];
        zs[0] = -l / 2.;
        rs[0] = 0.;
        ts[0] = 0.;
        for (i = 1; i <= layerNum; i += 1) {
            theta = M_PI / 2. * (3 + i / (double)layerNum);
            zs[i] = -l / 2. + r + r * sin(theta);
            rs[i] = r * cos(theta);
            ts[i] = (zs[i] + l / 2.) / l;
        }
        for (i = 0; i < layerNum; i += 1) {
            theta = M_PI / 2. * i / (double)layerNum;
            zs[layerNum + 1 + i] = l / 2. - r + r * sin(theta);
            rs[layerNum + 1 + i] = r * cos(theta);
            ts[layerNum + 1 + i] = (zs[layerNum + 1 + i] + l / 2.) / l;
        }
        zs[2 * layerNum + 1] = l / 2.;
        rs[2 * layerNum + 1] = 0.;
        ts[2 * layerNum + 1] = 1.;
        this->InitializeRevolution(2 * layerNum + 2, zs, rs, ts, sideNum);
        free(ts);
    }
}
};

/* Builds a mesh for a circular cylinder, centered at the origin, of radius r 
and length l. The sideNum parameter controls the fineness of the mesh. The 
attributes are XYZ position, ST texture, and NOP unit normal vector. The normals 
are smooth except where the side meets the ends. Don't forget to meshFinalize 
when finished. */
template<size_t sideNum>
struct Cylinder : public Mesh<4*sideNum, 4*sideNum+4, 8> {
    using Mesh<4*sideNum, 4*sideNum+4, 8>::SetTriangle;
    using Mesh<4*sideNum, 4*sideNum+4, 8>::SetVertex;

    Cylinder(double r, double l) {
    // int error = Initialize(triNum, vertNum, 3 + 2 + 3);
    // if (error != 0)
    //     return;
    double fraction, attr[3 + 2 + 3];
    /* Make the 2 * sideNum + 2 side vertices. */
    attr[7] = 0.;
    for (int i = 0; i < sideNum; i += 1) {
        fraction = (double)i / sideNum;
        attr[5] = cos(2. * M_PI * fraction);
        attr[6] = sin(2. * M_PI * fraction);
        attr[0] = r * attr[5];
        attr[1] = r * attr[6];
        attr[2] = -0.5 * l;
        attr[3] = fraction;
        attr[4] = 0.;
        SetVertex(2 * i, attr);
        attr[2] = 0.5 * l;
        attr[4] = 1.;
        SetVertex(2 * i + 1, attr);
    }
    attr[5] = cos(0.);
    attr[6] = sin(0.);
    attr[0] = r * attr[5];
    attr[1] = r * attr[6];
    attr[2] = -0.5 * l;
    attr[3] = 1.;
    attr[4] = 0.;
    SetVertex(2 * sideNum, attr);
    attr[2] = 0.5 * l;
    attr[4] = 1.;
    SetVertex(2 * sideNum + 1, attr); // same loop?  
    /* Make the sideNum + 1 top vertices. */
    attr[2] = 0.5 * l;
    attr[3] = 0.;
    attr[4] = 1.;
    attr[5] = 0.;
    //attr[6] = 0.;
    attr[7] = 1.;
    for (int i = 0; i < sideNum; i += 1) {
        attr[0] = r * cos(2. * M_PI * (double)i / sideNum);
        attr[1] = r * sin(2. * M_PI * (double)i / sideNum);
        SetVertex(2 * sideNum + 2 + i, attr);
    }
    attr[0] = 0.;
    attr[1] = 0.;
    SetVertex(2 * sideNum + 2 + sideNum, attr);
    /* Make the sideNum + 1 bottom vertices. */
    attr[2] = -0.5 * l;
    //attr[3] = 0.;
    attr[4] = 0.;
    //attr[5] = 0.;
    //attr[6] = 0.;
    attr[7] = -1.;
    for (int i = 0; i < sideNum; i += 1) {
        attr[0] = r * cos(2. * M_PI * (double)i / sideNum);
        attr[1] = r * sin(2. * M_PI * (double)i / sideNum);
        SetVertex(3 * sideNum + 3 + i, attr);
    }
    attr[0] = 0.;
    attr[1] = 0.;
    SetVertex(3 * sideNum + 3 + sideNum, attr);
    /* Make the 2 * sideNum side triangles. */
    for (int i = 0; i < sideNum; i += 1) {
        SetTriangle(2 * i, 2 * i, 2 * i + 2, 2 * i + 3);
        SetTriangle(2 * i + 1, 2 * i, 2 * i + 3, 2 * i + 1);
    }
    /* Make the sideNum top triangles. */
    for (int i = 0; i < sideNum - 1; i += 1)
        SetTriangle(2 * sideNum + i, 3 * sideNum + 2, 
            2 * sideNum + 2 + i, 2 * sideNum + 3 + i);
    SetTriangle(3 * sideNum - 1, 3 * sideNum + 2, 3 * sideNum + 1, 
        2 * sideNum + 2);
    /* Make the sideNum bottom triangles. */
    for (int i = 0; i < sideNum - 1; i += 1)
        SetTriangle(3 * sideNum + i, 4 * sideNum + 3, 
            3 * sideNum + 4 + i, 3 * sideNum + 3 + i);
    SetTriangle(4 * sideNum - 1, 4 * sideNum + 3, 3 * sideNum + 3, 
        4 * sideNum + 2);
}
};