#include <span>
#include <valarray>
#include <cstdio>
struct matrix : public std::valarray<double> {
	/*** 2 x 2 Matrices ***/
	matrix() : std::valarray<double>(4) {}

	/* Pretty-prints the given matrix, with one line of text per row of matrix. */
	void print() {
		int i, j;
		for (i = 0; i < 2; i += 1) {
			for (j = 0; j < 2; j += 1)
				printf("%f    ", (*this)[2*i+j]);
			puts("");
		}
	}

	/* Returns the determinant of the matrix m. If the determinant is 0.0, then the 
	   matrix is not invertible, and mInv is untouched. If the determinant is not 0.0, 
	   then the matrix is invertible, and its inverse is placed into mInv. The output 
	   CANNOT safely alias the input. */
	double invert() {
		return 0;
	}

	/* Multiplies a 2x2 matrix m by a 2-column v, storing the result in mTimesV. 
	   The output CANNOT safely alias the input. */
	std::valarray<double> multiply(std::span<const double> v) {
		std::valarray<double> p(2);
		p[0] = (*this)[0] * v[0] + (*this)[1] * v[1];
		p[1] = (*this)[2] * v[0] + (*this)[3] * v[1];
		return p;
	}

	/* Fills the matrix m from its two columns. The output CANNOT safely alias the 
	   input. */
	void mat22Columns(const double col0[2], const double col1[2], double m[2][2]) {

	}
};
/* The theta parameter is an angle in radians. Sets the matrix m to the 
   rotation matrix corresponding to counterclockwise rotation of the plane through 
   the angle theta. */
matrix mat22Rotation(double theta) {
	matrix m;
	m[0] = cos(theta);
	m[1] = -sin(theta);
	m[2] = sin(theta);
	m[3] = cos(theta);
	return m;
}

/* Multiplies the 3x3 matrix m by the 3x3 matrix n. The output CANNOT safely 
alias the input. */
void mat333Multiply(
        const double m[3][3], const double n[3][3], double (&mTimesN)[3][3]) {
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            mTimesN[i][j] = m[i][0]*n[0][j]+m[i][1]*n[1][j]+m[i][2]*n[2][j];
        }
    }
}

/* Multiplies the 3x3 matrix m by the 3x1 matrix v. The output CANNOT safely 
alias the input. */
void mat331Multiply(
        const double m[3][3], const double v[3], double (&mTimesV)[3]) {
    for(int i = 0; i < 3; ++i) {
        mTimesV[i] = m[i][0]*v[0]+m[i][1]*v[1]+m[i][2]*v[2];
    }
}

/* Builds a 3x3 matrix representing 2D rotation and translation in homogeneous 
coordinates. More precisely, the transformation first rotates through the angle 
theta (in radians, counterclockwise), and then translates by the vector t. */
void mat33Isometry(double theta, const double t[2], double (&isom)[3][3]) {
    isom[0][0] = cos(theta);
    isom[0][1] = -1*sin(theta);
    isom[1][0] = sin(theta);
    isom[1][1] = cos(theta);
    isom[0][2] = t[0];
    isom[1][2] = t[1];
    isom[2][0] = 0;
    isom[2][1] = 0;
    isom[2][2] = 1;
}

/* Adds the 3x3 matrix m to the 3x3 matrix n */
void mat33Add(
        const double m[3][3], const double n[3][3], double (&mPlusN)[3][3]) {
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            mPlusN[i][j] = m[i][j] + n[i][j];
        }
    }
}


/* Given a length-1 3D vector axis and an angle theta (in radians), builds the 
rotation matrix for the rotation about that axis through that angle. */
void mat33AngleAxisRotation(
        double theta, const double axis[3], double (&rot)[3][3]) {
            double u[3][3] = {{0,-1*axis[2],axis[1]},{axis[2],0,-1*axis[0]},{-1*axis[1],axis[0],0}};
            double uu[3][3];
            mat333Multiply(u, u, uu);
            for(int i = 0; i < 3; ++i) {
                for(int j = 0; j < 3; ++j) {
                    rot[i][j] = sin(theta) * u[i][j] + (1 - cos(theta)) * uu[i][j];
                    if(i == j) rot[i][j]++;
                }
            }
        }

/* Given two length-1 3D vectors u, v that are perpendicular to each other. 
Given two length-1 3D vectors a, b that are perpendicular to each other. Builds 
the rotation matrix that rotates u to a and v to b. */
void mat33BasisRotation(
        const double u[3], const double v[3], const double a[3], 
        const double b[3], double (&rot)[3][3]) {
    double w[3];
	w[0] = u[1]*v[2]-u[2]*v[1];
    w[1] = u[2]*v[0]-u[0]*v[2];
    w[2] = u[0]*v[1]-u[1]*v[0];
    // vec3Cross(u, v, w);
    double rt[3][3] = {{u[0],u[1],u[2]},{v[0],v[1],v[2]},{w[0],w[1],w[2]}};
    double c[3];
    c[0] = a[1]*b[2]-a[2]*b[1];
    c[1] = a[2]*b[0]-a[0]*b[2];
    c[2] = a[0]*b[1]-a[1]*b[0];
    // vec3Cross(a, b, c);
    double s[3][3] = {{a[0],b[0],c[0]},{a[1],b[1],c[1]},{a[2],b[2],c[2]}};
    mat333Multiply(s, rt, rot);
}

/* Computes the transpose M^T of the given matrix M. The output CANNOT safely 
alias the input. */
void mat44Transpose(const double m[4][4], double mT[4][4]) {
    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            mT[i][j] = m[j][i];
}

/* Multiplies m by n, placing the answer in mTimesN. The output CANNOT safely 
alias the input. */
void mat444Multiply(
        const double m[4][4], const double n[4][4], double (&mTimesN)[4][4]) {
    for(int i = 0; i < 4; ++i) {
        for(int j = 0; j < 4; ++j) {
            mTimesN[i][j] = m[i][0]*n[0][j]+m[i][1]*n[1][j]+m[i][2]*n[2][j]+m[i][3]*n[3][j];
        }
    }
}

/* Multiplies m by v, placing the answer in mTimesV. The output CANNOT safely 
alias the input. */
void mat441Multiply(
        const double m[4][4], const double v[4], double (&mTimesV)[4]) {
    for(int i = 0; i < 4; ++i) {
        mTimesV[i] = m[i][0]*v[0]+m[i][1]*v[1]+m[i][2]*v[2]+m[i][3]*v[3];
    }
}

/* Given a rotation and a translation, forms the 4x4 homogeneous matrix 
representing the rotation followed in time by the translation. */
void mat44Isometry(
        const double rot[3][3], const double trans[3], double (&isom)[4][4]) {
    isom[0][0] = rot[0][0]; isom[0][1] = rot[0][1]; isom[0][2] = rot[0][2]; isom[0][3] = trans[0];
    isom[1][0] = rot[1][0]; isom[1][1] = rot[1][1]; isom[1][2] = rot[1][2]; isom[1][3] = trans[1];
    isom[2][0] = rot[2][0]; isom[2][1] = rot[2][1]; isom[2][2] = rot[2][2]; isom[2][3] = trans[2];
    isom[3][0] = 0;         isom[3][1] = 0;         isom[3][2] = 0;         isom[3][3] = 1;
}

/* Sets its argument to the 4x4 zero matrix (which consists entirely of 0s). */
void mat44Zero(double (&m)[4][4]) {
    for(int i = 0; i < 16; ++i) {
        m[(int)i/4][i%4] = 0;
    }
}

/* Multiplies the transpose of the 3x3 matrix m by the 3x1 matrix v. To 
clarify, in math notation it computes M^T v. The output CANNOT safely alias the 
input. */
void mat331TransposeMultiply(
        const double m[3][3], const double v[3], double (&mTTimesV)[3]) {
    for(int i = 0; i < 3; ++i) {
        mTTimesV[i] = m[0][i]*v[0]+m[1][i]*v[1]+m[2][i]*v[2];
    }
}

/* Pretty-prints the 4x4 matrix, with one line of text per row of matrix. */
void mat44Print(const double m[4][4]) {
    int i, j;
    for (i = 0; i < 4; i += 1) {
        for (j = 0; j < 4; j += 1)
            printf("%f    ", m[i][j]);
        printf("\n");
    }
}

/* Builds a 4x4 matrix for a viewport with lower left (0, 0) and upper right 
(width, height). This matrix maps a projected viewing volume 
[-1, 1] x [-1, 1] x [-1, 1] to screen [0, w] x [0, h] x [0, 1] (each interval 
in that order). */
void mat44Viewport(double width, double height, double (&view)[4][4]) {
    view[0][0] = width/2; view[0][1] = 0; view[0][2] = 0; view[0][3] = width/2;
    view[1][0] = 0; view[1][1] = height/2; view[1][2] = 0; view[1][3] = height/2;
    view[2][0] = 0; view[2][1] = 0; view[2][2] = 0.5; view[2][3] = 0.5;
    view[3][0] = 0; view[3][1] = 0; view[3][2] = 0; view[3][3] = 1;
}

/* Inverse to mat44Viewport. */
void mat44InverseViewport(double width, double height, double (&view)[4][4]) {
    view[0][0] = 2/width; view[0][1] = 0; view[0][2] = 0; view[0][3] = -1;
    view[1][0] = 0; view[1][1] = 2/height; view[1][2] = 0; view[1][3] = -1;
    view[2][0] = 0; view[2][1] = 0; view[2][2] = 2; view[2][3] = -1;
    view[3][0] = 0; view[3][1] = 0; view[3][2] = 0; view[3][3] = 1;
}
