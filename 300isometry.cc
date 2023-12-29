/* Describes an isometry as a rotation followed by a translation. Can be used 
to describe the position and orientation of a rigid body. If the position is 
the translation, and the columns of the rotation are the local coordinate axes 
in global coordinates, then the isometry takes local coordinates to global. */

struct Isometry {
	double translation[3];
	double rotation[3][3];

/* Sets the rotation. */
void SetRotation(const double rot[3][3]) {
	std::copy_n((double *)rot, 9, (double *)(rotation));
}

/* Sets the translation. */
void SetTranslation(const double transl[3]) {
	std::copy_n(transl, 3, translation);
}

/* Applies the rotation and translation to a point. The output CANNOT safely 
alias the input. */
void TransformPoint(const double p[3], double (&isoP)[3]) {
	mat331Multiply(rotation, p, isoP);
	for(int i = 0; i < 3; ++i)
		isoP[i] += translation[i];
}

/* Applies the inverse of the isometry to a point. If you transform a point and 
then untransform the result, then you recover the original point. Similarly, if 
you untransform a point and then transform the result, then you recover the 
original point. The output CANNOT safely alias the input. */
void UntransformPoint(const double isoP[3], double (&p)[3]) {
	double untrans[3];
	untrans[0] = isoP[0]-translation[0],
	untrans[1] = isoP[1]-translation[1],
	untrans[2] = isoP[2]-translation[2];
	mat331TransposeMultiply(rotation, untrans, p);
}


/* Applies the rotation to a direction vector (typically unit). The output 
CANNOT safely alias the input. */
void RotateDirection(
        const double d[3], double (&rotD)[3]) {
	mat331Multiply(rotation, d, rotD);
}


/* Applies the inverse rotation to a direction vector (typically unit). The 
output CANNOT safely alias the input. */
void UnrotateDirection(
        const double rotD[3], double (&d)[3]) {
	mat331TransposeMultiply(rotation, rotD, d);
}


/* Fills homog with the homogeneous version of the isometry. */
void GetHomogeneous(double (&homog)[4][4]) {
	mat44Isometry(rotation, translation, homog);
}

/* Fills homog with the homogeneous version of the inverse isometry. That is, 
the product of this matrix and the one from isoGetHomogeneous is the identity 
matrix. */
void GetInverseHomogeneous(double (&homogInv)[4][4]) {
	double r[3][3], t[3];
	for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            r[j][i] = rotation[i][j];
	mat331Multiply(r, translation, t);
	t[0] *= -1;
	t[1] *= -1;
	t[2] *= -1; // can I combine and/or lambda
	mat44Isometry(r, t, homogInv);
}
};