struct Shading {
    int unifDim;
    int attrDim;
    int texNum;
	int varyDim;
    void (*shadeFragment) (
        const double[], Texture&, const double[], double(&)[4]);
    void (*shadeVertex) (
        const double[], const double[], double(&)[]);
};
