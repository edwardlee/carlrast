struct shaShading {
    int unifDim;
    int attrDim;
    int texNum;
	int varyDim;
    void (*shadeFragment) (
        int, const double[], int, const texTexture*[], 
        int, const double[], double(&)[4]);
    void (*shadeVertex) (
        int, const double[], int, const double[], 
        int, double(&)[]);
};
