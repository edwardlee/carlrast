#include <iostream>
#include <span>
#include <algorithm>
using namespace std;
/*** Public: For header file ***/

/* These are constants that are set at compile time. For example, whenever the 
compiler sees 'texLINEAR', it will substitute '0'. Let me emphasize: texLINEAR 
is not a variable. It does not occupy any memory in your running program, and 
your program cannot change its value. We use such constants to avoid having 
'magic numbers' sprinkled throughout our code. */
#define texLINEAR 0
#define texNEAREST 1
#define texREPEAT 2
#define texCLIP 3


/*** Private ***/

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STBI_FAILURE_USERMSG

struct Texture {
    int width, height;  /* do not have to be powers of 2 */
    int texelDim;       /* e.g. 3 for RGB textures */
    char filtering;      /* texLINEAR or texNEAREST */
    char topBottom;      /* texREPEAT or texCLIP */
    char leftRight;      /* texREPEAT or texCLIP */
    span<double> data;       /* width * height * texelDim doubles, row-major order */





/*** Public: Basics ***/

/* Sets all texels within the texture. Assumes that the texture has already 
been initialized. Assumes that texel has the same texel dimension as the 
texture. */
template<size_t D>
void ClearTexels(double (&texel)[D]) {
    for (int i = 0; i < data.size(); i += texelDim)
        ranges::copy(texel, &data[i]);
}

/* Initializes a texTexture struct to a given width and height and a solid 
color. The width and height do not have to be powers of 2. Returns 0 if no 
error occurred. The user must remember to call texFinalize when finished with 
the texture. */
template<size_t D>
Texture(int w, int h, int d, double (&&texel)[D]) {
    width = w;
    height = h;
    texelDim = d;
    data = span(new double[w * h * d], w * h * d);
    ClearTexels(texel);
}

/* Initializes a texTexture struct by loading an image from a file. Many image 
types are supported (using the public-domain STB Image library). The width and 
height do not have to be powers of 2. Returns 0 if no error occurred. The user 
must remember to call texFinalize when finished with the texture. */
/* WARNING: Currently there is a weird behavior, in which some image files show 
up with their rows and columns switched, so that their width and height are 
flipped. If that's happening with your image, then use a different image. */
Texture(const char *path) {
    /* Use the STB image library to load the file as unsigned chars. */
    stbi_set_flip_vertically_on_load(true);
    unsigned char *rawData = stbi_load(path, &width, &height, &texelDim, 0);
    if (!rawData) {
        cerr << "error: texInitializeFile: failed to load image " << path << endl
            << "    with STB Image reason: " << stbi_failure_reason() << endl;
        return;
    }
    data = span(new double[width * height * texelDim], width * height * texelDim);
    for (int i = 0; double &d : data)
        d = rawData[i++] / 255.;
    stbi_image_free(rawData);
}

/*
For image files with their rows and columns switched, we must use this code 
instead. I'm not sure how to detect this case. So only the other case is 
handled in the code above.
    rawData = stbi_load(path, &(height), &(width), &(texelDim), 0);
    ...
    newInd = texelDim * (x + width * y);
    oldInd = texelDim * (height * (width - x + 1) - y);
*/

/* Gets a single texel within the texture. Assumes that texel has the same texel 
dimension as the texture. Texel (s, t) = (0, 0) is in the lower left corner, 
texel (width - 1, 0) is in the lower right corner, etc. */
void GetTexel(const int s, int t, double (&texel)[]) {
    int k;
    for (k = 0; k < texelDim; k += 1)
        texel[k] = data[(s + width * t) * texelDim + k];
}

/* Sets a single texel within the texture. For details, see texGetTexel. */
template<size_t D>
void SetTexel(int x, int y, const double (&texel)[D]) {
    int index = texelDim * (x + width * y);
    if (index < data.size())
        ranges::copy(texel, begin(data) + index);
}



/*** Public: Higher-level sampling ***/

/* Samples from the texture, taking into account wrapping and filtering. The s 
and t parameters are texture coordinates. The texture itself is assumed to have 
texture coordinates [0, 1] x [0, 1], with (0, 0) in the lower left corner, (1, 
0) in the lower right corner, etc. Assumes that the texture has already been 
initialized. Assumes that sample has been allocated with (at least) texelDim 
doubles. Places the sampled texel into sample. */
void Sample(double s, double t, double (&sample)[]) {
    /* Handle clipping vs. repeating. */
    if (leftRight == texREPEAT)
        s -= floor(s);
    else {
        [[unlikely]] if (s < 0.)
            s = 0.;
        else if (s > 1.)
            s = 1.;
    }
    if (topBottom == texREPEAT)
        t -= floor(t);
    else {
        if (t < 0.)
            t = 0.;
        else if (t > 1.)
            t = 1.;
    }
    /* Scale to image space. */
    double u = s * (width - 1);
    double v = t * (height - 1);
    /* Handle nearest-neighbor vs. linear filtering. */
    if (filtering == texNEAREST)
        GetTexel(round(u), round(v), sample);
    else {
        // 0,0 , 0,1 , 1,0 , 1,1	
        double ff[texelDim], fc[texelDim], cf[texelDim], cc[texelDim];
        GetTexel(u, v, ff); 
        GetTexel(u, ceil(v), fc);	
        GetTexel(ceil(u), v, cf);	
        GetTexel(ceil(u), ceil(v), cc);	
        // m for mantissa	
        double um = u - floor(u);	
        double vm = v - floor(v);	
        // 1 - is because closer points have greater weight 
        for(double &n : ff) n *= (1-um) * (1-vm);
        for(double &n : fc) n *= (1-um) * vm;
        for(double &n : cf) n *= um * (1-vm);
        for(double &n : cc) n *= um * vm;
        // Take a weighted average for optimal smoothness
        for (int k = 0; k < texelDim; ++k)	
            sample[k] = ff[k] + fc[k] + cf[k] + cc[k];
    }
}


};
