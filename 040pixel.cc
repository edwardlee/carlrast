/*
On Ubuntu, compile with...
    g++ -c 040pixel.cc -std=c++20 -Ofast -march=native
...and then link with a main program by for example...
    cc main.c 040pixel.o -lglfw -lGL -lm -ldl
*/

/*
There is some missing error checking in this code.

I stopped revising midway through a revision that would allow width and height 
to be other than powers of 2. The idea was to make the OpenGL texture have 
width and height equal to pixPowerOfTwoCeil of the user width and height.
*/



/*** Private: infrastructure ***/

#include <iostream>
#include <span>
#include <GL/gl3w.h>
#include <GLFW/glfw3.h>

// Global variables.
GLFWwindow *pixWindow;
int pixWidth, pixHeight; // as initialized
GLuint pixTexture;
std::span<float> pixPixels;
bool pixNeedsRedisplay = true;
GLuint pixAttrBuffer, pixTriBuffer;
GLuint pixProgram;
GLint pixUnifLoc, pixAttrLoc;
double pixOldTime, pixNewTime;
void (*pixUserKeyDownHandler)(int);
void (*pixUserKeyUpHandler)(int);
void (*pixUserKeyRepeatHandler)(int);
void (*pixUserTimeStepHandler)(double, double);

GLuint pixBuildVertexFragmentProgram(const GLchar *vertexCode, 
    const GLchar *fragmentCode) {
    GLuint vertexShader, fragmentShader, program;
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexCode, 0);
    // !!check for errors from glCompileShader
    glCompileShader(vertexShader);
    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentCode, 0);
    // !!check for errors from glCompileShader
    glCompileShader(fragmentShader);
    program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    glLinkProgram(program);
    return program;
}



/*** Private: GLFW handlers ***/

void pixHandleError(int error, const char *description) {
    std::cerr << "GLFW code " << error << ", message...\n"
        << description << std::endl;
}

// key is GLFW_KEY_A, GLFW_KEY_B, etc. Or GLFW_KEY_UNKNOWN.
// scancode is a platform-dependent scan code for the key.
// action is GLFW_PRESS, GLFW_RELEASE, or GLFW_REPEAT.
// mods has bitwise masks for...
//     GLFW_MOD_SHIFT, GLFW_MOD_CONTROL, GLFW_MOD_ALT, or GLFW_MOD_SUPER
// ...which on macOS mean shift, control, option, command.
void pixHandleKey(GLFWwindow *window, int key, int scancode, int action,
        int mods) {
    if (action == GLFW_PRESS && pixUserKeyDownHandler)
        pixUserKeyDownHandler(key);
    if (action == GLFW_RELEASE && pixUserKeyUpHandler)
        pixUserKeyUpHandler(key);
    else if (pixUserKeyRepeatHandler)
        pixUserKeyRepeatHandler(key);
}

/*** Private: initializer ***/

// Initialize the user interface and OpenGL.
int pixInitGLFWGL3W(std::string_view name) {
    glfwSetErrorCallback(pixHandleError);
    // !!error if glfwInit() returns GL_FALSE
    glfwInit();
    glfwWindowHint(GLFW_RESIZABLE, 0);
    // This code is commented-out because 2.1 is the default.
    // glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    // glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    // glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    pixWindow = glfwCreateWindow(pixWidth, pixHeight, data(name), nullptr, nullptr);
    glfwMakeContextCurrent(pixWindow);
    gl3wInit();
    glfwSetKeyCallback(pixWindow, pixHandleKey);
    // The following code used to be...
    //glViewport(0, 0, pixWidth, pixHeight);
    int width, height;
    glfwGetFramebufferSize(pixWindow, &width, &height);
    glViewport(0, 0, width, height);
    return 0;
}

// Create the texture.
int pixInitTexture() {
    // !!error if malloc returns NULL
    pixPixels = std::span(new float[3 * pixWidth * pixHeight], 3 * pixWidth * pixHeight);
    // If we were using OpenGL 4.5, we might do this.
    glCreateTextures(GL_TEXTURE_2D, 1, &pixTexture);
    glTextureStorage2D(pixTexture, 1, GL_RGB32F, pixWidth, pixHeight);
    glBindTextureUnit(0, pixTexture);
    return glGetError();
}

/*
// Vertex shader.
uniform mat4 cameraMatrix;
attribute vec4 attribs;
varying vec2 texCoords;
void main(){
    gl_Position = cameraMatrix * vec4(attribs[0], attribs[1], 0., 1.);
    texCoords = vec2(attribs[2], attribs[3]);
}
// Fragment shader.
varying vec2 texCoords;
uniform sampler2D texture0;
void main(){
    gl_FragColor = texture2D(texture0, texCoords);
}
*/
// Build the shader program and leave it bound.
int pixInitShaders() {
    GLchar vertexCode[] = "uniform mat4 cameraMatrix;attribute vec4 attribs;varying vec2 texCoords;void main(){gl_Position = cameraMatrix * vec4(attribs[0], attribs[1], 0., 1.);texCoords = vec2(attribs[2], attribs[3]);}";
    GLchar fragmentCode[] = "varying vec2 texCoords;uniform sampler2D texture0;void main(){gl_FragColor = texture2D(texture0, texCoords);}";
    pixProgram = pixBuildVertexFragmentProgram(vertexCode, fragmentCode);
    glUseProgram(pixProgram);
    // Load identifiers for sampler units, for some reason before validating.
    glUniform1i(glGetUniformLocation(pixProgram, "texture0"), 0);
    //!!validateShaderProgram(self.drawingProgram)
    pixUnifLoc = glGetUniformLocation(pixProgram, "cameraMatrix");
    pixAttrLoc = glGetAttribLocation(pixProgram, "attribs");
    return 0;
}

// Make a rectangle from two triangles. Leave the OpenGL buffers bound.
int pixInitMesh() {
    float widthFrac = (float)pixWidth / pixWidth;
    float heightFrac = (float)pixHeight / pixHeight;
    GLfloat attributes[4 * 4] = {
        0., 0., 0., 0.,
        (float)pixWidth, 0., widthFrac, 0., 
        (float)pixWidth, (float)pixHeight, widthFrac, heightFrac,
        0., (float)pixHeight, 0., heightFrac};
    GLushort triangles[2 * 3] = {0, 1, 2, 0, 2, 3};
    glCreateBuffers(1, &pixAttrBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, pixAttrBuffer);
    glNamedBufferData(pixAttrBuffer, 4 * 4 * sizeof(GLfloat),
        (GLvoid *)attributes, GL_STATIC_DRAW);
    glCreateBuffers(1, &pixTriBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pixTriBuffer);
    glNamedBufferData(pixTriBuffer, 2 * 3 * sizeof(GLushort),
        (GLvoid *)triangles, GL_STATIC_DRAW);
    glVertexAttribPointer(pixAttrLoc, 4, GL_FLOAT, 0, 4 * sizeof(GLfloat), 0);
    glEnableVertexAttribArray(pixAttrLoc);
    return 0;
}



/*** Public: miscellaneous ***/

/* A pixel system program usually proceeds through these five steps:
    A. pixInitialize is invoked to set up certain resources.
    B. Other pixel system functions are invoked to configure the user interface. 
    C. pixRun is invoked. While it runs, user interface callbacks registered in 
       step B are invoked automatically.
    D. The user elects to quit the program, thus causing pixRun to terminate.
    E. pixFinalize is invoked to clean up the resources.
Sometimes a pixel system program makes multiple loops through these five steps. 
The important thing is that invocations of pixInitialize and pixFinalize come in 
non-overlapping pairs. The pixel system is designed to support only one active 
window at a time. */

/* Initializes the pixel system. This function must be called before any other 
pixel system functions. The width and height parameters describe the size of 
the window. They should be powers of 2. The name parameter is a string for the 
window's title. Returns an error code, which is 0 if no error occurred. Upon 
success, don't forget to call pixFinalize later, to clean up the pixel system. 
*/
int pixInitialize(int width, int height, const char *name) {
    pixWidth = width;
    pixHeight = height;
    pixInitGLFWGL3W(name);
    pixInitTexture();
    pixInitShaders();
    pixInitMesh();
    std::cerr << "OpenGL " << glGetString(GL_VERSION) << ", GLSL "
    << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
    return 0;
}

/* Runs the event loop. First, any pending user events are processed by their 
corresponding callbacks. Second, the time step callback is invoked. Third, if 
any drawing has occurred, then the screen is updated to reflect that drawing. 
When the user elects to quit, this function terminates. */
void pixRun() {
    while (glfwWindowShouldClose(pixWindow) == GL_FALSE) {
        glfwPollEvents();
        pixOldTime = pixNewTime;
        pixNewTime = glfwGetTime();
        if (pixUserTimeStepHandler)
            pixUserTimeStepHandler(pixOldTime, pixNewTime);
        if (pixNeedsRedisplay) {
            // In OpenGL 4.5 we might do this.
            glTextureSubImage2D(pixTexture, 0, 0, 0, pixWidth, 
               pixHeight, GL_RGB, GL_FLOAT, data(pixPixels));
            GLfloat matrix[] = {
                2.f / pixWidth, 0., 0., -1.,
                0., 2.f / pixHeight, 0., -1.,
                0., 0., -1., 0.,
                0., 0., 0., 1.};
            glUniformMatrix4fv(pixUnifLoc, 1, GL_TRUE, matrix);
            glDrawElements(GL_TRIANGLES, 3 * 2, GL_UNSIGNED_SHORT, 0);
            glfwSwapBuffers(pixWindow);
            pixNeedsRedisplay = false;
        }
    }
}

/* Deallocates the resources supporting the window. After this function is 
called, pixInitialize must be called again, before any further use of the pixel 
system. */
void pixFinalize() {    
    glDisableVertexAttribArray(pixAttrLoc);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glDeleteBuffers(1, &pixTriBuffer);
    glDeleteBuffers(1, &pixAttrBuffer);
    //glUseProgram(0);
    glDeleteProgram(pixProgram);
    glDeleteTextures(1, &pixTexture);
    delete[] pixPixels.data();
    glfwDestroyWindow(pixWindow);
    glfwTerminate();
}

/* Sets the pixel at coordinates (x, y) to the given RGB color. Coordinates are 
relative to the lower left corner of the window. */
void pixSetRGB(int index, const double (&rgb)[4]) {
    pixPixels[index] = rgb[0];
    pixPixels[index + 1] = rgb[1];
    pixPixels[index + 2] = rgb[2];
}

/* Sets all pixels to the given RGB color. */
void pixClearRGB(double red, double green, double blue) {
    for (int index = 0; index < pixPixels.size(); index += 3) {
        pixPixels[index] = red;
        pixPixels[index + 1] = green;
        pixPixels[index + 2] = blue;
    }
    pixNeedsRedisplay = true;
}

/* Inverse of pixCopyRGB. This function pastes the contents of the data array 
into the window. */
void pixPasteRGB(double (&&data)[]) {
    for (int i = 0; float &f : pixPixels)
        f = data[i++];
}



/*** Public: callbacks ***/

/* Sets a callback function for keys' being pressed. Invoked using something 
like
    pixSetKeyDownHandler(myKeyDownHandler);
where myKeyDownHandler is defined something like 
    void myKeyDownHandler(int key, int shiftIsDown, int controlIsDown,
        int altOptionIsDown, int superCommandIsDown);
The key parameter is GLFW_KEY_A, GLFW_KEY_B, ..., or GLFW_KEY_UNKNOWN, which 
are defined in GLFW/glfw3.h. The other parameters are flags describing the 
state of the modifier keys when the key was pressed. */
void pixSetKeyDownHandler(void (*handler)(int)) {
    pixUserKeyDownHandler = handler;
}

/* Sets a callback function for keys' being released. For details, see 
pixSetKeyDownHandler. */
void pixSetKeyUpHandler(void (*handler)(int)) {
    pixUserKeyUpHandler = handler;
}

/* Sets a callback function for keys' being held down. For details, see 
pixSetKeyDownHandler. */
void pixSetKeyRepeatHandler(void (*handler)(int)) {
    pixUserKeyRepeatHandler = handler;
}

/* Sets a callback function that fires once per animation frame, after all of 
the user interface callbacks. Invoked using something like 
    pixSetTimeStepHandler(myTimeStepHandler);
where myTimeStepHandler is defined something like
    void myTimeStepHandler(double oldTime, double newTime);
oldTime was the time for the previous frame; newTime is the time for the 
current frame. Both times are in seconds since the epoch (something like 1970). 
*/
void pixSetTimeStepHandler(void (*handler)(double, double)) {
    pixUserTimeStepHandler = handler;
}