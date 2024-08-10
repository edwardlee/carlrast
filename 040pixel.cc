/*
On Ubuntu, compile with...
    g++ -c 040pixel.cc -std=c++20 -Ofast -march=native
...and then link with a main program by for example...
    g++ 350mainClipping.cc -std=c++20 040pixel.o gl3w.o -lglfw -Ofast -march=native
*/


/*** Private: infrastructure ***/

#include <iostream>
#include <span>
#include <GL/gl3w.h>
#include <GLFW/glfw3.h>

// Global variables.
GLFWwindow *pixWindow;
int pixWidth, pixHeight;
GLuint pixTexture;
std::span<float> pixPixels;
bool pixNeedsRedisplay = false;
GLuint pixProgram;
double pixNewTime;
void (*pixUserKeyDownHandler)(int);
void (*pixUserKeyUpHandler)(int);
void (*pixUserKeyRepeatHandler)(int);
void (*pixUserTimeStepHandler)(double, double);

GLuint pixBuildVertexFragmentProgram(const GLchar *vertexCode, 
    const GLchar *fragmentCode) {
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexCode, nullptr);
    glCompileShader(vertexShader);
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentCode, nullptr);
    glCompileShader(fragmentShader);
    GLuint program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    glLinkProgram(program);
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    return program;
}



/*** Private: GLFW handlers ***/

void pixHandleError(int error, const char *description) {
    std::cerr << "GLFW code " << error << "...\n"
        << description << std::endl;
}

// key is GLFW_KEY_A, GLFW_KEY_B, etc. Or GLFW_KEY_UNKNOWN.
// scancode is a platform-dependent scan code for the key.
// action is GLFW_PRESS, GLFW_RELEASE, or GLFW_REPEAT.
// mods has bitwise masks for...
//     GLFW_MOD_SHIFT, GLFW_MOD_CONTROL, GLFW_MOD_ALT, or GLFW_MOD_SUPER
// ...which on macOS mean shift, control, option, command.
void pixHandleKey(GLFWwindow *window, int key, int, int action, int) {
    if (action == GLFW_PRESS && pixUserKeyDownHandler)
        pixUserKeyDownHandler(key);
    if (action == GLFW_RELEASE && pixUserKeyUpHandler)
        pixUserKeyUpHandler(key);
    else if (action != GLFW_RELEASE && pixUserKeyRepeatHandler)
        pixUserKeyRepeatHandler(key);
}

/*** Private: initializer ***/

// Initialize the user interface and OpenGL.
void pixInitGLFWGL3W(std::string_view name) {
    glfwSetErrorCallback(pixHandleError);
    glfwInit();
    glfwWindowHint(GLFW_RESIZABLE, 0);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    pixWindow = glfwCreateWindow(pixWidth, pixHeight, data(name), nullptr, nullptr);
    glfwMakeContextCurrent(pixWindow);
    glfwSetKeyCallback(pixWindow, pixHandleKey);
    int width, height;
    glfwGetFramebufferSize(pixWindow, &width, &height);
    gl3wInit();
    glViewport(0, 0, width, height);
}

// Create the texture.
void pixInitTexture() {
    pixPixels = std::span(new float[3 * pixWidth * pixHeight], 3 * pixWidth * pixHeight);
    glGenTextures(1, &pixTexture);
    glBindTexture(GL_TEXTURE_2D, pixTexture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, pixWidth, pixHeight, 0, 
        GL_RGB, GL_FLOAT, data(pixPixels));
}

GLchar vertexCode[] = R"(#version 410
out vec2 texCoords;
void main() {
    vec2 vertices[3] = vec2[3](vec2(-1,-1), vec2(3,-1), vec2(-1, 3));
    gl_Position = vec4(vertices[gl_VertexID], 0, 1);
    texCoords = 0.5 * gl_Position.xy + vec2(0.5);
})";
GLchar fragmentCode[] = R"(#version 410
out vec4 color;
in vec2 texCoords;
uniform sampler2D texture0;
void main() {
    color = texture(texture0, texCoords);
})";
// Build the shader program and leave it bound.
void pixInitShaders() {    
    pixProgram = pixBuildVertexFragmentProgram(vertexCode, fragmentCode);
    glUseProgram(pixProgram);
    GLuint emptyVAO;
    glGenVertexArrays(1, &pixProgram);
    glBindVertexArray(pixProgram);
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
void pixInitialize(int width, int height, std::string_view name) {
    pixWidth = width;
    pixHeight = height;
    pixInitGLFWGL3W(name);
    pixInitTexture();
    pixInitShaders();
    std::cerr << "OpenGL " << glGetString(GL_VERSION) << ", GLSL "
    << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
}

/* Runs the event loop. First, any pending user events are processed by their 
corresponding callbacks. Second, the time step callback is invoked. Third, if 
any drawing has occurred, then the screen is updated to reflect that drawing. 
When the user elects to quit, this function terminates. */
void pixRun() {
    while (glfwWindowShouldClose(pixWindow) == GL_FALSE) {
        glfwPollEvents();
        double pixOldTime = pixNewTime;
        pixNewTime = glfwGetTime();
        if (pixUserTimeStepHandler)
            pixUserTimeStepHandler(pixOldTime, pixNewTime);
        // if (pixNeedsRedisplay) {
            // pixNeedsRedisplay = false;
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, pixWidth, 
                            pixHeight, GL_RGB, GL_FLOAT, data(pixPixels));
            glDrawArrays(GL_TRIANGLES, 0, 3);
            glfwSwapBuffers(pixWindow);
        // }
    }
}

/* Deallocates the resources supporting the window. After this function is 
called, pixInitialize must be called again, before any further use of the pixel 
system. */
void pixFinalize() {    
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

#include <experimental/simd>
using namespace std::experimental;
/* Sets all pixels to the given RGB color. */
void pixClearRGB(const float (&f)[3]) {
    native_simd<float> a([&f](int i){return f[i % 3];});
    native_simd<float> b([&f](int i){return f[(i+2) % 3];});
    native_simd<float> c([&f](int i){return f[(i+1) % 3];});
    constexpr int k = 3 * a.size();
    for (int index = 0; index < pixPixels.size(); index += k) {
        a.copy_to(&pixPixels[index], element_aligned);
        b.copy_to(&pixPixels[index]+8, element_aligned);
        c.copy_to(&pixPixels[index]+16, element_aligned);
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
