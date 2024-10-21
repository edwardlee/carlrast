## Run me
#### Download [gl3w](https://github.com/skaslev/gl3w)

`python gl3w_gen.py`

`gcc -c gl3w.c`

Move gl3w.o to git root

#### Download [glfw](https://www.glfw.org/download)

Move include/GLFW folder to /usr/local/include

Move contents corresponding to your architecture to /usr/local/lib

Alternatively, installing with brew (glfw3) or apt (libglfw3 or libglfw3-dev) will probably work

#### Download [stb_image.h](https://github.com/nothings/stb/blob/master/stb_image.h)

Move it to git root

### Apple Silicon
My MacOS 15 SDK is missing `copy_to` in `std::experimental::simd<float, std::experimental::simd_abi::__vec_ext<4>>`. That seems very silly but it wouldn't be the first time Apple outsillied me. Compile against 14 instead.

`g++ -c 040pixel.cc -std=c++20 -Ofast -march=native -isysroot/Library/Developer/CommandLineTools/SDKs/MacOSX14.sdk -I/usr/local/include`

`g++ 350mainClipping.cc -std=c++20 040pixel.o gl3w.o -lglfw3 -Ofast -march=native -framework Cocoa -framework IOKit`

### WSL1
Install [VcXsrv](https://github.com/marchaesen/vcxsrv/releases)

Launch it

Maybe add `export DISPLAY=:0` in your .bashrc

`g++ -c 040pixel.cc -std=c++20 -Ofast -march=native`

`g++ 350mainClipping.cc -std=c++20 040pixel.o gl3w.o -lglfw -Ofast -march=native`
