# math submodule

This repository is a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules)
providing C++ math utilities for other projects, including:

* <tt>QuadraticPolynomial</tt> : A quadratic polynomial.
* <tt>CubicPolynomial</tt> : A cubic polynomial.
* <tt>AnalyzedCubic</tt> : Cached analysis of a CubicPolynomial: one of the local extremes and the inflection point.
* <tt>Direction</tt> : A unit vector.
* <tt>Line</tt> : A line.
* <tt>LinePiece</tt> : Two points representing the beginning and end of a line piece.
* <tt>Point</tt> : A point in a 2D plane.
* <tt>Polynomial</tt> : A polynomial with double coefficients (not as specialized as QuadraticPolynomial and CubicPolynomial).
* <tt>Vector</tt> : A vector in a 2D plane.
* <tt>bracket\_zero</tt> : Utility to find the zero of a given function.

The root project should be using
[cmake](https://cmake.org/overview/)
[cwm4](https://github.com/CarloWood/cwm4) and
[libcwd](https://github.com/CarloWood/libcwd).

## Checking out a project that uses the math submodule.

To clone a project example-project that uses math simply run:

    git clone --recursive <URL-to-project>/example-project.git
    cd example-project
    AUTOGEN_CMAKE_ONLY=1 ./autogen.sh

The ``--recursive`` is optional because ``./autogen.sh`` will fix
it when you forgot it.

When using [GNU autotools](https://en.wikipedia.org/wiki/GNU_Autotools) you should of course
not set ``AUTOGEN_CMAKE_ONLY``. Also, you probably want to use ``--enable-mainainer-mode``
as option to the generated ``configure`` script.

In order to use ``cmake`` configure as usual, for example to build with 16 cores a debug build:

    mkdir build_debug
    cmake -S . -B build_debug -DCMAKE_MESSAGE_LOG_LEVEL=DEBUG -DCMAKE_BUILD_TYPE=Debug -DCMAKE_VERBOSE_MAKEFILE=ON -DEnableDebugGlobal:BOOL=OFF
    cmake --build build_debug --config Debug --parallel 16

Or to make a release build:

    mkdir build_release
    cmake -S . -B build_release -DCMAKE_BUILD_TYPE=Release
    cmake --build build_release --config Release --parallel 16

## Adding the math submodule to a project

To add this submodule to a project, that project should already
be set up to use [utils](https://github.com/CarloWood/ai-utils).

Then simply execute the following in a directory of that project
where you want to have the ``math`` subdirectory (the
root of the project is recommended as that is the only thing
I've tested so far):

    git submodule add https://github.com/CarloWood/math.git

This should clone math into the subdirectory ``math``, or
if you already cloned it there, it should add it.

Checkout [math-testsuite](https://github.com/CarloWood/math-testsuite)
for an example of a project that uses this submodule.
