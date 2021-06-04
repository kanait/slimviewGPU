# slimviewGPU

slimviewGPU is a research software for displaying SLIM surfaces on programmable GPUs, which includes an implementation of the following paper:

Takashi Kanai, Yutaka Ohtake, Hiroaki Kawata, Kiwamu Kase: ``GPU-based Rendering of Sparse Low-degree IMplicit Surfaces'', 4th International Conference on Computer Graphics and Interactive Techniques in Australasia and the Southeast Asia (GRAPHITE 2006), pp.165-171, 2006.

*[SLIM](https://dl.acm.org/doi/10.5555/1281920.1281944) (Sparse Low-degree IMplicits) is an implicit surface representation delivering an accurate approximation to a set of points scattered over a smooth surface.*

This software was originally developed in 2005-2006 and was renovated in 2021 so as to build successfully by Visual Studio 2019.

## Getting Started

This software can run only on Windows. 
At first, you may try to use binary release (x64), 
which is available from [here](https://github.com/kanait/slimviewGPU/releases/tag/v1.0).
Uncompress zip file and then execute run_pai.x64.bat.

## Prerequisites

The following libraries are required for successfully compiling this software.

### [Cg Toolkit](https://developer.nvidia.com/cg-toolkit/)

Our code uses Cg Toolkit library from nVIDIA. The latest version (Cg 3.1 release) is supported.

### [zlib](https://zlib.net/)

### [libpng](http://www.libpng.org/pub/png/libpng.html)

### [glew](http://glew.sourceforge.net/)

### glut (included in Cg Toolkit library)

### [vecmath-cpp](https://github.com/yuki12/vecmath-cpp)
### [render](https://github.com/kanait/render)

When you execute "git clone" with "--recursive" option, you will also get vecmath-cpp and render libraries as a submodule "external/vecmath-cpp" and "external/render":

```
git clone https://github.com/kanait/slimviewGPU.git --recursive
```

## Authors

* **[Takashi Kanai](https://graphics.c.u-tokyo.ac.jp/hp/en/)** - The University of Tokyo

## License

This software is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
