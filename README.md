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

## How To Use

SLIM surface rendering can be done by executing:

```
% slimviewGPU.exe in.slim2t
```

SLIM (.slim2t) data can be created by [mesh2slim](https://github.com/kanait/mesh2slim). Please check it out.

On a window, rendering modes can be changed by pressing keys:

- Press "1": Render by billboards, LOD disable (rendered by leaf nodes)
- Press "2": Render by point sprites, LOD disable
- Press "3": Render by billboards, LOD by depth of a slim tree
- Press "4": Render by point sprites, LOD by depth of a slim tree
- Press "5": Render by billboards, LOD by error
- Press "6": Render slim balls, LOD disable
- Press "7": Renderp slim balls, LOD by error
- Press "m": Down LOD levels (finer levels)
- Press "n": Up LOD levels (coaser levels)

## Prerequisites

The following libraries are required for successfully compiling this software.

### [Cg Toolkit](https://developer.nvidia.com/cg-toolkit/)

Our code uses Cg Toolkit library from nVIDIA. The latest version (Cg 3.1 release) is supported.

### [zlib](https://zlib.net/)

### [libpng](https://www.libpng.org/pub/png/libpng.html)

### [glew](https://glew.sourceforge.net/)

### [glfw](https://www.glfw.org/)

### [vecmath-cpp](https://github.com/yuki12/vecmath-cpp)
### [render](https://github.com/kanait/render)

When you execute "git clone" with "--recursive" option, you will also get vecmath-cpp and render libraries as a submodule "external/vecmath-cpp" and "external/render":

```
git clone https://github.com/kanait/slimviewGPU.git --recursive
```

## Authors

* **[Takashi Kanai](https://graphics.c.u-tokyo.ac.jp/hp/en/)** - The University of Tokyo

## License

This software is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
