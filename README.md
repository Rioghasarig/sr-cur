Paper: https://math.berkeley.edu/~oekenta/Spectrum_Revealing_CUR__MLDM_.pdf

## SR-CUR 
Spectrum-Revealing CUR (SR-CUR) is a tool for the low rank CUR decomposition of sparse matrices. SR-CUR  implements a novel CUR decomposition method that allows it to scale to larger matrices and produce more accurate factorizations than previous methods. It employs a variant of LU pivoting with complete pivoting to select the appropriate rows and columns. In particular, it makes use of a random sketch of the Schur complement to speed up the search for desirable rows and columns. It is highly efficient, outperforming  Matlab's svds by a factor of 10 on large-scale matrices.  For more details see the following [paper][CHEN20] and [presentation][BEBOP] .
  
  [CHEN20]: https://proceedings.mlr.press/v108/chen20a.html
  [BEBOP]: https://math.berkeley.edu/~oekenta/Sparse_CUR_Talk.pdf

We make use of the [LUSOL][LUSOL] package to compute the sparse LU factorizations.
## LUSOL

[LUSOL][LUSOL] maintains LU factors of a square or rectangular sparse matrix.

This repository provides [LUSOL][LUSOL] source code and a Matlab interface.

The code is distributed under the terms of the MIT License or the BSD License.

  [LUSOL]: http://web.stanford.edu/group/SOL/software/lusol/

## Contents

* `gen/`: code generation scripts and specification files
* `matlab/`: Matlab interface code
* `src/`: LUSOL Fortran code
* `LICENSE`: [Common Public License][CPL]
* `makefile`: GNU Make file to build interface

## Download and install

Installation simply requires adding the `matlab` subdirectory to your Matlab
path.  This may be done with Matlab's [`addpath`][ADDPATH] function.

Please follow the directions below to build the interface

  [ADDPATH]: http://www.mathworks.com/help/matlab/ref/addpath.html

## Usage Examples



Stable-CUR

```
% A - large sparse matrix
% l - oversampling parameter
% k - rank of the approximation. k <= l  

mycur = srcur_obj(A, 'pivot', 'TRP', 'rank', l); 
f = 2; maxswaps = l;
mycur.srlu(f,maxwaps); 
[C,U,R] = mycur.stable_cur(k); 
%Compute error
err = mycur.cur_error(); % Equivalent to norm(mylu.Apq-C*U*R,'fro')/norm(mylu.Apq,'fro') 

```

## Build

### Environment

The build environments as of 2016-01-26 are:

- Fedora 21 & Matlab 2013b
- Mac OS X 10.11 & Matlab 2015b

Building the LUSOL interface in other environments may require modification of
`makefile` and `matlab/lusol_build.m`.

### Requirements

Linux:

* `make`
* `gcc`
* `gfortran`
* Matlab

Mac:

* [Xcode][XC] for `clang` and `make`
* `gfortran` (possibly via [Homebrew][HB])
* Matlab

  [HB]: http://brew.sh/
  [XC]: http://itunes.apple.com/us/app/xcode/id497799835

Notes:

* The `matlab` binary must be on the system `PATH`.
* `python3` is required to generate the interface code.  However, the interface
  code is pre-generated and included in the repository.
* It is necessary to launch Xcode and accept the license agreement before
  building the interface.
* The smaller Xcode Command Line Tools package does not appear to work with
  Matlab 2015b.  The full Xcode install is required.

### Setup `mex`

Matlab's `mex` compiler driver must be configured to use the appropriate `C`
compiler.  This can be achieved by executing `mex -setup` from the Matlab prompt
or operating system terminal.  On Linux the selected compiler must be the
correct version of `gcc`.  On Mac OS X 10.9 the selected compiler must be
`clang`.  It is a good idea to match the compiler suggested on the Mathworks
[supported compilers][MC] page.  See this [note][MEX-XCODE-7] on Matlab
compatibility with Xcode 7.

  [MC]: http://www.mathworks.com/support/compilers/
  [MEX-XCODE-7]: http://www.mathworks.com/matlabcentral/answers/246507-why-can-t-mex-find-a-supported-compiler-in-matlab-r2015b-after-i-upgraded-to-xcode-7-0#answer_194526

### Install `gfortran` on Mac OS X

1. Install [Homebrew][HB]
3. `$ brew install gcc`

### Steps

From the base directory:

```
$ make
$ make matlab
```



### Notes

The basic requirements to build LUSOL are GNU `make`, `gfortran`, a C compiler,
and Matlab.  The build works with `gcc` on Linux and `clang` on Mac OS X.  It
may be necessary to modify the compiler variables in the `makefile` (`CC`,
`F90C`, and `F77C`) depending on the operating system and environment.

The `matlab` executable must be on the system path.  On Mac OS X with Matlab
R2015b this is achieved with:

```
$ export PATH=/Applications/MATLAB_R2015b.app/bin:$PATH
```

The `makefile` may have to be modified on Mac OS X depending on versions of
Matlab and `gfortran`.


