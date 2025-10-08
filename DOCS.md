# Documentation Generation

`magpie-matlab` is following the [MATLAB's guidelines for documentation](https://uk.mathworks.com/help/matlab/matlab_prog/display-custom-documentation.html) for documentation structure and generation. Documentation files are written using the [publishing markup](https://uk.mathworks.com/help/matlab/matlab_prog/marking-up-matlab-comments-for-publishing.html) system and the [`publish`](https://uk.mathworks.com/help/matlab/ref/publish.html) command. 

Files which are documented are contained in the `docs/` directory and are ending in `_help.m`. The html files generated are th ones that will be used with the [`docs`](https://uk.mathworks.com/help/matlab/ref/doc.html) commmand when the toolbox is available.

## Generating Documentation

Documentation is generated with the `docs/build_docs.m` script

In MATLAB, change directory and run `build_docs.m`. This will generate a contents `index.html` page and an entry for each function ending in `_help.m` in the `docs/html/` directory.

```matlab
cd docs
build_docs
```

## Documentation Style

Documentation comments are written _above_ the function decleration. Using `%%` as a heading, the documents should adhere to the following struture

- function name
  - normal text summary
- Syntax
  - list Function arguments as a code block
- Description
  - List of each function overload and what it returns
- Example
  - A short working example
- Input Arguments
  - List of arguments and there meaning
- See Also
  - Pipe seperated list of relevant functions

e.g.
```
%% fidimat 
%
% Generate Finite Difference Spatial Sparcity Matrices
%
%% Syntax
%
%   fdmat = fidimat(l,ord)
%   fdmat = fidimat(l,m,ord)
%   fdmat = fidimat(l,m,ord,bctype)
%
%% Description
%
% |fdmat = fidimat(l,ord)| generates 1D stencil |fdmat| with order given by |ord|
%
% |fdmat = fidimat(l,m,ord)| generate a stencil |fdmat| with order given by |ord|
% for a 2D system of size |l|-by-|m|. Boundary condition defaults to simply supported
%
% |fdmat = fidimat(l,m,ord,bctype)| generate stencil |fdmat| with order given
% by |ord| for a 2D system of size |l|-by-|m| with specified boundary conditions bctype.
%
%% Example
%
%   Nx = 100;
%   Ny = 100;
%   XXYY = fidimat(Ny,Nx,'xxyy', 1);  % 2D second order matrix
%
%% Input Arguments
% * m       % number of total grid points X axis
% * l       % number of total grid points Y axis
% * ord     % order of the matrix (string)
% * bctype  % boundary condition type: 1: simply supported, 2: clamped
%
%   % Valid orders
%     ['x-','x+','x.','xx','xxxx',
%     'y-','y+','y.','yy','yyyy',
%     'grad','xy','xxyy','laplace','biharm','I'];
% 
%% See Also
% <./bhmat_help.html  |bhmat|>
```