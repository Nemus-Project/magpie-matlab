# MAGPIE
> *An open-source plate vibration simulation and analysis framework*

M. Ducceschi`*`, M. Hamilton`*`, A. Mousseau`^`, S. Duran`*`

`*` NEMUS Lab - University of Bologna

`^` LAUM - Laboratoire d'Acoustique de l'Université du Mans


## About

MAGPIE is an open-source framework for simulating the dynamics of elastic plates. The plate equation is derived from the *orthotropic* Kirchhoff-Love plate theory. It can then effectively simulate the dynamics of moderately thin wooden and composite plates. The geometry is restricted to rectangular, and the thickness is assumed constant. These limitations may be overcome in future releases. MAGPIE can simulate boundary conditions of classic type (free, simply-supported and clamped) and generally elastic boundary support, controlling the flexural displacement and boundary rotation through applied elastic forces and moments. 

MAGPIE is a framework allowing both *direct* and *inverse modelling* of orthotropic sheet materials, such as wood and composites, displaying orthotropic behaviour. Note that the isotropic case can be recovered as a special case. 

** Direct modelling consists of supplying MAGPIE with the appropriate elastic parameters (two Young's moduli, one shear modulus and one Poisson ratio), material parameters (density) and geometric parameters (edge lengths and thickness) plus the values of the elastic boundary constants. The output of the direct modelling can be:

1. Eigenvalues and eigenvectors (i.e., modal shapes)
2. A time-domain simulation of a struck plate using a simple raised cosine (simulating an impact) or sinusoidal input 

** Inverse modelling consists of supplying MAGPIE with a number of measured modal frequencies from an experimental plate, along with the plate's known material and geometric constants (density, edge lengths, thickness) to compute the elastic constants.

## Table of Contents
**[Features](#features)**<br>
**[The Orthotropic Kirchhoff Model](#the-orthotropic-kirchhoff-model)**<br>
**[Discretisation](#discretisation)**<br>
**[Braces and static loads ](#braces-and-static-loads)**<br>
**[Structure](#structure)**<br>
**[Direct Modelling](#direct-modelling)**<br>


## Features

Here is a list of features currently supported by MAGPIE:

1. orthotropic wave propagation (includes isotropic as a special case)
2. rectangular geometry
3. constant thickness
4. generally elastic boundary conditions
5. possibility of adding ribs (Euler-Bernoulli beams)
6. possibility of adding lumped static loads and stiffeners
7. frequency domain analysis
8. time domain simulation using two input forces: impulse and sines
9. arbitrary amount of output points
10. stress computation, plus output displacement, velocity and acceleration
11. Rayleigh damping
12. sound synthesis of displacement, velocity and acceleration at output points
13. Inverse modelling to compute the elastic constants from a set of measured experimental frequencies

## The Orthotropic Kirchhoff Model

The vibration of an orthotropic plate can be described via the Kirchhoff-Love model. This dynamical model describes the time evolution of the flexural displacement $u = u(x,y,t) : \mathcal{V} \times \mathbb{R}^+_0$. Here,  ${\bf x} \in \mathcal{V} := [0,L_x] \times [0,L_y]$ is the domain of definition, a rectangle with side lengths $L_x, L_y$. Here, and in what follows, $x$ denotes the longitudinal orthotropic direction, $y$ is radial, and $z$ is tangential. In quarter-sawing, $z$ is the direction along the thickness of the board. With this notation, the system reads:

$$\rho \zeta \partial_t^2 u(x,y,t) = -D_1\partial_x^4 u(x,y,t) - (D_2+D_4)\partial_y^2\partial_x^2 u(x,y,t) - D_3\partial_y^4 u(x,y,t) := -\mathcal{B}u(x,y,t)$$

Here, $B$ denotes the orthotropic spatial differential operator generalising the biharmonic in the isotropic case. The model depends on four rigidity constants, denoted by $D_\circ$, defined in terms of Young's moduli $E_\circ$, the shear modulus $G_{xy}$ and the Poisson's ratios $\nu_\circ$ as follows:

$$ D_1 := \frac{E_x\zeta^3}{12(1-\nu_{x}\nu_{y})} \quad D_3 := \frac{E_y\zeta^3}{12(1-\nu_{x}\nu_{y})} \quad D_4 := \frac{G_{xy}\zeta^3}{3} \quad D_2 := \frac{\nu_{yx}E_{x}\zeta^3}{6(1-\nu_{y}\nu_{y})} $$

Note that, because of the symmetry of the compliance matrix, one out of five elastic constants is fixed. Usually, one selects

$$\nu_{y} = \nu_{x}E_{y}E_{x}^{-1}$$

thus leaving the following elastic constants to be specified at the input: $E_x,E_y,G_{xy},\nu_{x}$. Note that the model also depends on the thickness $\zeta$ (assumed constant throughout $\mathcal V$), the side lengths $L_x,L_y$ and the density $\rho$ (also assumed constant).


```

                          [K_{xL_y},R_{xL_y}]
                  ---------------------                       ^
                 |                     |                      |
                 |                     |                      | 
 [K_{0y},R_{0y}] |                     |                      |
                 |                     |                       L_y
                 |                     |[K_{L_xy},R_{L_xy}]   |
                 |                     |                      |
                 |                     |                      |
                  ---------------------                       ˇ
               [K_{x0},R_{x0}]

                 <-------- L_x -------->

```

The equation must be supplied with appropriate boundary conditions. Here, conditions of elastic type are assumed. Referring to the diagram above, the boundary conditions are:

- at $x = 0$
  * Balance of Forces:      $$K_{0y} u = -D_1 \left( \partial_x^3 u + \left( \frac{D_4}{D_1} + \nu_{y} \right) \partial_x\partial_y^2 u\right)$$
  * Balance of Moments:     $$R_{0y} \partial_x u = D_1 \left(\partial_x^2 u + \nu_y \partial_{yy} u \right)$$


- at $x = L_x$
  * Balance of Forces:      $$K_{L_x y} u = D_1 \left( \partial_x^3 u + \left( \frac{D_4}{D_1} + \nu_{y} \right) \partial_x\partial_y^2 u\right)$$
  * Balance of Moments:     $$R_{L_x y} \partial_x u = -D_1 \left(\partial_x^2 u + \nu_y \partial_{yy} u \right)$$


- at  $y = 0$
  * Balance of Forces:      $$K_{x0} u = -D_3 \left( \partial_y^3 u + \left( \frac{D_4}{D_3} + \nu_{x} \right) \partial_y\partial_x^2 u\right)$$
  * Balance of Moments:     $$R_{x0} \partial_y u = D_3 \left(\partial_y^2 u + \nu_x \partial_{xx} u \right)$$


- at  $y = L_y$
  * Balance of Forces:      $$K_{x L_y} u = D_3 \left( \partial_y^3 u + \left( \frac{D_4}{D_3} + \nu_{x} \right) \partial_y\partial_x^2 u\right)$$
  * Balance of Moments:     $$R_{x L_y} \partial_y u = -D_3 \left(\partial_y^2 u + \nu_x \partial_{xx} u \right)$$

Note a clamped edge is recovered by setting $(K_\circ,R_\circ)$ to very large values, yielding $u \approx 0, \partial_n u \approx 0$, where $n$ is the direction normal to the boundary. A simply-supported edge is recovered by setting $K_\circ$ to a large value, and $R_\circ$ to zero (vanishing applied moment). A free edge is obtained by setting $K_\circ = R_\circ = 0$ (vanishing applied force and moment). 

A further condition arises at the corners, namely

$$\partial_x\partial_y u = 0$$  at a corner. 

## Discretisation

Discretisation of the equation of motion is performed on a two-dimensional grid of points. Let $h_x$ be the grid spacing along $x$, and let $h_y$ be the grid spacing along $y$. One may approximate the continuous function $u(x,y)$ using the grid function $u_{l,m} \approx u(lh_x,mh_y)$. Then, the following are used:

$$\delta_{x\cdot} u(x,y) := (2h_x)^{-1}(u(x+h_x,y)-u(x-h_x,y)) \approx \partial_x u(x,y)$$

$$\delta_{xx} u(x,y) := (h_x)^{-2}(u(x+h_x,y)-2u(x,y) + u(x-h_x,y)) \approx \partial^2_{x} u(x,y)$$

$$\delta_{xxxx} u(x,y) := \delta_{xx}\delta_{xx} u(x,y) \approx \partial^4_{x} u(x,y)$$

When applied to the grid function $u_{l,m}$, the difference operators become:

$$\delta_{x\cdot} u_{l,m} := (2h_x)^{-1}(u_{l+1,m}-u_{l-1,m}) $$

$$\delta_{xx} u_{l,m} := (h_x)^{-2}(u_{l+1,m}-2u_{l,m} + u_{l-1,m}) $$

$$\delta_{xxxx} u_{l,m} := (h_x)^{-4}(u_{l+2,m}-4u_{l+1,m} + 6u_{l,m} -4u_{l-1,m} + u_{l-2,m}) $$

Similar definitions are used for the differences along the $y$ direction. A discrete version of the equation of motion in terms of the difference operators is obtained as:


$$\rho \zeta \ddot u_{l,m}(t) = -D_1\delta_{xxxx} u_{l,m}(t) - (D_2+D_4)\delta_{xx}\delta_{yy} u_{l,m} - D_3\delta_{yyyy} u_{l,m}(t)$$

Note that overdots now indicate total time derivatives. Note as well that the spatial indices are integer numbers restricted to the intervals $0 \leq l \leq N_x = \frac{L_x}{h_x}$ and $0 \leq m \leq N_y = \frac{L_y}{h_y}$, yielding $N_x+1$ points along $x$, and $N_y+1$ points along $y$, including grid point on the boundary. 

Boundary conditions follow as a direct discretisation of the continuous conditions above. Thus

- at $l = 0$
  * Balance of Forces:      $$K_{0y} w_{0,m} = -D_1 \left( \delta_{x\cdot}\delta_{xx} u_{0,m} + \left( \frac{D_4}{D_1} + \nu_{y} \right) \delta_{x\cdot}\delta_{yy} u_{0,m}\right)$$
  * Balance of Moments:     $$R_{0y} \delta_{x\cdot} u_{0,m} = D_1 \left(\delta_{xx} u_{0,m} + \nu_y \delta_{yy} u_{0,m} \right)$$
 
- at $l = N_x$
  * Balance of Forces:      $$K_{L_x y} w_{N_x,m} = -D_1 \left( \delta_{x\cdot}\delta_{xx} u_{N_x,m} + \left( \frac{D_4}{D_1} + \nu_{y} \right) \delta_{x\cdot}\delta_{yy} u_{N_x,m}\right)$$
  * Balance of Moments:     $$R_{L_x y} \delta_{x\cdot} u_{N_x,m} = D_1 \left(\delta_{xx} u_{N_x,m} + \nu_y \delta_{yy} u_{N_x,m} \right)$$
 
- at  $m = 0$
  * Balance of Forces:      $$K_{x0} u_{l,0} = -D_3 \left( \delta_{y\cdot}\delta_{yy} u_{l,0} + \left( \frac{D_4}{D_3} + \nu_{x} \right) \delta_{y\cdot}\delta_{xx} u_{l,0}\right)$$
  * Balance of Moments:     $$R_{x0} \delta_{y\cdot} u_{l,0} = D_3 \left(\delta_{yy} u_{l,0} + \nu_x \delta_{xx} u_{l,0} \right)$$


- at  $y = L_y$
  * Balance of Forces:      $$K_{x L_y} u_{l,N_y} = D_3 \left( \delta_{y\cdot}\delta_{yy} u_{l,N_y} + \left( \frac{D_4}{D_3} + \nu_{x} \right) \delta_{y\cdot}\delta_{xx} u_{l,N_y}\right)$$
  * Balance of Moments:     $$R_{x L_y} \delta_{y\cdot} u_{l,N_y} = -D_3 \left(\delta_{yy} u_{l,N_y} + \nu_x \delta_{xx} u_{l,N_y} \right)$$

with analogous discretisations holding for the boundary conditions along the other edges. An important condition arises at the corners, as:

$$\delta_{x\cdot}\delta_{y\cdot}u_{\circ,\circ} = 0$$ at a corner.

### Stencil 

The semidiscrete system above may be recast in a convenient matrix-vector form. To that end, a state vector  $\mathbf{u}$ is constructed by stacking consecutive strips of points as per the image below.

<p align="center">
<img src="/img/stackedGrid1.png" width="300" />
</p>

Each coloured strip is composed of $N_y+1$ points. The vector is then of size $(N_x+1)(N_y+1)\times 1$. The semidiscrete equation of motion given above is then written compactly as:

$$\rho \zeta \ddot {\mathbf u}(t)  =  {\mathbf B}{\mathbf u}(t)$$

where the matrix $\mathbf{B}$ is of size $(N_x+1)(N_y+1)\times (N_x+1)(N_y+1)$. The matrix is clearly sparse, since this represents the "locality" of the difference operators. A sparsity pattern is given in the figure below. 

<p align="center">
<img src="/img/sparsitypattern.png" width="200" />
</p>

At a generic inner point in the domain $(l,m)$ such that $(l\pm 2,m\pm 2)$ are all within the domain or on the boundary, the difference operator takes on the particular stencil below, a 13-point stencil generalising the biharmonic difference operator in the isotropic case:

<p align="center">
<img src="/img/OrthoGrid1.png" width="500" />
</p>

In the above: $$a_x := \frac{D_1}{h_x^4}$$, $$a_y := \frac{D_3}{h_y^4}$$, $$a_{xy} := \frac{(D_2+D_4)}{h_x^2h_y^2}$$. Note that, in the above, the grid spacings $h_x$, $h_y$ are chosen equal, but this need not be the case in general. When the difference operator acts on points nearby the boundary, the stencil looks outside the plate grid. An example is represented by corners, such as in the image below.

<p align="center">
<img src="/img/BCs.png" width="300" />
</p>

The values of the points in shaded "ghost" region must be set according to the numerical boundary conditions above, including the important corner condition. Doing so, one is able to revert the value of the difference operator at each point on the domain using points defined on the plate grid. 

**The expressions for the resulting coefficients are cumbersome**. They appear in the appropriate Matlab functions. 

The block form of the matrix is as per the image below. Blocks of the same colour have the same structure. All middle blocks, represented by a shaded blue area, have the same structure. 
<p align="center">
<img  width="600" src="/img/BiharmStructure1.png" />
</p>

The diagram above shows the overall structure of the difference operator. Yellow blocks are pentadiagonal, blue blocks are tridiagonal, and red blocks are diagonal. 


## Braces and static loads 

MAGPIE allows the coupling of the plate structure to stiffening ribs and pointwise loads. The equations are modified as follows:

$$ \rho \zeta \partial_t^2 u(x,y,t) + \sum_j M_j \delta_j(x,y) (\partial_t^2 u(x_j,y_j,t)) = -B u(x,y,t) - \sum_{i}\theta_i(x,y) f_i $$

$$ \rho_i A_i \partial_t^2 w_i(z_i,t) = -E_iI_i \partial_{z_i}^4 w_i(z_i,t) + f_i $$

where $j \in [1,N_l]$ and $i \in [1,N_r]$ with $N_l$ being the total number of static loads, and $N_r$ being the total number of braces. The braces are modelled as Euler-Bernoulli beams. In the above, $\delta_j := \delta(x-x_j)\delta(y-y_j)$ is a two-dimensional Dirac delta. Furthermore, $\theta_i$ is a kind of projector of $x,y$ onto $z_i$, such that $\int \theta_i g(x,y) dxdy = \int_{0}^L g(z_i) dz_i$. Besides the above equation, a further kinematic condition arises, namely that the displacement of the ribs equals the displacement of the plate along $z_i$: 

$$\int \theta_i u dxdy = w_i(z_i)$$

These elements change the vibrational modes and the response of the plate to applied dynamic (i.e. time-dependent) loads. 

## Structure

All source code is found in the `src/` directory of the repository. 

## Direct Modelling

Direct modelling can be run in both the frequency and the time domains. 

### Frequency Domain Analysis
For the frequency domain, the simulation parameters are specified in the [FreqDomainLauncher.m](https://github.com/Nemus-Project/magpie-matlab/blob/b9f2f20d2ebe36d25f70e85deda9e8486661e0a5/orthotropic/src/FreqDomainLauncher.m). 

User-editable parameters appear below this line 

```
%-------------------------------------------------------------------------
% CUSTOM PARAMETERS (user can change these)

```

and up to 

```
% END CUSTOM PARAMETERS 
%-------------------------------------------------------------------------
```

First, global simulation parameters are set 

```
%-- general parameters
fmax = 2000 ; % maximum frequency to be computed
ppw  = 20 ; % points per wavelength at maximum frequency. Choose 3 <= ppw 
Nmodes = 27 ; % select total number of modes to be computed. If Nmodes = [], all possible modes are computed
%--------------------
```

Here `fmax` is the largest frequency to be computed. `ppw` is the number of points per wavelength at such frequency. The grid spacings in the x and y directions are set using these two parameters as 

$$ h_\circ = v^\phi_\circ(\text{fmax})/\text{fmax}/\text{ppw} $$

The number of grid points along x and y is set via such grid spacings, as 

$$N_\circ = L_\circ / h_\circ$$

with $\circ$ either x or y, and $v^{\phi}$ being the phase velocity (itself a function of frequency since waves on a plate are dispersive). Finally, `Nmodes` is the total number of modes to be computed by the eigenvalue routine. The largest mode number is `Nmax = (Nx+1)*(Ny+1)`, but in most cases one is interested in just a handful of modes given by `Nmodes` -- note that considerable speedups are so obtained. 

Plate parameters, rather self-explanatory, come next:

```
%--------------------
%-- plate parameters
rho = 390 ; % density [kg / m^3]
nux = 0.39 ; % poisson ratio
Ex = 10.4e9 ; % young's mod along x [Pa]
Ey = 0.994e9 ; % young's mod along y [Pa]
Gxy = 0.526e9 ; % shear mod [Pa]
Lx = 0.4 ; % edge lenght x [m]
Ly = 0.6 ; % edge length y [m]
Lz = 3e-3 ; % thickness [m]

%-- elastic constants around the edges
KRmat = [1e13,1e13; %Kx0 Rx0 
    1e13, 1e13; % K0y R0y
    1e13, 1e13; % KxL RxL
    1e13, 1e13] ; % KLy RLy
%--------------------
```
Note that, in order to simulate an "infinite" stiffness at the boundary, a rather large value must be selected. Usually, this can be set as 1e13. 

Next, braces and lumped load parameters must be set:

```
%--------------------
%-- braces parameters
Nribs = 8 ; % number of braces 

Eb = [10e9,10e9,10e9,10e9,10e9,10e9,10e9,10e9].' ; % youngs moduli
Lzb = [3e-3,3e-3,3e-3,3e-3,3e-3,3e-3,3e-3,3e-3].' ; % thicknesses 
bb = [3e-2,3e-2,3e-2,3e-2,3e-2,3e-2,3e-2,3e-2].' ; % width cross section
rhob = [400,400,400,400,400,400,400,400].' ; % densities

% rib coordinates along x (start and end) AS A FRACTION OF Lx
x_beam_coord = ...
    [0.2,0.8;
    0.2,0.8;
    0.2,0.8;
    0.2,0.8;
    0.2,0.8;
    0.2,0.8;
    0.2,0.8;
    0.2,0.8] ;

% rib coordinates along y (start and end) AS A FRACTION OF Ly
y_beam_coord = ...
    [0.1,0.2;
    0.2,0.3;
    0.3,0.4;
    0.4,0.5;
    0.5,0.6;
    0.6,0.7;
    0.7,0.8;
    0.8,0.9] ;
%--------------------

%--------------------
%-- static loads and stiffeners parameters

Nlump = 3 ; % number of lumped elements
x_lump_coord = [0.12,0.47,0.91].' ; % x coordinates of lumped elements 
y_lump_coord = [0.75,0.4,0.38].' ; % y coordinates of lumped elements 
Mlump  = [0.5,0.01,0.01].' ;
%--------------------
```

Again, these are rather self-explanatory. Note that the lengths of the rib arrays must be consistent with `Nribs`. So, the array containing the rib Young's moduli, `Eb`, must be of length `Nribs`. All rib parameters are given in SI units (Pa for Young's moduli, m for the thickness and width of the braces' cross sections, and Kg/m^3 for the densities). The ribs' x and y coordinates are specified in two matrices, `x_beam_coord` and `y_beam_coord`, of size `(Nribs) X 2`. Each brace's initial and final x and y coordinates are given as fractions of `Lx` and `Ly`. Similar parameters are set for the static loads. Note that `Mlump` contains the loads' masses (in Kg). 

Finally, the plot section completes the set of the user-editable parameters. These include the choice of the colormap and the flag `absPlot` to plot the absolute value of the modal shapes. 

```
%--------------------
%-- plot parameters parameters
cmap = cmaps(4) ; % select colormap 1 = RedBlue, 2 = GreenPurple, 3 = OrangeGreen, 4 = PurpleOrange
NN = 9 ; % first mode number to be plotted 
Nplots = 9 ; % select 3,6 or 9. If another number is selected, it is defaulted to 3. Displayed plots are NN + (0:Nplots - 1)
%--------------------
```
### Time Domain Analysis

The file [TimeDomainLauncher.m](https://github.com/Nemus-Project/magpie-matlab/blob/30c0e8bfd3d4026de2bd9c634764ddf6f7481469/orthotropic/src/TimeDomainLauncher.m) allows setting the parameters for the time domain simulations. Most of these are analogous to the frequency domain case above. The time domain simulation requires setting extra parameters controlling the input/output locations, the forcing type and parameters, and the decay parameters. 

The I/O setup is straightforward: the input `x_f` is a one-by-two vector including the x and y coordinates (again, scaled by the side lengths). The output `outMat` is similar, but it can contain an arbitrary amount of rows each representing an virtual contact microphone on the plate's surface.  

```
%--------------------
%-- input / output locations
x_f = [0.513,0.678] ; % frac of Lx Ly 
outMat = [0.51,0.52; 0.12,0.76] ; % frac of Lx Ly, each row represents an output point
%--------------------
```

The forcing parameters consist of a flag `forceType` and a one-by-three vector. Currently, MAGPIE supports two force types: a raised cosine simulating an impulse and a sinusoid. When the raised cosine is selected, the three parameters appearing in the `forceParams` array correspond to contact duration (in seconds), the largest forcing amplitude (in Newtons) and a noise modulation parameter (this can be useful to give the resulting synthesis a brighter tone).  In sinusoidal mode, the three forcing parameters are the sinusoidal input frequency, peak amplitude, and noise modulation. 

```
%--------------------
%-- forcing
forceType = 1 ; % 1 = impulse , 2 = sinusoid
if forceType == 1 
forceParams = [0.0007,50,0.5] ; % [time of contact (s), max Amplitude (N), noise modulation]
else
    forceParams = [100,50,0.2] ; % [frequency of sine (Hz), max Amplitude (N), noise modulation]
end
%--------------------
```

The decay-time parameters are expressed via another one-by-three array. Currently, MAGPIE handles a Rayleigh-type decay profile, such that the damping matrix is 

$$C = \alpha M + \beta K$$

where $M$ and $K$ are the mass and stiffness matrices of the system. Two free positive parameters $\alpha$ and $\beta$ are available to the user. In the frequency domain, Rayleigh's damping is quadratic in frequency, with $\alpha$ controlling the decay time at DC and $\beta$ multiplying the quadratic term (there is no linear term). Rather than setting the two parameters as are, it is more convenient to set the decay times $\tau_{60}$ at two known frequencies. One is fixed at DC, and the second is user-selectable. These compose the three parameters in the `dampVec` array.

```
%--------------------
%-- loss parameters
dampVec = [0.9, 0.4, 500] ; %[t60_0Hz (s),t60_f1 (s), f1 (Hz)] NB you MUST use t60_f1 < t60_0Hz for stability
%--------------------
```

Finally, there are extra plot parameters compared to the frequency domain. The `LivePlot` flag allows plotting the animations in real-time. If selected, six plots are updated at a rate corresponding to `Refreshrate`. Finally, `FilmRec` allows recording to file the resulting animation. 

```
%--------------------
%-- plot parameters parameters
cmap = cmaps(4) ; % select colormap 1 = RedBlue, 2 = GreenPurple, 3 = OrangeGreen, 4 = PurpleOrange
LivePlot = 1 ; % 1 : live plot on
RefreshRate = 1 ; % 1 = play all frames, 2 = play one out of two frames, etc
absPlot = 0 ; % 1 = plots absolute value (colormap will adjust accordingly)
FilmRec = 0 ; % 1= record video to file
%--------------------
```

## Inverse Modelling

MAGPIE allows the estimation of the elastic constants $E_x$, $E_y$, and $G_{xy}$ of a real orthotropic plate starting from a set of $N_{meas}$ measured frequencies and modal shapes. Note: **the $N_{meas}$ modal shapes must be known along with the corresponding frequencies**. The reference script is [InverseModelling.m](https://github.com/Nemus-Project/magpie-matlab/blob/e715bbf5cf5e27d54e0c0f392357c9037e2b2cd0/orthotropic/src/InverseModelling.m). The top of the script allows setting the plate parameters. 

```
%--------------------
%-- plate parameters
rho      = 473.9 ;
Lx       = 0.223 ;
Ly       = 0.114 ;
Lz       = 0.003 ;

%-- ballpark values 
Ex0      = 10.7e9 ;   
Ey0      = 716e6 ;
Gxy0     = 500e6 ;
nux0     = 0.51 ;

%-- elastic constants around the edges
KRmat = [0e13,0e13; %Kx0 Rx0
    1e13, 1e13; % K0y R0y
    0e13, 0e13; % KxL RxL
    0e13, 0e13] ; % KLy RLy

%-- measured experimental frequencies (Hz)
ExpFreqs = [52
    98
    311
    337
    398
    637] ;
%--------------------
```

The plate parameters here are drawn from the experimental plate. Of course, the values of the elastic constants are unknown at this stage. 
**The only requirement here is that the first $N_{meas}$ modal shapes correspond to the numerical eigenshapes returned by MAGPIE using the ballpark elastic constant values**. Checking that the experimental boundary conditions are correctly implemented is also fundamental at this stage (if one or more sides are clamped, you must ensure that sufficient pressure is applied on the plate's edges by the clamps). The experimentalist must ensure that these requirements are met before estimating the elastic constants via MAGPIE. 

Note that the $N_{meas}$ experimental frequencies corresponding to the first $N_{meas}$ mode shapes appear here in the preamble. 

The rest of the preamble is self-explanatory. 

```
%--------------------
%-- general parameters
fmax = 2000 ; %-- maximum frequency to be computed
ppw  = 100 ; % points per wavelength at maximum frequency. Choose 3 <= ppw
Nmeas = 6 ; % total number of modes required
%--------------------

%--------------------
%-- braces parameters
Nribs = 0 ; % number of braces
Eb = [].' ; % youngs moduli
Lzb = [].' ; % thicknesses
bb = [].' ; % width cross section
rhob = [].' ; % densities

% rib coordinates along x (start and end) AS A FRACTION OF Lx
x_beam_coord = ...
    [] ;

% rib coordinates along y (start and end) AS A FRACTION OF Ly
y_beam_coord = ...
    [] ;

%--------------------
%-- static loads and stiffeners parameters

Nlump = 0 ; % number of lumped elements
x_lump_coord = [].' ; % x coordinates of lumped elements
y_lump_coord = [].' ; % y coordinates of lumped elements
KLump  = [].' ;
MLump  = [].' ;
%--------------------

%--------------------
%-- plot parameters parameters
cmap = cmaps(1) ; % select colormap 1 = RedBlue, 2 = GreenPurple, 3 = OrangeGreen, 4 = PurpleOrange
absPlot = 0 ;
NN = 1 ; % mode number to plot (plot will display modes between NN and NN + 8)
%--------------------

%--------------------
%-- MAC parameters
MAC_threshold = 0.99 ;
%--------------------

% END CUSTOM PARAMETERS
%-------------------------------------------------------------------------
```

The only new parameter to set here is `MAC_threshold`. MAGPIE performs modal identification by computing the MAC (modal assurance criterion). A batch of "training" numerical plates is computed, and the modal shapes are assessed against the reference modal shapes computed using the ballpark elastic constants. When the average MAC across all $N_{meas}$ is below `MAC_threshold`, the training plate is excluded. This ensures that the training is operated only on a batch of plates returning consistent modal shapes. 
