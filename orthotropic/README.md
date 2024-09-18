# MAGPIE
> *An open-source plate vibration simulation and analysis package*

M. Ducceschi`*`, M. Hamilton`*`, A. Mousseau`^`, S. Duran`*`

`*` NEMUS Lab - University of Bologna

`^` LAUM - Laboratoire d'Acoustique de l'Université du Mans


## About

MAGPIE is an open-source framework for simulating plate vibration with elastic boundary conditions. The plate equation is derived from the *orthotropic* Kirchhoff-Love plate theory. It can then effectively simulate the dynamics of moderately thin wooden and composite plates. The geometry is restricted to rectangular, and the thickness is assumed constant. These limitations may be overcome in future releases. MAGPIE can simulate boundary conditions of classic type (free, simply-supported and clamped) and generally elastic boundary support, controlling the flexural displacement and boundary rotation through applied elastic forces and moments. 

MAGPIE is a framework allowing both *direct* and *inverse modelling* of orthotropic sheet materials, such as wood and composites, displaying orthotropic behaviour. Note that the isotropic case can be recovered as a special case. 

** Direct modelling consists of supplying MAGPIE with the appropriate elastic parameters (two Young's moduli, one shear modulus and one Poisson ratio), material parameters (density) and geometric parameters (edge lengths and thickness) plus the values of the elastic boundary constants. The output of the direct modelling can be:

1. Eigenvalues and eigenvectors (i.e., modal shapes)
2. A time-domain simulation of a struck plate using a simple raised cosine (simulating an impact) or sinusoidal input 

** Inverse modelling consists of supplying MAGPIE with a number of measured modal frequencies from an experimental plate, along with the plate's known material and geometric constants (density, edge lengths, thickness) to compute the elastic constants.

## The Orthotropic Kirchhoff-Love Model

The vibration of an orthotropic plate can be described via the Kirchhoff-Love model. This dynamical model describes the time evolution of the flexural displacement $u = u(x,y,t) : \mathcal{V} \times \mathbb{R}^+_0$. Here,  ${\bf x} \in \mathcal{V} := [0,L_x] \times [0,L_y]$ is the domain of definition, a rectangle with side lengths $L_x, L_y$. Here, and in what follows, $x$ denotes the longitudinal orthotropic direction, $y$ is radial, and $z$ is tangential. In quarter-sawing, $z$ is the direction along the thickness of the board. With this notation, the system reads:

$$\rho h \partial_t^2 u(x,y,t) = -D_1\partial_x^4 u(x,y,t) - (D_2+D_4)\partial_y^2\partial_x^2 u(x,y,t) - D_3\partial_y^4 u(x,y,t)$$

The model depends on four rigidity constants, denoted by $D_\circ$, defined in terms of Young's moduli $E_\circ$, the shear modulus $G_{xy}$ and the Poisson's ratios $\nu_\circ$ as follows:

$$ D_1 := \frac{E_xh^3}{12(1-\nu_{x}\nu_{y})} \quad D_3 := \frac{E_yh^3}{12(1-\nu_{x}\nu_{y})} \quad D_4 := \frac{G_{xy}h^3}{3} \quad D_2 := \frac{\nu_{yx}E_{x}h^3}{6(1-\nu_{y}\nu_{y})} $$

Note that, because of the symmetry of the compliance matrix, one out of five elastic constants is fixed. Usually, one selects

$$\nu_{y} = \nu_{x}E_{y}E_{x}^{-1}$$

thus leaving the following elastic constants to be specified at the input: $E_x,E_y,G_{xy},\nu_{xy}$. Note that the model also depends on the thickness $h$ (assumed constant throughout $\mathcal V$), the side lengths $L_x,L_y$ and the density $\rho$ (also assumed constant).


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


$$\rho h \partial_t^2 u_{l,m}(t) = -D_1\delta_{xxxx} u_{l,m}(t) - (D_2+D_4)\delta_{xx}\delta_{yy} u_{l,m} - D_3\delta_{yyyy} u_{l,m}(t)$$

Note that the spatial indices are integer numbers restricted to the intervals $0 \leq l \leq N_x = \frac{L_x}{h_x}$ and $0 \leq m \leq N_y = \frac{L_y}{h_y}$, yielding $N_x+1$ points along $x$, and $N_y+1$ points along $y$, including grid point on the boundary. 

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

with analogous discretisations holding for the boundary conditions along the other edges. 

<picture>
  <source srcset="orthotropic/img/BCs.png">
</picture>

![Cacca](/img/BCs.png?raw=true "Dio Cane")

![Alt text](relative%20path/to/img.jpg?raw=true "Title")

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="[magpie-matlab/img/plateBCs1.pdf](https://github.com/Nemus-Project/magpie-matlab/blob/b505f5980a970319211afdcf859e7b1dc65353b2/orthotropic/img/plateBCs1.pdf)">
</picture>
 
## Structure

All source code is found in the `src/` directory of the repository. In `src/` you can find a directory for each supported language as well as a `data/` directory which contains any datasets shared across implementations.

## MAGPIE function

The MAGPIE function accepts:

- plate material properties: density, elasticity (Young's modulus), Poisson ratio
- plate dimensions (rectangular geometry: edge lengths and thickness)
- elastic boundary constants (2 per edge

The user can stipulate an accuracy coefficient and how many modes they wish to calculate.

The output from `magpie` is:

- `Q`: A list of eigenvectors, one for each mode. This can be used to visualize the mode shape.
- `Om`: The angular frequency for the corresponding eigenvector.
- `N`:  Number of grid points in the $x$ and $y$ dimensions for the plate. A smaller `h` value will result in a larger number of grid points
- `biharm`: The biharmonic used for deriving eigenvectors and modal frequencies. The biharmonic can also be used for a finite difference difference time domain scheme of the plate.

## References

- Howard, & Angus, J. A. S. (2009). _Acoustics and psychoacoustics_ (4th ed..).
- Ashby, (October 2021) [_Material property data for engineering materials_](https://www.ansys.com/content/dam/amp/2021/august/webpage-requests/education-resources-dam-upload-batch-2/material-property-data-for-eng-materials-BOKENGEN21.pdf) Department of Engineering, University of Cambridge (5th edition)
