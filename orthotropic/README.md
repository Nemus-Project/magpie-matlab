## About


MAGPIE is an open-source framework for simulating plate vibration with elastic boundary conditions. The plate equation is derived from the *orthotropic* Kirchhoff-Love plate theory. It can then effectively simulate the dynamics of moderately thin wooden and composite plates. The geometry is restricted to rectangular, and the thickness is assumed constant. These limitations may be overcome in future releases. MAGPIE can simulate boundary conditions of classic type (free, simply-supported and clamped) and generally elastic boundary support, controlling the flexural displacement and boundary rotation through applied elastic forces and moments. 

MAGPIE is a framework allowing both *direct* and *inverse modelling* of orthotropic sheet materials, such as wood and composites, displaying orthotropic behaviour. Note that the isotropic case can be recovered as a special case. 

** Direct modelling consists of supplying MAGPIE with the appropriate elastic parameters (two Young's moduli, one shear modulus and one Poisson ratio), material parameters (density) and geometric parameters (edge lengths and thickness) plus the values of the elastic boundary constants. The output of the direct modelling can be:

1. Eigenvalues and eigenvectors (i.e., modal shapes)
2. A time-domain simulation of a struck plate using a simple raised cosine (simulating an impact) or sinusoidal input 

** Inverse modelling consists of supplying MAGPIE with a number of measured modal frequencies from an experimental plate, along with the plate's known material and geometric constants (density, edge lengths, thickness) to compute the elastic constants.

## The Orthotropic Kirchhoff-Love Model

The vibration of an orthotropic plate can be described via the Kirchhoff-Love model. This dynamical model describes the time evolution of the flexural displacement $u = u({\bf x},t) : \mathcal{V} \times \mathbb{R}^+_0$. Here,  ${\bf x} \in \mathcal{V} := [0,L_x] \times [0,L_y]$ is the domain of definition, a rectangle with side lengths $L_x$. Here, and in what follows, $x$ denotes the longitudinal orthotropic direction, $y$ is radial, and $z$ is tangential. In quarter-sawing, $z$ is the direction along the thickness of the board. With this notation, the system reads:

$$\rho h \partial_t^2 w = -D_1\partial_x^4 w - (D_2+D_4)\partial_y^2\partial_x^2 w - D_3\partial_y^4 w$$

The model depends on four rigidity constants, denoted by $D_\circ$, defined in terms of Young's moduli $E_\circ$, the shear modulus $G_{xy}$ and the Poisson's ratios $\nu_\circ$ as follows:

$$ D_1 := \frac{E_xh^3}{12(1-\nu_{xy}\nu_{yx})} \quad D_3 := \frac{E_yh^3}{12(1-\nu_{xy}\nu_{yx})} \quad D_4 := \frac{G_{xy}h^3}{3} \quad D_2 := \frac{\nu_{yx}E_{x}h^3}{6(1-\nu_{xy}\nu_{yx})} $$

Note that, because of the symmetry of the compliance matrix, one out of five elastic constants is fixed. Usually, one selects

$$\nu_{yx} = \nu_{xy}E_{y}E_{x}^{-1}$$

thus leaving the following elastic constants to be specified at the input: $E_x,E_y,G_{xy},\nu_{xy}$. Note that the model also depends on the thickness $h$ (assumed constant throughout $\mathcal V$), the side lengths $L_x,L_y$ and the density $\rho$ (also assumed constant).


```
Bottom-left corner: 
                 |
                 |
 [K_{0y},R_{0y}] |
                 |
                 |
                 |
                 |
                  ---------------------
                      [K_{x0},R_{x0}]
```

The equation must be supplied with appropriate boundary conditions. Here, conditions of elastic-type are assumed. The boundary force balance is expressed as:



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
