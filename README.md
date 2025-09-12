# Structure-Preserving High-Order Methods for the Compressible Euler Equations in Potential Temperature Formulation for Atmospheric Flows

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17106781.svg)](https://zenodo.org/doi/10.5281/zenodo.17106781)

This repository contains information and code to reproduce the results presented in the article

```bibtex
@online{artiano2025structure,
  title={Structure-Preserving High-Order Methods for the
         Compressible {E}uler Equations in Potential Temperature 
         Formulation for Atmospheric Flows},
  author={Artiano, Marco and Knoth, Oswald and Spichtinger, Peter 
          and Ranocha, Hendrik},
  year={2025},
  month={09},
  eprint={TODO},
  eprinttype={arxiv},
  eprintclass={TODO}
}
```

If you find these results useful, please cite the article mentioned above. If you use the implementations provided here, please also cite this repository as

```bibtex
@misc{artiano2025structureRepo,
  title={Reproducibility repository for
         "{S}tructure-Preserving High-Order Methods for the
         Compressible {E}uler Equations in Potential Temperature 
         Formulation for Atmospheric Flows"},
  author={Artiano, Marco and Knoth, Oswald and Spichtinger, Peter 
            and Ranocha, Hendrik},
  year={2025},
  howpublished={\url{https://github.com/MarcoArtiano/2025_structure_potential_temperature}},
  doi={10.5281/zenodo.17106781}
}
```

## Abstract
We develop structure-preserving numerical methods for the compressible Euler equations, employing potential temperature as a prognostic variable.We construct three numerical fluxes designed to ensure the conservation of entropy and total energy within the discontinuous Galerkin framework on general curvilinear meshes.Furthermore, we introduce a generalization for the kinetic energy preservation property and total energy conservation in the presence of a gravitational potential term. To this end, we adopt a flux-differencing approach for the discretization of the source term, treated as non-conservative product. We present well-balanced schemes for different constant background states for both formulations (total energy and potential temperature) on curvilinear meshes. Finally, we validate the methods by comparing the potential temperature formulation with the traditional Euler equations formulation across a range of classical atmospheric scenarios.

## Numerical experiments
To reproduce the numerical experiments presented in this article, you need to install Julia. The numerical experiments presented in this article were performed using Julia v1.11.6.

First, you need to download this repository, e.g., by cloning it with git or by downloading an archive via the GitHub interface. Then, you need to start Julia in the code directory of this repository and follow the instructions described in the README.md file therein.

## Authors
- Marco Artiano
- Oswald Knoth
- Peter Spichtinger
- [Hendrik Ranocha](https://ranocha.de/) (Johannes Gutenberg University Mainz, Germany)

## License
The code in this repository is published under the MIT license, see the LICENSE file.

## Disclaimer
Everything is provided as is and without warranty. Use at your own risk!
