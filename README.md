# EchelleCCFs [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://RvSpectML.github.io/EchelleCCFs.jl/stable) [![Build Status](https://github.com/RvSpectML/EchelleCCFs.jl/workflows/CI/badge.svg)](https://github.com/RvSpectML/EchelleCCFs.jl/actions) [![Coverage](https://codecov.io/gh/RvSpectML/EchelleCCFs.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/RvSpectML/EchelleCCFs.jl)


The EchelleCCFs package is designed to compute Cross Correlation Functions of stellar spectra from high-resolution Echelle spectragraphs.  It is part of the [RvSpectML family](https://rvspectml.github.io/RvSpectML-Overview/) of Julia packages.  
Examples of using it with other RvSpectML packages can be found in the examples and/or notebooks folders of the [RvSpectML](https://github.com/eford/RvSpectML.jl) and [RvSpectMLPlots](https://github.com/RvSpectML/RvSpectMLPlots.jl) packages.  

## Examples
The EchelleCCFs package includes a few simplistic tests and one example of applying it to simulated solar spectra from [Gilbertson, Ford & Dumusque 2020](https://doi.org/10.5281/zenodo.3753253) which is described in the associated [Research Note of the AAS](https://ui.adsabs.harvard.edu/link_gateway/2020RNAAS...4...59G/doi:10.3847/2515-5172/ab8d44).  
The input file for this example is large, so it is not part of the GitHub repository, but is downloaded after installing the pacakge.
To install the package, download the simulated data and perform a CCF calculation, you can run

```julia
import Pkg
Pkg.add("EchelleCCFs")  # Installs EchelleCCFs and required dependancies
using EchelleCCFs       # Loads EchelleCCFs
cd(joinpath(pkgdir(EchelleCCFs),"examples"))  # Change into examples folder
Pkg.activate(".")      # Specify to use Project.toml in examples folder
Pkg.instantiate()      # Installs extra packages for example
include("../deps/download.jl")  # Force download of simulated dataset
include("soap_ccf.jl") # Run example CCF calculation and make plots
```

Good luck.
