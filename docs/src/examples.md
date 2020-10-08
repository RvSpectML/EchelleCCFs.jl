# Example Data Files

The EchelleCCFs package includes a few tests and one example of applying it to simulated solar spectra from [Gilbertson, Ford & Dumusque 2020](https://doi.org/10.5281/zenodo.3753253) which is described in the associated [Research Note of the AAS](https://ui.adsabs.harvard.edu/link_gateway/2020RNAAS...4...59G/doi:10.3847/2515-5172/ab8d44).  
The input file for this example is large, so it is not part of the GitHub repository, but is downloaded as part of the build process for the the package.
This is the first spectrum from
`res-1000-1years_full_id1.h5` (with the wavelengths converted from air to vacuumb).
To install the package, and perform an example CCF calculation, you can run

```julia
import Pkg
Pkg.add("EchelleCCFs")  # Installs EchelleCCFs and required dependancies
using EchelleCCFs       # Loads EchelleCCFs
cd(joinpath(pkgdir(EchelleCCFs),"examples"))  # Change into examples folder
Pkg.activate(".")      # Specify to use Project.toml in examples folder
Pkg.instantiate()      # Installs extra packages for example
include("soap_ccf.jl") # Run example CCF calculation and make plots
```

To analyze more spectra from the [Gilbertson, Ford & Dumusque 2020 dataset](https://doi.org/10.5281/zenodo.3753253), you'll need to download even larger file(s), as illustrated in `deps/download_soap_example_from_gilbertson_etal_2020.jl`.
