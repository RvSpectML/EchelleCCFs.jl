# Should be run from EchelleCCFs/deps or EchelleCCFs/examples or EchelleCCFs/test
# Since this May be called before EchelleCCFs is installed, so this doesn't use pkgdir.

include("download.jl")

download_url = "https://zenodo.org/record/3753254/files/res-1000-1years_full_id1.h5?download=1"
download_filename = joinpath("..","data","spectra","res-1000-1years_full_id1.h5")
download_md5 = "5659082144cd093d617bb54dca937ad9"

@warn "This is a large download.  Be prepared to be patient."
download_and_check_md5sum(download_url, download_filename, md5_goal=download_md5)
