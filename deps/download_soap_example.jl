# Should be run from EchelleCCFs/deps or EchelleCCFs/examples or EchelleCCFs/test
# Since this May be called before EchelleCCFs is installed, so this doesn't use pkgdir.

include("download.jl")

download_url = "https://psu.box.com/shared/static/nikva44r7xai9h0ulfj2pxewzv8ez4fo.h5"
download_filename = joinpath("..","data","spectra","soap_demo.h5")
download_md5 = "18a6291039d7684ef9187762e4b2dfd5"

download_and_check_md5sum(download_url, download_filename, md5_goal=download_md5)
