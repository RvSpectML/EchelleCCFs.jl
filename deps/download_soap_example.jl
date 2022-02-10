# Should be run from EchelleCCFs/deps or EchelleCCFs/examples or EchelleCCFs/test
# Since this May be called before EchelleCCFs is installed, so this doesn't use pkgdir.

include("download.jl")

# download_url = "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/ebf11_psu_edu/EYyYhOHEjM5RqvZn1D-4I58BDbpbcYg6QZKhd0KQ3N7YoQ?e=Is4PcG"
download_url = "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/ebf11_psu_edu/EYyYhOHEjM5RqvZn1D-4I58BDbpbcYg6QZKhd0KQ3N7YoQ?download=1"
download_filename = joinpath("..","data","spectra","soap_demo.h5")
download_md5 = "18a6291039d7684ef9187762e4b2dfd5"

download_and_check_md5sum(download_url, download_filename, md5_goal=download_md5)
