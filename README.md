# mira-sparco: the SPARCO plugin for the MiRA image reconstruction algorithm

SPARCO (Semi-Parametric Approach for Reconstruction of Chromatic Objects) is
an approach to reconstruct chromatic images from optical interferometric data
by adding a geometrical model to the reconstructed image.
This plugin is designed to work with [MiRA2](https://github.com/emmt/MiRA) directly from the command line.

## Installation

Either install it by copying it to MIRA_HOME or use the option
`-plugin=/the/path/to/sparco/mira2_plugin_sparco.i`

TODO

## Usage

add the following options to mira2:

`-sparco_model=`      the model you want to use (See Sect. Implemented Models)

`-sparco_params=`     a vector of parameters corresponding to your sparco model (e.g. `-params=0.4,1` )

`-sparco_w0=`         central wavelengths for flux power laws computation for chromaticity (See Sect. Implemented Models)

`-sparco_star_type=`  Type of spectral behavior for the star. Can be either: "pow" for a power law, "BB" for a blackbody

`-sparco_star_index=` if `-sparco_star_type=pow` then the stellar spectral index needs to be specified (default: `-sparco_star_index=-4`)

`-sparco_star_temp=` if `-sparco_star_type=BB` then the stellar black body temperature needs to be specified (default: `-sparco_star_temp=10000`)

## Implemented Models

* **STAR:** It adds a point source at the center of the image using the stellar-to-total flux ratio (fs0) and the spectral index of the environment (denv)

fs = fs0 * stellar_spectrum

fd = (1-fs0) * (lambda/lambda0)^denv

Vtot = fd * Vimg + fs

Vtot /= fd +fs

* **BINARY:** It adds two point sources (one being at the center of the image) using the flux ratios (fs0, fbin0), the spectral index of the environment (denv) and the positions of the secondary (xbin amnd ybin).

fs = fs0 * stellar_spectrum

fbin = fbin0 * (lambda/lambda0)^-4

fd = (1-fs0) * (lambda/lambda0)^denv

Vtot = fd * Vimg + fs + fbin * Vbin

Vtot /= fd + fs + fbin

* **UD:** It adds a Uniform Disk at the center of the image using the stellar-to-total flux ratio (fs0) and the spectral index of the environment (denv)

fs = fs0 * stellar_spectrum

fd = (1-fs0) * (lambda/lambda0)^denv

Vtot = fd * Vimg + fs * V_UD

Vtot /= fd +fs


## Contact

If you have problems or you want to implement a specific model to be used with SPARCO,
please drop me an email at: jacques.kluska@kuleuven.be

## References

* Kluska, J. et al.: *"SPARCO : a semi-parametric approach for image reconstruction of chromatic objects. Application to young stellar objects"* in Astronomy & Astrophysics, Volume 564, id.A80, 11 pp. [DOI](https://ui.adsabs.harvard.edu/link_gateway/2014A&A...564A..80K/doi:10.1051/0004-6361/201322926)

## Examples of astrophysical results using SPARCO

* Kluska, J. et al.: *"A disk asymmetry in motion around the B[e] star MWC158"* in
    Astronomy & Astrophysics, Volume 591, id.A82, 15 pp. [DOI](https://ui.adsabs.harvard.edu/link_gateway/2016A&A...591A..82K/doi:10.1051/0004-6361/201527924)

* Hillen, M., Kluska, J. et al.: *"Imaging the dust sublimation front of a circumbinary disk"* in Astronomy & Astrophysics, Volume 588, id.L1, 6 pp. [DOI](https://ui.adsabs.harvard.edu/link_gateway/2016A&A...588L...1H/doi:10.1051/0004-6361/201628125)
