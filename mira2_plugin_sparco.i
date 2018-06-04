/*
 * mira2_plugin_sparco.i
 *
 * Implement the SPARCO plugin to mira2
 *
 */

//MIRA_PLUGDIR = "~/apps/src/mira2_plugin_sparco";


func mira_plugin_sparco_init(nil) {
  /*DOCUMENT mira_plugin_sparco_init(nil);

  Initialisation of the plugin from the parameters read in the
  command line.
  */

  inform, "Loading \"sparco\" plugin...";

  SPARCO_options = _lst("\nSparco specific options",
    _lst("sparco_model", "star", "NAME", OPT_STRING,
         "Name of the SPARCO model used"),
    _lst("sparco_params", [], "VALUES", OPT_REAL_LIST,
         "Parameters used with SPARCO"),
    _lst("sparco_w0", [], "VALUE", OPT_REAL,
         "Central wavelength (in microns) for SPARCO"),
    _lst("sparco_image", [], "NAME", OPT_STRING,
         "Name of the fits file that will be used for SPARCO/IMAGE"),
    _lst("sparco_star_type", [], "NAME", OPT_STRING,
         "Type of spectrum that the central star has"),
    _lst("sparco_star_file", [], "NAME", OPT_STRING,
         "Name of the ascii file where the stellar spectrum is"),
    _lst("sparco_star_index", [], "VALUE", OPT_REAL,
         "Stellar spectral index"),
    _lst("sparco_star_temp", [], "VALUE", OPT_REAL,
         "Stellar temperature for blackbody")
         );

  plugin = mira_new_plugin(options = SPARCO_options,
                parse_options = parse_options,
                tweak_visibilities=tweak_visibilities,
                tweak_gradient=tweak_gradient,
                add_keywords=add_keywords,
                add_extensions=add_extensions
                  );

  return plugin;
}

func parse_options(plugin, opt)
/* DOCUMENT parse_options(plugin, opt)

  plugin is a hash_table with following arguments:

  plugin.model = name of the model in SPARCO. It can be:
                       "star", "binary", "UD".
  plugin.params = vector of necessary parameters for the
                         models
  plugin.w0 = central wavelenghts for computation of chromaticity

  plugin.image = a fits-file with the image that will be used as a
                        SPARCO model

*/

{
  local sparco, params, w0, image;

  h_set, opt, flags=opt.flags | MIRA_KEEP_WAVELENGTH;

  /* Modify options if needed, for example xmin xmax
  h_set, opt, an_option=my_needed_value;

  /* Read SPARCO settings */
  sparco = opt.sparco_model;
  params = opt.sparco_params;
  w0 = opt.sparco_w0;
  image = opt.sparco_image;
  startype = opt.sparco_star_type;
  starfile = opt.sparco_star_file;
  startemp = opt.sparco_star_temp;
  starindex = opt.sparco_star_index;

  if ( !is_void(w0) ) {
    w0 *= 1e-6;
  };

  if (!is_void(sparco)) {
    if (sparco=="star") {
      if (numberof(params)!=2 ) {
        opt_error, "sparco/star takes a 2-parameters vector for ´-sparco_params´: [fs,denv].";
      };
      if (params(1) >=1. | params(1) <0.) {
        opt_error, "uncorrect value for the first sparco parameter fs. Should be: 0 <= fs < 1";
      }
    } else if (sparco == "binary") {
      if (numberof(params)!=5 | params(1) >=1 | params(1) <0 | params(3) >=1 | params(3) <0 | params(1) + params(3) >= 1) {
        opt_error, "sparco/binary takes a 5-parameters vector for ´-sparco_params´: [fs,denv,fbin,xbin,ybin].";
      };
    } else if (sparco == "UD") {
      if (numberof(params)!=3 | params(1) >=1 | params(1) <0) {
        opt_error, "sparco/UD takes a 3-parameters vector for ´-sparco_params´: [fs,denv,UD].";
      };
    } else if (sparco == "image") {
      opt_error, "not implemented yet...";
    } else if (sparco == "imageBB") {
      opt_error, "not tested and validated";
      if (numberof(params)!=3 | params(1) >=1 | params(1) <0 | params(2)<=0 | params(3)<=0 ) {
        opt_error, "sparco/UD takes a 3-parameters vector for ´-sparco_params´: [fim0,T0,Timg].";
      } else if ( !is_string(image) ) {
        opt_error, "sparco/imageBB takes an string with the name of the fits file";
      } else {
        /* Read initial image. */
        img = mira_read_image(image);
        naxis = img.naxis;
        if (naxis < 2 || naxis > 3) {
          opt_error, "Expecting a 2D/3D sparco model image";
        }
        naxis1 = img.naxis1;
        naxis2 = img.naxis2;
        naxis3 = (naxis >= 3 ? img.naxis3 : 1);
        if (naxis == 3) {
          if (naxis3 != 1) {
            warn, "Converting 3D image into a 2D grayscaled image";
          }
          h_set, img, arr = img.arr(,,avg), naxis=2;
          h_pop, img, "naxis3";
          h_pop, img, "crpix3";
          h_pop, img, "crval3";
          h_pop, img, "cdelt3";
          h_pop, img, "cunit3";
          h_pop, img, "ctype3";
          naxis = 2;
        }

        /* Maybe resample the SPARCO image. */
        siunit = "radian"; // pixelsize and FOV are in SI units
        cdelt1 = mira_convert_units(img.cunit1, siunit)*img.cdelt1;
        cdelt2 = mira_convert_units(img.cunit2, siunit)*img.cdelt2;
        if (is_void(pixelsize)) {
          pixelsize = min(cdelt1, cdelt2);
        }
        if (is_void(dims)) {
          if (is_void(fov)) {
            fov1 = cdelt1*naxis1;
            fov2 = cdelt2*naxis2;
            fov = max(fov1, fov2);
            dim1 = lround(fov1/pixelsize);
            dim2 = lround(fov2/pixelsize);
          } else {
            dim1 = lround(fov/pixelsize);
            dim2 = dim1;
          }
          dims = [2, dim1, dim2];
        }
      cunit = "mas";
      cdelt = mira_convert_units(siunit, cunit)*pixelsize;
      crpix1 = (dim1/2) + 1;
      crpix2 = (dim2/2) + 1;
      crval = 0.0;
      img = mira_resample_image(img, pad=0, norm=1n,
                              naxis1=dim1,   naxis2=dim2,
                              crpix1=crpix1, crpix2=crpix2,
                              crval1=crval,  crval2=crval,
                              cdelt1=cdelt,  cdelt2=cdelt,
                              cunit1=cunit,  cunit2=cunit);
      eq_nocopy, image, img.arr;
    };
  } else {
    opt_error, "the " + sparco + " model is not implemented in sparco";
  };
} else {
    opt_error, "no sparco model specified ";
};

inform, "SPARCO will be run with model "+sparco;
inform, "SPARCO has these parameters "+pr1(params);
inform, "SPARCO will use this w0 "+pr1(w0);

if (!is_void(startype)) {
  if (startype=="BB" | startype=="blackbody" | startype=="bbody" | startype=="bb") {
    if ( !is_void(startemp) )  {
      if (startemp > 1000.)
      inform, "The star has a black body spectrum with T=%g K \n", startemp;
      startype = "BB";
    } else {
      throw, "If you want to use a blackbody spectrum for the star, you must specify a temperature by using the option `-sparco_star_temp=` "
    };
  } else if (startype=="spectrum") {
    throw, "FIXME: implement the loading of the ascii file";
  } else if (startype=="pow") {
    if (!is_void(starindex)) {
      inform, "The star has a spectral index of %g", starindex;
    }
  }
} else {
  startype = "pow";
  starindex = -4.;
  inform, "The star has a spectral index of %g", starindex;
}

h_set, plugin, model=sparco,
               params=params,
               w0=w0,
               image=image,
               startype=startype,
               starindex=starindex,
               startemp=startemp,
               starfile=starfile,
               starinit=1n;

}

func _get_spectrum(type, w, w0, index=, temp=, file=)
/* DOCUMENT _get_spectrum(type, w, w0, index=, temp=, file=);

  DO NOT USE OUTSIDE THE SPARCO PLUGIN

  Defines the spectrum of a component with respect of
  what is specified
  w    = the wavelengths at which the spectrum will be computed
         (meters)
  w0   = central wavelength (meters)
  type = "pow",  "BB",        "spectrum"
  index= powerlaw index
  temp = black body temperature (K)
  file = ascii file with the spectrum

*/
{
  local star;

  if (type == "pow") {
  star = (w/w0)^index;
} else if (plugin.startype == "BB") {
  star = _BB(w,temp)/_BB(w0,temp);
} else {
  throw, "Spectrum not implemented yet";
}

  return star
};

func tweak_visibilities (master, vis)
/* DOCUMENT tweak_complex_visibilities (master, vis);

  Takes the complex visibilities from the image and
  adds the selected sparco model to return complex
  visibilities of image + SPARCO model

*/
{
 plugin = mira_plugin(master);

 /* init of stellar spectrum */
 if (plugin.starinit) {
   w = mira_model_wave(master);
   w0 = plugin.w0;
   star = _get_spectrum(plugin.startype, w, w0,
      index=plugin.starindex, temp=plugin.startemp,
      file=plugin.starfile);
   h_set, plugin, star_spectrum = star,
                  starinit=0n;
 }


 if (plugin.model == "star") {
   vis = mira_sparco_star(master, vis);
 } else if (plugin.model == "binary") {
   vis = mira_sparco_binary(master, vis);
 } else if (plugin.model == "UD") {
   vis = mira_sparco_UD(master, vis);
 } else if (plugin.model == "starBB") {
   vis = mira_sparco_starBB(master, vis);
 } else if (plugin.model == "image") {
   throw, "not implemented yet...";
   vis = mira_sparco_image(master, vis);
 } else {
   throw, "do not understand the sparco model used: "+plugin.model;

 };

return vis;
}


func tweak_gradient (master, grd)
/* DOCUMENT tweak_complex_gradient (master, grd);

  Takes the complex gradient on the image and
  normalise it to take into account SPARCO

  SEE ALSO: mira_sparco_star, mira_sparco_UD, mira_sparco_binary.
*/
{
  plugin = mira_plugin(master);
  /* Gradient modification for SPARCO */
  if (plugin.model == "star" | plugin.model == "UD" | plugin.model == "image") {
    fs0 = plugin.params(1);
    denv = plugin.params(2);
    w = mira_model_wave(master);
    w0 = plugin.w0;

    fs = fs0 * plugin.star_spectrum;
    fd = (1.-fs0) * (w/w0)^denv;
    ftot = fs + fd;

    grd_re = grd(1,..);
    grd_im = grd(2,..);

    grd_re *= ftot / fd;
    grd_im *= ftot / fd;

    grd = [grd_re, grd_im];
    grd = transpose(grd);

  } else if ( plugin.model == "binary") {

    fs0 = plugin.params(1);
    denv = plugin.params(2);
    fbin0 = plugin.params(3);
    w = mira_model_wave(master);
    w0 = plugin.w0;

    fs = fs0 * plugin.star_spectrum;
    fbin = fbin0 * (w/w0)^-4;
    fd = (1-fs0-fbin0) * (w/w0)^denv;
    ftot = fs + fd + fbin;

    grd_re = grd(1,..);
    grd_im = grd(2,..);

    grd_re *= ftot / fd;
    grd_im *= ftot / fd;

    grd = [grd_re, grd_im];
    grd = transpose(grd);
  };

  return grd

}


func mira_sparco_star(master, vis)
  /* DOCUMENT mira_sparco_star(master, vis);

     Compute the total complex visibilities by adding a point source at
     the center of the image using the stellar-to-total flux ratio (fs0) and
     the spectral index of the environment (denv)
     fs = fs0 * (lambda/lambda0)^-4
     fd = (1-fs0) * (lambda/lambda0)^denv
     Vtot = fd*Vimg + fs
     Vtot /= fd +fs

     SEE ALSO: mira_sparco_UD, mira_sparco_binary.
   */
{
  local vis, vis_re, vis_im, vis_amp, vis_phi, fs0, denv;

  vis_re = vis(1,..);
  vis_im = vis(2,..);

  plugin = mira_plugin(master);
  params = plugin.params;


  fs0 = params(1);
  denv = params(2);
  w = mira_model_wave(master);
  w0 = plugin.w0;

  fs = fs0 * plugin.star_spectrum;
  fd = (1-fs0) * (w/w0)^denv;
  ftot = fs + fd;

  vis_re = vis_re * fd + fs;
  vis_im = vis_im * fd;

  vis_re /= ftot;
  vis_im /= ftot;

  vis = [vis_re, vis_im];
  vis = transpose(vis);

  return vis;
};

func mira_sparco_starBB(master, vis)
  /* DOCUMENT mira_sparco_starBB(master, vis);

     Compute the total complex visibilities by adding a point source at
     the center of the image using the stellar-to-total flux ratio (fs0) and
     the spectral index of the environment (denv)
     fs = star
     fd = fd0 * BB(w,T)/BB(w0,T)
     Vtot = fd*Vimg + fs
     Vtot /= fd +fs

     SEE ALSO: mira_sparco_UD, mira_sparco_binary.
   */
{
  local vis, vis_re, vis_im, vis_amp, vis_phi, fs0, denv;

  vis_re = vis(1,..);
  vis_im = vis(2,..);

  plugin = mira_plugin(master);
  params = plugin.params;


  fs0 = params(1);
  T = params(2);
  w = mira_model_wave(master);
  w0 = plugin.w0;

  fs = fs0 * plugin.star_spectrum;
  fd = (1-fs0) * _BB(T, w) / _BB(T, w0);
  ftot = fs + fd;

  vis_re = vis_re * fd + fs;
  vis_im = vis_im * fd;

  vis_re /= ftot;
  vis_im /= ftot;

  vis = [vis_re, vis_im];
  vis = transpose(vis);

  return vis;
};


func mira_sparco_binary(master, vis)
  /* DOCUMENT mira_sparco_binary(master, vis);

     Compute the total complex visibilities by adding two point sources at
     (one being at the center of the image) using the flux ratios (fs0, fbin0),
     the spectral index of the environment (denv) and the positions of the
     secondary (xbin amnd ybin).
     fs = fs0 * (lambda/lambda0)^-4
     fbin = fbin0 * (lambda/lambda0)^-4
     fd = (1-fs0) * (lambda/lambda0)^denv
     Vtot = fd*Vimg + fs + fbin*Vbin
     Vtot /= fd + fs + fbin

     SEE ALSO: mira_sparco_star, mira_sparco_binary.
   */
{
  local vis, vis_re, vis_im, vis_amp, vis_phi, fs0, denv;

  vis_re = vis(1,..);
  vis_im = vis(2,..);

  plugin = mira_plugin(master);

  fs0 = plugin.params(1);
  denv = plugin.params(2);
  fbin0 = plugin.params(3);
  xbin = plugin.params(4) * MIRA_MAS;
  ybin = plugin.params(5) * MIRA_MAS;
  w = mira_model_wave(master);
  w0 = plugin.w0;

  fs = fs0 * plugin.star_spectrum;
  fbin = fbin0 * (w/w0)^-4;
  fd = (1-fs0-fbin0) * (w/w0)^denv;
  ftot = fs + fd + fbin;

  u = mira_model_u(master) / w;
  v = mira_model_v(master) / w;

  Vbin_re = cos( -2*pi*(xbin*u + ybin*v) );
  Vbin_im = sin( -2*pi*(xbin*u + ybin*v) );

  vis_re = vis_re * fd + fs + fbin*Vbin_re;
  vis_im = vis_im * fd + fbin*Vbin_im;

  vis_re /= ftot;
  vis_im /= ftot;

  vis = [vis_re, vis_im];
  vis = transpose(vis);

  return vis;
};

func mira_sparco_UD(master, vis)
  /* DOCUMENT mira_sparco_UD(master, vis);

     Compute the total complex visibilities by adding a Uniform Disk at
     the center of the image using the stellar-to-total flux ratio (fs0) and
     the spectral index of the environment (denv)
     fs = fs0 * (lambda/lambda0)^-4
     fd = (1-fs0) * (lambda/lambda0)^denv
     Vtot = fd*Vimg + fs*V_UD
     Vtot /= fd +fs

     SEE ALSO: mira_sparco_star, mira_sparco_binary.
   */
{
  local vis, vis_re, vis_im, vis_amp, vis_phi, fs0, denv, B, UD, w;

  vis_re = vis(1,..);
  vis_im = vis(2,..);

  plugin = mira_plugin(master);

  fs0 = plugin.params(1);
  denv = plugin.params(2);
  UD = plugin.params(3) * MIRA_MAS;
  w = mira_model_wave(master);
  w0 = plugin.w0;

  fs = fs0 * plugin.star_spectrum;
  fd = (1-fs0) * (w/w0)^denv;
  ftot = fs + fd;

  u = mira_model_u(master) / w;
  v = mira_model_v(master) / w;
  B = abs(u,v);

  V_UD = 2*bessj1(pi * B * UD) / ( pi * B * UD);

  vis_re = vis_re * fd + fs * V_UD;
  vis_im = vis_im * fd;

  vis_re /= ftot;
  vis_im /= ftot;

  vis = [vis_re, vis_im];
  vis = transpose(vis);

//  h_set, master, model_vis_re = vis_re, model_vis_im = vis_im, model_vis = vis;

  return vis;
};

func mira_sparco_imageBB(master, vis)
  /* DOCUMENT mira_sparco_imageBB(master, vis);

     Compute the total complex visibilities by adding a predefined image (im0)
     the reconstructed image using the im0-to-total flux ratio (fim0) and
     the Black Body temperatures of the two images (T0, Timg)
     fim0 = fim0 * BB(T0, lambda) / BB(T0, lam0)
     fd = (1-fim0) * BB(Timg, lambda) / BB(Timg, lam0)
     Vtot = fd*Vimg + fim0*Vimg0
     Vtot /= fd +fim0

     SEE ALSO: mira_sparco_star, mira_sparco_binary.
   */
{
  local vis, vis_re, vis_im, vis_amp, vis_phi, fs0, denv, B;

  vis_re = vis(1,..);
  vis_im = vis(2,..);

  plugin = mira_plugin(master);

  fim0 = plugin.params(1);
  T0 = plugin.params(2);
  Tim = plugin.params(3);
  w = mira_model_wave(master);
  w0 = plugin.w0;
  img0 = plugin.image;
  BB = []; //TODO Find a function for blackbody


  fim = fim0 * plugin.star_spectrum;
  fd = (1-fim0) * _BB(Tim, w) / _BB(Tim, w0);
  ftot = fim + fd;

  Vimg0 = master.xform(img0);
  Vimg0_re = Vimg0(1,..);
  Vimg0_im = Vimg0(2,..);

  vis_re = vis_re * fd + fim * Vim0_re;
  vis_im = vis_im * fd + fim * Vim0_im;

  vis_re /= ftot;
  vis_im /= ftot;

  vis = [vis_re, vis_im];
  vis = transpose(vis);

//  h_set, master, model_vis_re = vis_re, model_vis_im = vis_im, model_vis = vis;

  return vis;
};

func add_keywords (master, fh)
/* DOCUMENT add_keywords (master, fh);

    This function adds the sparco plugin add_keywords
    to the finale saved fits file.

    SEE ALSO: add_extension.
*/
{
  plugin = mira_plugin(master);

  fits_set, fh, "SMODEL",  plugin.model,  "Model used in SPARCO";
  fits_set, fh, "SWAVE0",  plugin.w0,  "Central wavelength (mum) for chromatism";
  fits_set, fh, "STTYPE",  plugin.startype,  "Spectral model for star";
  if (plugin.startype == "BB") {
    fits_set, fh, "STTEMP",  plugin.startemp,  "Temperature of the star";
  } else if (plugin.startype == "pow") {
    fits_set, fh, "STINDE",  plugin.starindex,  "Spectral index of the star";
  }

}

func add_extension (master, fh)
/* DOCUMENT add_extension (master, fh);

    This function adds the sparco plugin add_keywords
    to the finale saved fits file.

    SEE ALSO: add_keywords.
*/
{
  plugin = mira_plugin(master);


 //FIXME: add stellar spectrum if used by sparco
 // fits_new_hdu, fh, "IMAGE", "SPARCO adds an extension";
 // fits_write_array, fh, random(128,128);
 // fits_pad_hdu, fh;

}

func _BB(lambda, T)
    /* DOCUMENT _BB(lambda, T);

       DESCRIPTION
       Computatiof black body radiation function with a given temperature T
       at a specific wavelength lambda.
       (Inspired from yocoAstroBBodyLambda)

       PARAMETERS
       - lambda : wavelength (m)
       - T      : temperature (K)

       RETURN VALUES
       Return the energy radiated (W/m2/m)
    */
{
    local c, h, kb;
    c = 2.99792453e8;
    h = 6.626070040e-34;
    kb = 1.38064852e-23;
    mask = abs(h*c / (kb*T*lambda) ) > 700;
    if(numberof(where(mask))==0)
        flambda = 2*h*c^2 / lambda^5 / (exp(h*c / (kb*T*lambda)) - 1);
    else
    {
        flambda = lambda;
        flambda(where(mask)) = 0.0;
        flambda(where(!mask)) = 2*h*c^2 / lambda^5 / (exp(h*c / (kb*T*lambda)) - 1);
    }

    return flambda;
}
