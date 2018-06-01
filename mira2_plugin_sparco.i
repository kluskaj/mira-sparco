/*
 * mira2_plugin_sparco.i
 *
 * Implement the SPARCO plugin to mira2
 *
 */

// MIRA_PLUGDIR = "~/apps/src/mira2_plugin_sparco";


func mira_plugin_sparco_init(nil) {

  inform, "Loading \"sparco\" plugin...";

  SPARCO_options = _lst("\nSparco specific options",
    _lst("sparco_model", "star", "NAME", OPT_STRING,
         "Name of the SPARCO model used"),
    _lst("sparco_params", [], "VALUES", OPT_REAL_LIST,
         "Parameters used with SPARCO"),
    _lst("sparco_w0", [], "VALUE", OPT_REAL,
         "Central wavelength (in microns) for SPARCO"),
    _lst("sparco_image", [], "NAME", OPT_STRING,
         "Name of the fits file that will be used for SPARCO/IMAGE")
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
  plugin.func = function to linearly add the SPARCO model to
                            complex visibilities from the image.

  plugin.image = a fits-file with the image that will be used as a
                        SPARCO model

  Creation of plugin hash_table in this function ?
*/

{
  local sparco, params, w0, image;

  // h_set, opt, flags=opt.flags | MIRA_KEEP_WAVELENGTH;

  /* Read SPARCO settings */
  sparco = opt.sparco_model;
  params = opt.sparco_params;
  w0 = opt.sparco_w0;
  image = opt.sparco_image;

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

h_set, plugin, model=sparco,
               params=params,
               w0=w0,
               image=image;

}


func tweak_visibilities (master, vis)
/* DOCUMENT tweak_complex_visibilities (master, vis);

  Takes the complex visibilities from the image and
  adds the selected sparco model to return complex
  visibilities of image + SPARCO model

*/
{
 plugin = mira_plugin(master);

 if (plugin.model == "star") {
   vis = mira_sparco_star(master, vis);
 } else if (plugin.model == "binary") {
   vis = mira_sparco_binary(master, vis);
 } else if (plugin.model == "UD") {
   vis = mira_sparco_UD(master, vis);
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

    fs = fs0 * (w/w0)^-4.;
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

    fs = fs0 * (w/w0)^-4;
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
  /* DOCUMENT mira_sparco_star(master);

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

  fs = fs0 * (w/w0)^-4;
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

func mira_sparco_binary(master, vis)
  /* DOCUMENT mira_sparco_binary(master);

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
  fbin = plugin.params(3);
  xbin = plugin.params(4) * MIRA_MAS;
  ybin = plugin.params(5) * MIRA_MAS;
  w = mira_model_wave(master);
  w0 = plugin.w0;

  fs = fs0 * (w/w0)^-4;
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
  /* DOCUMENT mira_sparco_UD(master);

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

  fs = fs0 * (w/w0)^-4;
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
  /* DOCUMENT mira_sparco_imageBB(master);

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


  fim = fim0 * BB(T0, w) / BB(T0, w0);
  fd = (1-fim0) * BB(Tim, w) / BB(Tim, w0);
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
/*

*/
{
  plugin = mira_plugin(master);

  fits_set, fh, "SMODEL",  plugin.model,  "Model used in SPARCO";
  fits_set, fh, "SWAVE0",  plugin.w0,  "Central wavelength (mum) for chromatism";

}

func add_extension (master, fh)
/*

*/
{
  plugin = mira_plugin(master);

  fits_new_hdu, fh, "IMAGE", "SPARCO adds an extension";

}
