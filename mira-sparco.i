/*
 * mira-sparco.i
 *
 * Implement the SPARCO plugin to mira2
 *
 */


func parse_arguments (plugin, opt)
/* DOCUMENT parse_arguments(plugin, opt)

  plugin is a hash_table with following arguments:

  plugin.sparco = name of the model in SPARCO. It can be:
                       "star", "binary", "UD".
  plugin.sparco_params = vector of necessary parameters for the
                         models
  plugin.sparco_w0 = central wavelenghts for computation of chromaticity
  plugin.sparco_func = function to linearly add the SPARCO model to
                            complex visibilities from the image.

  plugin.sparco_image = to be implemented...

  Creation of plugin hash_table in this function ?
*/

{

  /* Read SPARCO settings */
  sparco = opt.sparco;
  params = opt.params;

  if (!is_void(sparco)) {
    if ( !is_void(w0) ) {
      w0 *= 1e-6;
    };
    if (is_void(sparco) ) {
      opt_error, "no sparco model specified ";
    } else if (sparco=="star") {
      if (numberof(params)!=2 | params(1) >=1. | params(1) <0.) {
        opt_error, "sparco/star takes a 2-parameters vector for ´-params´: [fs,denv].";
      };
    } else if (sparco == "binary") {
      if (numberof(params)!=5 | params(1) >=1 | params(1) <0 | params(3) >=1 | params(3) <0 | params(1) + params(3) >= 1) {
        opt_error, "sparco/binary takes a 5-parameters vector for ´-params´: [fs,denv,fbin,xbin,ybin].";
      };
    } else if (sparco == "UD") {
      if (numberof(params)!=3 | params(1) >=1 | params(1) <0) {
        opt_error, "sparco/UD takes a 3-parameters vector for ´-params´: [fs,denv,UD].";
      };
    } else if (sparco == "image") {
      opt_error, "not implemented yet...";
    } else {
      opt_error, "the " + sparco + " model is not implemented in sparco";
    };
  };

h_set, plugin, w0=w0;

}


func tweak_complex_visibilities (master, vis)
/* DOCUMENT tweak_complex_visibilities (master, vis);

  Takes the complex visibilities from the image and
  adds the selected sparco model to return complex
  visibilities of image + SPARCO model

*/
{

 if (master.plugin.sparco == "star") {
   vis = mira_sparco_star(master, vis);
 } else if (master.plugin.sparco == "binary") {
   vis = mira_sparco_binary(master, vis);
 } else if (master.plugin.sparco == "UD") {
   vis = mira_sparco_UD(master, vis);
 } else if (master.plugin.sparco == "image") {
   throw, "not implemented yet...";
   vis = mira_sparco_image(master, vis);
 };

return vis;
}


func tweak_complex_gradient (master, grd)
/* DOCUMENT tweak_complex_gradient (master, grd);

  Takes the complex gradient on the image and
  normalise it to take into account SPARCO

*/
{

  /* Gradient modification for SPARCO */
  if (master.plugin.sparco == "star" | master.plugin.sparco == "UD" | master.plugin.sparco == "image") {

    fs0 = master.plugin.params(1);
    denv = master.plugin.params(2);
    w = master.coords.unique.w;
    w0 = master.plugin.w0;

    fs = fs0 * (w/w0)^-4;
    fd = (1-fs0) * (w/w0)^denv;
    ftot = fs + fd;

    grd_re = grd(1,..);
    grd_im = grd(2,..);

    grd_re *= ftot / fd;
    grd_im *= ftot / fd;

    grd = [grd_re, grd_im];
    grd = transpose(grd);

  } else if ( master.plugin.sparco == "binary") {

    fs0 = master.plugin.params(1);
    denv = master.plugin.params(2);
    fbin0 = master.plugin.params(3);
    w = master.coords.unique.w;
    w0 = master.plugin.w0;

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


func mira_sparco_star(master)
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

  fs0 = master.plugin.params(1);
  denv = master.plugin.params(2);
  w = master.coord.unique.w;
  w0 = master.plugin.w0;

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

func mira_sparco_binary(master)
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

  fs0 = master.plugin.params(1);
  denv = master.plugin.params(2);
  fbin = master.plugin.params(3);
  xbin = master.plugin.params(4) * MIRA_MAS;
  ybin = master.plugin.params(5) * MIRA_MAS;
  w = master.coords.unique.w;
  w0 = master.plugin.w0;

  fs = fs0 * (w/w0)^-4;
  fbin = fbin0 * (w/w0)^-4;
  fd = (1-fs0-fbin0) * (w/w0)^denv;
  ftot = fs + fd + fbin;

  u = master.coords.unique.u / w;
  v = master.coords.unique.v / w;

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

  fs0 = master.plugin.params(1);
  denv = master.plugin.params(2);
  UD = master.plugin.params(3) * MIRA_MAS;
  w = master.coords.unique.w;
  w0 = master.plugin.w0;

  fs = fs0 * (w/w0)^-4;
  fd = (1-fs0) * (w/w0)^denv;
  ftot = fs + fd;

  u = master.coords.unique.u / w;
  v = master.coords.unique.v / w;
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
