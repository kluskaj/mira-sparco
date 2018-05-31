/*
 * mira-sparco.i
 *
 * Implement the SPARCO plugin to mira2
 *
 */


func parse_arguments (plugin, opt)
/* DOCUMENT parse_arguments(plugin, opt)
*/

{


}


func tweak_complex_visibilities (master, vis)
/* DOCUMENT tweak_complex_visibilities (master, vis);

  Takes the complex visibilities from the image and
  adds the selected sparco model to return complex
  visibilities of image + SPARCO model

*/
{


}


func tweak_complex_gradient (master, grd)
/* DOCUMENT tweak_complex_gradient (master, grd);

  Takes the complex gradient on the image and
  normalise it to take into account SPARCO

*/
{


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


  eq_nocopy, vis, mira_get_model_vis(master);
  eq_nocopy, vis_re, mira_get_model_vis_re(master);
  eq_nocopy, vis_im, mira_get_model_vis_im(master);

  fs0 = master.params(1);
  denv = master.params(2);
  w = master.w;
  w0 = master.w0;
  fs = fs0 * (w/w0)^-4;
  fd = (1-fs0) * (w/w0)^denv;
  ftot = fs + fd;

  vis_re = vis_re * fd + fs;
  vis_im = vis_im * fd;

  vis_re /= ftot;
  vis_im /= ftot;

  vis = [vis_re, vis_im];
  vis = transpose(vis);

  h_set, master, model_vis_re = vis_re, model_vis_im = vis_im, model_vis = vis;

  return master;
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


  eq_nocopy, vis, mira_get_model_vis(master);
  eq_nocopy, vis_re, mira_get_model_vis_re(master);
  eq_nocopy, vis_im, mira_get_model_vis_im(master);

  fs0 = master.params(1);
  denv = master.params(2);
  fbin = master.params(3);
  xbin = master.params(4) * MIRA_MAS;
  ybin = master.params(5) * MIRA_MAS;
  w = master.w;
  w0 = master.w0;
  fs = fs0 * (w/w0)^-4;
  fbin = fbin0 * (w/w0)^-4;
  fd = (1-fs0-fbin0) * (w/w0)^denv;
  ftot = fs + fd + fbin;

  u = master.u;
  v = master.v;

  Vbin_re = cos( -2*pi*(xbin*u + ybin*v) );
  Vbin_im = sin( -2*pi*(xbin*u + ybin*v) );

  vis_re = vis_re * fd + fs + fbin*Vbin_re;
  vis_im = vis_im * fd + fbin*Vbin_im;

  vis_re /= ftot;
  vis_im /= ftot;

  vis = [vis_re, vis_im];
  vis = transpose(vis);

  h_set, master, model_vis_re = vis_re, model_vis_im = vis_im, model_vis = vis;

  return master;
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

  vis_re = vis_re(1,..);
  vis_im = vis_re(2,..);

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
