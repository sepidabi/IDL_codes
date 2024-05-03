pro calib_data_greg

  invdir = '/home/seki2695/INV/stic/example_me/'
  invdir_test =  '/home/seki2695/INV/stic/example/'
  pxn = 2

  restore, '/scratch/gviss/tmp/sepideh/inv/ca8_cak_fe1_20160915.idlsave'

  fe_me = readfits(invdir+'testf_feobs.fits') & fe = fe_me
  ca8_me = readfits(invdir+'testf_ca8obs.fits') & ca8 = ca8_me
  cak_me = readfits(invdir+'testf_cakobs.fits') & cak = cak_me

  fe_test = readfits(invdir_test+'feobs.fits')
  ca8_test = readfits(invdir_test+'ca8obs.fits')
  cak_test = readfits(invdir_test+'cakobs.fits')


  fe[*,0] = fe1_wav
  fe[*,1:4] = fe1_dat[0,pxn,*,*]

  ca8[*,0] = ca8_wav
  ca8[*,1:4] = ca8_dat[0,pxn,*,*]

  cak[0:20,0] = cak_wav[0:20]
  cak[0:20,1] = cak_dat[0,pxn,0:20]

  cak[-1,0] = cak_wav[-1]
  cak[-1,1] = cak_dat[0,pxn,-1]

  cgwindow
  cgplot, fe[*,0], fe[*,1], color = 'blue', /add;,psym = 16
  cgplot, fe_me[*,0], fe_me[*,1], color = 'red', /add,/over;,psym = 16
  cgplot, fe_test[*,0], fe_test[*,1], color = 'green', /add,/over;,psym = 16
  cglegend, title = ['Gregal','example','me'], color = ['blue','green','red'], /add;,psym = 16

  cgwindow
  cgplot, ca8[*,0], ca8[*,1], color = 'blue', /add;,psym = 16
  cgplot, ca8_me[*,0], ca8_me[*,1], color = 'red', /add,/over;,psym = 16
  cgplot, ca8_test[*,0], ca8_test[*,1], color = 'green', /add,/over;,psym = 16
  cglegend, title = ['Gregal','example','me'], color = ['blue','green','red'], /add,psym = 16

  cgwindow
  cgplot, cak[0:20,0], cak[0:20,1], color = 'blue', /add;,psym = 16
  cgplot, cak_me[*,0], cak_me[*,1], color = 'red', /add,/over;,psym = 16
  cgplot, cak_test[*,0], cak_test[*,1], color = 'green', /add,/over;,psym = 16
  cglegend, title = ['Gregal','example','me'], color = ['blue','green','red'], /add;,psym = 16

  writefits, invdir+'testf_feobs_g.fits', fe
  writefits, invdir+'testf_ca8obs_g.fits', ca8
  writefits, invdir+'testf_cakobs_g.fits', cak
  stop
end
