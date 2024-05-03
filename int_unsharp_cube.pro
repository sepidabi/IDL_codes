pro int_unsharp_cube

  ; DATA info
  datadir = '/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/'
  
  ; integration range
  wi_cak = 9 & wf_cak = 17
  wi_ca8 = 6 & wf_ca8 =16
  
  ; smoothing factors
  rad = 10. ; radius of smoothness
  cc = 4.   ; plotting scale factor
  
  cak = datadir+'crispex_3950_2018-07-22T08:23:58_scans=0-212_time-corrected.fcube'
  cak_sp = datadir+'crispex_3950_2018-07-22T08:23:58_scans=0-212_time-corrected_sp.fcube'
  ;h = 'crispex.6563.08:49:51.time_corrected_CHROMIS.fcube'
  ;h_sp = 'crispex.6563.08:49:51.time_corrected_CHROMIS_sp.fcube'
  ca8 = datadir+'crispex.stokes.8542.08:23:57.time_corrected_CHROMIS.fcube'
  ca8_sp = datadir+'crispex.stokes.8542.08:23:57.time_corrected_CHROMIS_sp.fcube'

  mapdata,cak,cak_sp,cub_k,nx,ny,nt,nw_cak
  mapdatastokes,ca8,ca8_sp,cub_8,nx,ny,nt,nw_ca8

  int_cubcak = fltarr(nx,ny,nt)
  int_cubcak_un = fltarr(nx,ny,nt)
  ;int_cubh = fltarr(nx,ny,nt)
  ;int_cubh_un = fltarr(nx,ny,nt)
  int_cubca8 = fltarr(nx,ny,nt)
  int_cubca8_un = fltarr(nx,ny,nt)
;  window,xsize = nx/cc,ysize = ny/cc,0
;  window,xsize = nx/cc,ysize = ny/cc,1


  
  for t = 0L, nt-1 do begin

     int_cubcak[*,*,t] = total((cub_k[t])[*,*,wi_cak:wf_cak],3)

     int_cubca8[*,*,t] = total((cub_k[t])[*,*,wi_ca8:wf_ca8,0],3)
;stop

     int_cubcak_un[*,*,t] = unsharp_mask(int_cubcak[*,*,t], radius = rad)
     ;int_cubh_un[*,*,t] = unsharp_mask(int_cubh[*,*,t], radius = rad)
     int_cubca8_un[*,*,t] = unsharp_mask(int_cubca8[*,*,t], radius = rad)
     print, t
;     wset,0
 ;    tvscl,congrid(int_cubcak[*,*,t],nx/cc,ny/cc)
 ;    wset,1
 ;    tvscl,congrid(int_cubcak_un[*,*,t],nx/cc,ny/cc)
     ;stop
  endfor

  ;Creating the integrated and unsharped cubes
  ;===========================================
  lp_write,int_cubcak, datadir+'cak_int.fcube'
  lp_write,int_cubcak_un, datadir+'cak_int_un.fcube'

  ;lp_write,int_cubh, 'h_int.fcube'
  ;lp_write,int_cubh_un, 'h_int_un.fcube'

  lp_write,int_cubca8, datadir+'ca8_int.fcube'
  lp_write,int_cubca8_un, datadir+'ca8_int_un.fcube'

stop
  ;General spect profile
  ;=====================
  restore,'spectfile.3950.idlsave'
  cgwindow
  cgplot,spect_pos[0:20]-3934,norm_spect[0:20]/max(norm_spect[0:20]),yrange = [min(norm_spect[0:20]/max(norm_spect[0:20])),1],xtitle ='$\lambda$-$\lambda$$\downCa,0$',ytitle = 'Normalised Intensity',/addc
  cgplot,spect_pos[6:14]-3934,norm_spect[6:14]/max(norm_spect[0:20]),color = 'red',/addcm,/over
  cgplot,/over,/addcm,psym = 16, color = 'red', [spect_pos[6],spect_pos[14]]-3934,[norm_spect[6],norm_spect[14]]/max(norm_spect[0:20])
  ;cgtext, /plac, /addcm,color = 'red', strtrim(spect_pos[6]-3934,2)
  cgcontrol,output = '/scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/cak_spect.eps',ps_char = 0.8,/ps_en,resize = [150,55]

  cgdelete

  restore,'spectfile.6563.idlsave'
  cgwindow
  cgplot,spect_pos[0:15]-6563,norm_spect[0:15]/max(norm_spect[0:15]),yrange = [min(norm_spect[0:15]/max(norm_spect[0:15])),1],xtitle ='$\lambda$-$\lambda$$\downHa,0$',ytitle = 'Normalised Intensity',/addc
  cgplot,spect_pos[4:10]-6563,norm_spect[4:10]/max(norm_spect[0:15]),/addc,/over,color = 'red'
  cgplot,/over,/addcm,psym = 16, color = 'red', [spect_pos[4],spect_pos[10]]-6563,[norm_spect[4],norm_spect[10]]/max(norm_spect[0:15])

  cgcontrol,output = '/scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/Ha_spect.eps',ps_char = 0.8,/ps_en,resize = [150,55]

  cgdelete

  stop

end

  
