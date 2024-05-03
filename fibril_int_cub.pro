pro fibril_int_cub

  cak = 'crispex_3950_2016-09-15T08:49:53_scans=24-57_time-corrected.fcube'
  cak_sp = 'crispex_3950_2016-09-15T08:49:53_scans=24-57_time-corrected_sp.fcube'
  h = 'crispex.6563.08:49:51.time_corrected_CHROMIS.fcube'
  h_sp = 'crispex.6563.08:49:51.time_corrected_CHROMIS_sp.fcube'

  lp_header, cak, nx = nx,ny = ny
  lp_header, cak_sp, nx = nw, ny = nt

  lp_header, h_sp, nx = nwh

  rad = 10.
  cc = 4.

  ;t = 20
  cub = lp_read(cak)
  cubh = lp_read(h)
  int_cub = fltarr(nx,ny,nt)
  int_cub_un = fltarr(nx,ny,nt)
  int_cubh = fltarr(nx,ny,nt)
  int_cubh_un = fltarr(nx,ny,nt)
;  window,xsize = nx/cc,ysize = ny/cc,0
;  window,xsize = nx/cc,ysize = ny/cc,1
  for t = 0L, nt-1 do begin
     for i = 0L,8 do begin
        int_cub[*,*,t] = int_cub[*,*,t]+cub[*,*,t*nw+(i+6)]
     endfor

     for j = 0L, 6 do begin
        int_cubh[*,*,t] = int_cubh[*,*,t]+cubh[*,*,t*nwh+(j+4)]
     endfor

     int_cub_un[*,*,t] = unsharp_mask(int_cub[*,*,t], radius = rad)
     int_cubh_un[*,*,t] = unsharp_mask(int_cubh[*,*,t], radius = rad)
     
;     wset,0
 ;    tvscl,congrid(int_cub[*,*,t],nx/cc,ny/cc)
 ;    wset,1
 ;    tvscl,congrid(int_cub_un[*,*,t],nx/cc,ny/cc)
     ;stop
  endfor
;  stop

  ;Creating the integrated and unsharped cubes
  ;===========================================
  lp_write,int_cub, 'cak_int.fcube'
  lp_write,int_cub_un, 'cak_int_un.fcube'

  lp_write,int_cubh, 'h_int.fcube'
  lp_write,int_cubh_un, 'h_int_un.fcube'

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

  
