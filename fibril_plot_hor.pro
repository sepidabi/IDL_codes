pro fibril_plot_hor

  cgdelete,/all

  ;TO BE DEFINED
  w_h = 7                       ;h alpha line core position
  w_h_cont = 0
  w_ca = 7
  w_ca_cont = 42
  fn = 01                       ;desired frame

  cak_cube = 'crispex_3950_2016-09-15T08:49:53_scans=24-57_time-corrected.fcube'
  caksp_cube = 'crispex_3950_2016-09-15T08:49:53_scans=24-57_time-corrected_sp.fcube'
  cakint_cube = 'cak_int.fcube'
  cakintun_cube = 'cak_int_un.fcube'

  hsp_cube = 'crispex.6563.08:49:51.time_corrected_CHROMIS_sp.fcube'
  hint_cube = 'h_int.fcube'
  hintun_cube = 'h_int_un.fcube'
  
  datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
  fibdir = 'fr'+strtrim(long(fn),2)+'/'
  savedir = '/scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/'
  file_ca = file_search(datadir+fibdir+'*3950*.csav')
  file_h = file_search(datadir+fibdir+'*6563*.csav')
  res = 0.0375

  ;Image Dimensions
  lp_header, datadir+cak_cube, nx = nxx, ny = nyy ;image dimensions
  lp_header, datadir+caksp_cube, nx = nw_ca, ny = ntt ;frame number and Ca K wavelength
  lp_header,  datadir+hsp_cube, nx = nw_h             ;H alpha wavelength

  Ic_cak = lp_get(datadir+cak_cube,fn*nw_ca+w_ca_cont)
  Ic_cak_m = max(Ic_cak)
  Ic_h = lp_get(datadir+cak_cube,fn*nw_h+w_h_cont)
  Ic_h_m = max(Ic_h)
;stop
  cgwindow, wxsize = 650, wysize = n_elements(file_ca)*300/2
  
for n = 0L, n_elements(file_ca)/2-1 do begin

  ;fibril

  ;Ca K  
  restore, file_ca(2*n)
  F_ca = loop_slab
  F_ca_x = x_coords
  F_ca_y = y_coords
  F_ca_l = loop_size
  ca_scl = scaling_factor

  ;H alpha
  restore, file_h(2*n)
  F_h = loop_slab
  F_h_x = x_coords
  F_h_y = y_coords
  F_h_l = loop_size

  
  ;fibril background

  ;Ca K
  restore, file_ca(2*n+1)
  B_ca = loop_slab
  B_ca_x = x_coords
  B_ca_y = y_coords
  B_ca_l = loop_size

  ;H alpha
  restore, file_h(2*n+1)
  B_h = loop_slab
  B_h_x = x_coords
  B_h_y = y_coords
  B_h_l = loop_size
  h_scl = scaling_factor
  
  ;Ca K
  Ftot_ca =fltarr(21)
  Btot_ca = fltarr(21)
  Fmean_ca=fltarr(21)
  Bmean_ca = fltarr(21)
  
  for i = 0L, F_ca_l-1 do begin
     Ftot_ca(*) = Ftot_ca(*) + F_ca[i,fn,0:20]
  endfor
  for j = 0L, B_ca_l-1 do begin
     Btot_ca(*) = Btot_ca(*) + B_ca[j,fn,0:20]
  endfor
  
  Fmean_ca = Ftot_ca/F_ca_l
  Bmean_ca = Btot_ca/B_ca_l
  
  restore,datadir+'spectfile.3950.idlsave'
  ca_sp_pos = spect_pos(0:20)   ;Ca K spectral positions
  
  ;H alpha
  Ftot_h =fltarr(nw_h)
  Btot_h = fltarr(nw_h)
  Fmean_h=fltarr(nw_h)
  Bmean_h = fltarr(nw_h)

  for i = 0L, F_h_l-1 do begin
     Ftot_h(*) = Ftot_h(*) + F_h[i,fn,0:nw_h-1]
  endfor
  
  for j = 0L, B_h_l-1 do begin
     Btot_h(*) = Btot_h(*) + B_h[j,fn,0:nw_h-1]
  endfor
  Fmean_h = Ftot_h/F_h_l
  Bmean_h = Btot_h/B_h_l
  
  restore,datadir+'spectfile.6563.idlsave'
  h_sp_pos = spect_pos(0:nw_h-1) ;H K spectral positions

 
  ;plot averaged signal vs. wavelength
  pos = cgLayout([2,n_elements(file_ca)/2], OXMargin=[7,10], OYMargin=[0.5,0.5], XGap=1.,ygap=0.1)
  
  if n eq n_elements(file_ca)/2-1 then begin
       ;Ca K
     cgplot,ca_sp_pos-3934,Bmean_ca*10e8,color = 'dodger blue',/addcmd,xrange=[min(ca_sp_pos)-3934,max(ca_sp_pos)-3934],position = pos(*,2*n),/noerase,thick = 0.5,xthick = 0.5,ythick = 0.5,ytitle = 'I$\downmean$ (nW m$\up-3$ sr$\up-1$)', yrange = [1,6], xtitle = '$\lambda$ - $\lambda$$\downCa K$ ($\Angstrom$)'
     cgplot,ca_sp_pos-3934,Fmean_ca*10e8,color = 'red',/addcm,/overplot,position = pos(*,2*n),thick = 0.5
     
       ;H alpha
     cgplot,h_sp_pos-6563,Bmean_h/(h_scl),color = 'dodger blue',/addcmd,xrange=[min(h_sp_pos)-6563,max(h_sp_pos)-6563],position = pos(*,2*n+1),/noerase,thick = 0.5,xthick = 0.5,ythick = 0.5,ytickformat="(A1)",yrange = [0.1,1],ystyle = 8, xtitle = '$\lambda$ - $\lambda$$\downHa$ ($\Angstrom$)'
     cgAxis, YAxis=1,yrange=[0.1,1],/window,ythick = 0.5,ytitle = 'Scaled I$\downmean$'
     cgplot,h_sp_pos-6563,Fmean_h/(h_scl),color = 'red',/addcm,position = pos(*,2*n+1),/noerase,/overplot,thick = 0.5

  endif else begin
       ;Ca K
       cgplot,ca_sp_pos-3934,Bmean_ca*10e8,color = 'dodger blue',/addcmd,xrange=[min(ca_sp_pos)-3934,max(ca_sp_pos)-3934],position = pos(*,2*n),/noerase,thick = 0.5,xthick = 0.5,ythick = 0.5,xtickformat="(A1)",ytitle = 'I$\downmean$ (nW m$\up-3$ sr$\up-1$)', yrange = [1,6]
       cgplot,ca_sp_pos-3934,Fmean_ca*10e8,color = 'red',/addcm,/overplot,position = pos(*,2*n),thick = 0.5
  
       ;H alpha
       cgplot,h_sp_pos-6563,Bmean_h/(h_scl),color = 'dodge rblue',/addcmd,xrange=[min(h_sp_pos)-6563,max(h_sp_pos)-6563],position = pos(*,2*n+1),/noerase,thick = 0.5,xthick = 0.5,ythick = 0.5,yrange = [0.1,1],xtickformat="(A1)",ytickformat="(A1)",ystyle = 8
       cgAxis, YAxis=1,yrange=[0.1,1],/window,ythick = 0.5,ytitle = 'Scaled I$\downmean$'
       cgplot,h_sp_pos-6563,Fmean_h/(h_scl),color = 'red',/addcm,position = pos(*,2*n+1),/noerase,/overplot,thick = 0.5
    endelse
  
endfor
;stop
cgcontrol, output = savedir+'w_plot_'+strtrim(fn,2)+'.eps',/ps_encap,ps_char = 0.1,resize = [800,n_elements(file_ca)*300/2];,/adjustsi
cgdelete

stop
end
