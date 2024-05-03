pro fibril_path_part

  cgdelete,/all
  ;TO BE DEFINED
  w_h = 7                       ;h alpha line core position
  w_ca = 7
  fn = 01                       ;desired frame

  ;coords of the specific region/fibril
  x = [1010,1055]
  y = [1025,1100]

  n = 36                        ;fibril no.

  nxx = x(1)-x(0)
  nyy = y(1)-y(0)

  cak_cube = 'crispex_3950_2016-09-15T08:49:53_scans=24-57_time-corrected.fcube'
  caksp_cube = 'crispex_3950_2016-09-15T08:49:53_scans=24-57_time-corrected_sp.fcube'
  cakint_cube = 'cak_int.fcube'
  cakintun_cube = 'cak_int_un.fcube'
  hsp_cube = 'crispex.6563.08:49:51.time_corrected_CHROMIS_sp.fcube'
  hint_cube = 'h_int.fcube'
  hintun_cube = 'h_int_un.fcube'
  
  datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
  savedir = '/scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/'
  fibdir = 'fr'+strtrim(long(fn),2)+'/'
  file_ca = file_search(datadir+fibdir+'*3950*.csav')
  file_h = file_search(datadir+fibdir+'*6563*.csav')
  res = 0.0375

  ;Image Dimensions
;  lp_header, datadir+cak_cube, nx = nxx, ny = nyy ;image dimensions
;  lp_header, datadir+caksp_cube, nx = nw_ca, ny = ntt ;frame number and Ca K wavelength
;  lp_header,  datadir+hsp_cube, nx = nw_h ;H alpha wavelength

  img_ca = lp_get(datadir+cakintun_cube,fn)
  img_h = lp_get(datadir+hintun_cube,fn)

  ;display fibril paths on the data
  cgwindow, wxsize =2*nxx,wysize = nyy;,wmulti = [0,2,1] ;,woxmargin = [1,0],woymargin = [4,0.5]
  pos1 = cgLayout([2,1], OXMargin=[4,0], OYMargin=[10,2], XGap=0)
  axis_format1 = {xticklen:-0.02,yticklen:-0.02,xtitle:'[arcsec]',ytitle:'[arcsec]'}
  cgimage,img_ca(x[0]:x[1],y[0]:y[1]),/scale,/keep_,/axes,/addcmd,axkeywords = axis_format1,xrange = [0,nxx*res],yrange = [0,nyy*res],position = pos1[*,0]

;stop
;for n = 0L, n_elements(file_ca)/2-1 do begin

  ;fibril

  ;Ca K  
  restore, file_ca(n)
  F_ca = loop_slab
  F_ca_x = x_coords
  F_ca_y = y_coords
  F_ca_l = loop_size

  ;H alpha
  restore, file_h(n)
  F_h = loop_slab
  F_h_x = x_coords
  F_h_y = y_coords
  F_h_l = loop_size

  
  ;fibril background

  ;Ca K
  restore, file_ca(n+1)
  B_ca = loop_slab
  B_ca_x = x_coords
  B_ca_y = y_coords
  B_ca_l = loop_size
;stop
  cgplot, F_ca_x*res,F_ca_y*res,/addcmd,/noerase,color = 'red',thick = 2.,position = pos1[*,0],/overplot,aspect=2./3;, linestyle = 2
  cgplot, B_ca_x*res,B_ca_y*res,/addcmd,/noerase,color = 'dodger blue',thick = 0.9,position = pos1[*,0],/overplot

;  cgtext,/addcmd,/data,n, F_ca_x(0)*res, F_ca_y(0)*res, color = 'red', charsize = 0.5
;  print,n
;  stop
;endfor

  ;display background paths on the data
  axis_format2 = {xticklen:-0.02,yticklen:-0.02,xtitle:'[arcsec]',ytickformat:"(A1)"}
  cgimage, img_h(x[0]:x[1],y[0]:y[1]),/scale,/keep_,/axes,/addcmd,axkeywords = axis_format2,xrange = [0,nxx*res],yrange = [0,nyy*res],position = pos1[*,1]

;stop
;extracting each fibril spectral info
;for n = 0L, n_elements(file_ca)/2-1 do begin

  ;fibril

  ;Ca K  
  restore, file_ca(n)
  F_ca = loop_slab
  F_ca_x = x_coords
  F_ca_y = y_coords
  F_ca_l = loop_size

  ;H alpha
  restore, file_h(n)
  F_h = loop_slab
  F_h_x = x_coords
  F_h_y = y_coords
  F_h_l = loop_size

  
  ;fibril background

  ;Ca K
  restore, file_ca(n+1)
  B_ca = loop_slab
  B_ca_x = x_coords
  B_ca_y = y_coords
  B_ca_l = loop_size

  ;H alpha
  restore, file_h(n+1)
  B_h = loop_slab
  B_ca_l = loop_size

  cgoplot, F_ca_x*res,F_ca_y*res,/addcmd,/noerase,color = 'red',thick = 2.,position = pos1[1,1];, linestyle = 2
  cgoplot, B_ca_x*res,B_ca_y*res,/addcmd,/noerase,color = 'dodger blue',thick = 0.9,position = pos1[1,1]

;endfor

  ;save the path map
  cgcontrol,output = savedir+'path_'+strtrim(fn,2)+'_un_f'+strtrim(n,2)+'.eps',/ps_encap, ps_char = 0.7,/adjustsize
  cgdelete
  
  stop
end
