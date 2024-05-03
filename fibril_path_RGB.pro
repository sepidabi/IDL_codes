pro fibril_path_RGB

  cgdelete,/all
  ;TO BE DEFINED
  w_h = 7                       ;h alpha line core position
  w_ca = 7
  fn = 28                       ;desired frame
  ;crop crisp edges with [x1,x2,y1,y2]
  crop = [245,1678,10,1130]
  x1 = crop[0] & x2 = crop[1] & y1 = 20 & y2 = crop[3]


  cak_cube = 'crispex_3950_2016-09-15T08:49:53_scans=24-57_time-corrected.fcube'
  caksp_cube = 'crispex_3950_2016-09-15T08:49:53_scans=24-57_time-corrected_sp.fcube'
  cakint_cube = 'cak_int.fcube'
  cakintun_cube = 'cak_int_un.fcube'
  hsp_cube = 'crispex.6563.08:49:51.time_corrected_CHROMIS_sp.fcube'
  hint_cube = 'h_int.fcube'
  hintun_cube = 'h_int_un.fcube'
  ca8sp_cube = 'crispex.stokes.8542.08:49:51.time_corrected_CHROMIS_sp.fcube'
  ca8int_cube = 'ca8_int.fcube'
  ca8intun_cube = 'ca8_int_un.fcube'
  
  datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
  savedir = '/nadir-scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/'
  fibdir = 'fr'+strtrim(long(fn),2)+'/'
  file_ca = file_search(datadir+fibdir+'*3950*.csav')
  file_h = file_search(datadir+fibdir+'*6563*.csav')
  res = 0.0375

  ;Image Dimensions
  lp_header, datadir+cak_cube, nx = nxx, ny = nyy ;image dimensions
  lp_header, datadir+caksp_cube, nx = nw_ca, ny = ntt ;frame number and Ca K wavelength
  lp_header,  datadir+hsp_cube, nx = nw_h ;H alpha wavelength

  def = fltarr(nxx,nyy)
  def8 = fltarr(nxx,nyy)
  ;int, un Ca II K
  img_ca = cgGmaScl(lp_get(datadir+cakintun_cube,fn),min = 7.034e-9,max = 5.043e-8,gamma = 0.8)
  ;int, un H alpha
  img_h = cgGmaScl(lp_get(datadir+hintun_cube,fn),min = 7040.089,max = 26250.090,gamma = 0.8)
  ;int, un Ca 8542
  img_ca8 = cgGmaScl(lp_get(datadir+ca8intun_cube,fn),min = 350,max = 1200,gamma = 0.8)
  ;...cropping edges
  def(crop[0]:crop[1],crop[2]:crop[3]) = img_h(crop[0]:crop[1],crop[2]:crop[3])
  img_h = def
  def8(crop[0]:crop[1],crop[2]:crop[3]) = img_ca8(crop[0]:crop[1],crop[2]:crop[3])
  img_ca8 = def8

  ;stop

  ;display fibril paths on the data
  cgwindow;, wxsize =2*nxx,wysize = nyy+140;,wmulti = [0,2,1] ;,woxmargin = [1,0],woymargin = [4,0.5]
  pos1 = cgLayout([2,2], OXMargin=[5.5,1], OYMargin=[3.5,1], XGap=2, Ygap = 0)
  ;maps
  map = [img_ca8,img_ca,img_h,img_ca]
  ;axis formats
  ax0 = {xticklen:-0.02,yticklen:-0.02,xtitle:'',xtickformat:"(A1)",ytickformat:"(I0)",ytitle: 'y [arcsec]'}
  ax1 = {xticklen:-0.02,yticklen:-0.02,xtitle:'',xtickformat:"(A1)",ytickformat:"(A1)",ytitle: ''}
  ax2 = {xticklen:-0.02,yticklen:-0.02,xtitle:'x [arcsec]',xtickformat:"(I0)",ytickformat:"(I0)",ytitle:'y [arcsec]'}
  ax3 = {xticklen:-0.02,yticklen:-0.02,xtitle:'x [arcsec]',xtickformat:"(I0)",ytickformat:"(A1)",ytitle: ''}
  axis = [ax0, ax1, ax2, ax3]

  
  for p = 0L, n_elements(map(*,0))/nxx -1 do begin
   cgloadct,3, /addcmd
  ;plotting background maps
  cgimage,map(p*nxx:(p+1)*nxx-1,*),/scale,/keep_,/axes,/addcmd,axkeywords = axis(p),xrange = [0,nxx*res],yrange = [0,nyy*res],position = pos1[*,p]

  
  ;overplotting fibrils
  for n = 0L, n_elements(file_ca)/2-1 do begin

  ;fibril

  ;Ca K  
  restore, file_ca(2*n)
  F_ca = loop_slab
  F_ca_x = x_coords
  F_ca_y = y_coords
  F_ca_l = loop_size

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
;stop

  ;if p eq 2 then begin
   ;  cgplot, F_ca_x*res,F_ca_y*res,/addcmd,/noerase,color = 'red',thick = 1.5,position = pos1[*,0],/overplot,aspect=2./3
     ;cgplot, B_ca_x*res,B_ca_y*res,/addcmd,/noerase,color = 'gray',thick = 1.,position = pos1[*,0],/overplot , linestyle = 1
  ;endif
  
  if p eq 3 then begin
     cgplot, F_ca_x*res,F_ca_y*res,/addcmd,/noerase,color = 'cyan',thick = 1.5,position = pos1[*,0],/overplot,aspect=2./3 ;, linestyle = 2
     ;cgplot, B_ca_x*res,B_ca_y*res,/addcmd,/noerase,color = 'black',thick = 1.,position = pos1[*,0],/overplot, linestyle = 1
  endif

  if p eq 1 then begin
     if n eq 95 then begin
        ;cgplot, F_ca_x*res,F_ca_y*res,/addcmd,/noerase,color = 'red',thick = 1.25,position = pos1[*,0],/overplot,aspect=2./3, linestyle = 1
     endif
  endif
  
         

  
;  cgtext,/addcmd,/data,n, F_ca_x(0)*res, F_ca_y(0)*res, color = 'red', charsize = 0.5
;  print,n
;  stop
endfor
  
  ;caption = ['Ca II K $\lambda$-integrated & unsharp-masked','H$\alpha$ $\lambda$-integrated & unsharp-masked']
  caption = ['Ca II K','H$\alpha$', 'Ca II 8542 $\AA$']

  if p eq 2 then cgtext,pos1(0,p)+1, pos1(1,p)+1, 'c) H$\alpha$', /add, color = 'white'
  if p eq 1 then cgtext,pos1(0,p)+1, pos1(1,p)+1, 'b) Ca II K', /add, color = 'white'
  if p eq 0 then cgtext,pos1(0,p)+1, pos1(1,p)+1, 'a) Ca II 8542 $\Angstrom$', /add, color = 'white'
  if p eq 3 then cgtext,pos1(0,p)+1, pos1(1,p)+1, 'd) Ca II K', /add, color = 'white'

  

  if p eq 1 then begin

         ; the ROI box
         xmin = 1064. & xmax = 1353.
         ymin = 487. & ymax = 639.
         
         ;cgplot, [xmin,xmax]*res, [ymin,ymin]*res, /addcmd,/noerase,color = 'white',thick = 0.9,position = pos1[*,p],/overplot, linestyle = 0
         ;cgplot, [xmin,xmax]*res, [ymax,ymax]*res, /addcmd,/noerase,color = 'white',thick = 0.9,position = pos1[*,p],/overplot, linestyle = 0
         ;cgplot, [xmin,xmin]*res, [ymin,ymax]*res, /addcmd,/noerase,color = 'white',thick = 0.9,position = pos1[*,p],/overplot, linestyle = 0
         ;cgplot, [xmax,xmax]*res, [ymin,ymax]*res, /addcmd,/noerase,color = 'white',thick = 0.9,position = pos1[*,p],/overplot, linestyle = 0

         ; the CRISP box
         cgplot, [x1,x2]*res, [y1,y1]*res, /addcmd,/noerase,color = 'white',thick = 1.5,position = pos1[*,p],/overplot, linestyle = 1
         cgplot, [x1,x2]*res, [y2,y2]*res, /addcmd,/noerase,color = 'white',thick = 1.5,position = pos1[*,p],/overplot, linestyle = 1
         cgplot, [x1,x1]*res, [y1,y2]*res, /addcmd,/noerase,color = 'white',thick = 1.5,position = pos1[*,p],/overplot, linestyle = 1
         cgplot, [x2,x2]*res, [y1,y2]*res, /addcmd,/noerase,color = 'white',thick = 1.5,position = pos1[*,p],/overplot, linestyle = 1


         
  endif


  
  

endfor
  ; save the path map
cgcontrol,output = savedir+'path_'+strtrim(fn,2)+'_un_updated_RGB.eps',/ps_encap, ps_char = 0.7,/adjustsize, resize = [2*nxx,2*nyy+110]

  cgdelete
  
  stop
end
