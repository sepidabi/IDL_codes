pro fibril_map_backup

  cgdelete,/all
  ;TO BE DEFINED
  w_h = 7                       ;h alpha line core position
  w_ca = 7
  fn = 01                       ;desired frame
  sm_fact = 1.

  ;crop crisp edges with [x1,x2,y1,y2]
  crop = [245,1678,10,1130]

  cak_cube = 'crispex_3950_2016-09-15T08:49:53_scans=24-57_time-corrected.fcube'
  caksp_cube = 'crispex_3950_2016-09-15T08:49:53_scans=24-57_time-corrected_sp.fcube'
  cakint_cube = 'cak_int.fcube'
  cakintun_cube = 'cak_int_un.fcube'
  h_cube = 'crispex.6563.08:49:51.time_corrected_CHROMIS.fcube'
  hsp_cube = 'crispex.6563.08:49:51.time_corrected_CHROMIS_sp.fcube'
  hint_cube = 'h_int.fcube'
  hintun_cube = 'h_int_un.fcube'

  fe_cube = 'crispex.stokes.6302.08:49:22.time_corrected_CHROMIS.fcube'
  fesp_cube = 'crispex.stokes.6302.08:49:22.time_corrected_CHROMIS_sp.fcube'
  
  ca8_cube = 'crispex.stokes.8542.08:49:51.time_corrected_CHROMIS.fcube'
  ca8sp_cube = 'crispex.stokes.8542.08:49:51.time_corrected_CHROMIS_sp.fcube'
  
  datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
  savedir = '/scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/'
  fibdir = 'fr'+strtrim(long(fn),2)+'/'
  file_ca = file_search(datadir+fibdir+'*3950*.csav')
  file_h = file_search(datadir+fibdir+'*6563*.csav')
  res = 0.0375

  ;Image Dimensions
  lp_header, datadir+cak_cube, nx = nxx, ny = nyy ;image dimensions
  lp_header, datadir+caksp_cube, nx = nw_ca, ny = ntt ;frame number and Ca K wavelength
  lp_header,  datadir+hsp_cube, nx = nw_h             ;H alpha wavelength
  lp_header,  datadir+fesp_cube, nx = nw_fe           ;Fe wavelength
  lp_header,  datadir+ca8sp_cube, nx = nw_ca8
  nss = 4

  ;==========================================================
                          ;Stokes maps
  ;==========================================================

  ;continuum points
  Ic = lp_get(datadir+fe_cube, fn*nw_fe*nss + 0*nw_fe + 9)
  Qc = lp_get(datadir+fe_cube, fn*nw_fe*nss + 1*nw_fe + 9)
  Uc = lp_get(datadir+fe_cube, fn*nw_fe*nss + 2*nw_fe + 9)
  Vc = lp_get(datadir+fe_cube, fn*nw_fe*nss + 3*nw_fe + 9)
  ;==========================================================
                          ;Stokes maps
  ;==========================================================

  ;Initial Value of LP and CP
   LPt_l = fltarr(nxx,nyy)
   LPt_r = fltarr(nxx,nyy)
   CPt_l = fltarr(nxx,nyy)
   CPt_r = fltarr(nxx,nyy)
   StokesI = fltarr(nxx,nyy)
  ;==========================================================
                          ;Stokes maps
  ;==========================================================

   ;Total CP over the left-side of "selected" spectral positions
   for l = 0L, 3 do begin
      del = 1.
      StV_l = lp_get(datadir+fe_cube, fn*nw_fe*nss + 3*nw_fe + l)
      CPt_l = CPt_l + del*smooth(reform(StV_l), [sm_fact,sm_fact], /edge_truncate)

      StQ = lp_get(datadir+fe_cube, fn*nw_fe*nss + 1*nw_fe + l)
      StU = lp_get(datadir+fe_cube, fn*nw_fe*nss + 2*nw_fe + l)
      LPt_l = LPt_l + sqrt(smooth(reform(StQ), [sm_fact,sm_fact], /edge_truncate)^2 + smooth(reform(StU), [sm_fact,sm_fact], /edge_truncate)^2)

   endfor
  ;==========================================================
                          ;Stokes maps
  ;==========================================================

   ;Total CP over the right-side of "selected" spectral positions
   for r = 0L, 3 do begin
      del = -1.
      StV_r = lp_get(datadir+fe_cube, fn*nw_fe*nss + 3*nw_fe + r+5)
      CPt_r = CPt_r + del*smooth(reform(StV_r), [sm_fact,sm_fact], /edge_truncate)

      StQ = lp_get(datadir+fe_cube, fn*nw_fe*nss + 1*nw_fe + l)
      StU = lp_get(datadir+fe_cube, fn*nw_fe*nss + 2*nw_fe + l)
      LPt_r = LPt_r + sqrt(smooth(reform(StQ), [sm_fact,sm_fact], /edge_truncate)^2 + smooth(reform(StU), [sm_fact,sm_fact], /edge_truncate)^2)

   endfor
  ;==========================================================
                          ;Stokes maps
  ;==========================================================

   ;Averaging over the 8 wavelengths around the 1st line core (630.1 nm)
   CP = float(CPt_l+CPt_r)/float(8.*Ic)
   LP = float(LPt_l+LPt_r)/float(8.*Ic)
   CPcont = float(Vc)/float(Ic)
   LPcont = float(sqrt(Qc^2+Uc^2))/float(Ic)
   ;sigmaCP = stddev(reform(CPcont[x1:x2,y1:y2,i]))/sqrt(8.)
   ;sigmaLP = stddev(reform(LPcont[x1:x2,y1:y2,i]))/sqrt(8.)

;=============================================================================
;=============================================================================
  ;stop
  ;maps

  def = fltarr(nxx,nyy)
  ;int, un Ca II K
  img_ca = cgGmaScl(lp_get(datadir+cakintun_cube,fn),min = 1.034e-9,max = 6.043e-8,gamma = 0.8)
  ;int, un H alpha
  img_h = cgGmaScl(lp_get(datadir+hintun_cube,fn),min = 7040.089,max = 26250.090,gamma = 0.8)
  ;...cropping edges
  def(crop[0]:crop[1],crop[2]:crop[3]) = img_h(crop[0]:crop[1],crop[2]:crop[3])
  img_h = def
  
  ;cont 400 nm
  img_cont = cgGmaScl(lp_get(datadir+cak_cube,fn*nw_ca+nw_ca-1),min=1.194e-8,max=3.620e-8,gamma=0.65)

  ;CP
  img_cp = cgGmaScl(CP, gamma = 2.)
  ;...cropping edges
  def(crop[0]:crop[1],crop[2]:crop[3]) = img_cp(crop[0]:crop[1],crop[2]:crop[3])
  img_cp = def

  ;Ca 8542 core
  img_ca8 = cgGmaScl(lp_get(datadir+ca8_cube,fn*nw_ca8*nss+nw_ca8*0+10),gamma = 0.85, min = 47.607, max = 187.330)
  def(crop[0]:crop[1],crop[2]:crop[3]) = img_ca8(crop[0]:crop[1],crop[2]:crop[3])
  img_ca8 = def

  
  ;LP
  img_lp = cggmascl(LP,min = 0.0003,max = 0.01,gamma = 1.1)
  ;...cropping edges
  def(crop[0]:crop[1],crop[2]:crop[3]) = img_lp(crop[0]:crop[1],crop[2]:crop[3])
  img_lp = def

  ;Ca II K line-core
  img_ca_sing = cgGmaScl(lp_get(datadir+cak_cube,fn*nw_ca+11),min = 7.343e-10,max = 8.000e-9,gamma = 0.6)

  ;H alpha line-core
  img_ha_sing = cgGmaScl(lp_get(datadir+h_cube,fn*nw_h+8),min = 1050.66,max = 2903.823,gamma = 0.8)
  ;...cropping edges
  def(crop[0]:crop[1],crop[2]:crop[3]) = img_ha_sing(crop[0]:crop[1],crop[2]:crop[3])
  img_ha_sing = def

  ;PLOTS...
  ;===========
  cgwindow ;, wxsize =2*nxx,wysize = 2*nyy;,wmulti = [0,2,1] ;,woxmargin = [1,0],woymargin = [4,0.5]
  pos1 = cgLayout([2,4], OXMargin=[6.5,1], OYMargin=[2,1], XGap=3, Ygap=0)
  
  ;maps
  map = [img_ca,img_h,img_cp,img_lp,img_ca8,img_cont,img_ca_sing,img_ha_sing]

  ;captions
  caption = ['g) Ca II K $\lambda$-integrated & unsharp-masked','h) H$\alpha$ $\lambda$-integrated & unsharp-masked','e) CP at photosphere', 'f) LP at photosphere','c) Ca II 8542 $\angstrom$ core','d) continuum 4000 $\angstrom$','a) Ca II K core','b) H$\alpha$ core']

  ;caption = [' Ca II K $\lambda$-integrated & unsharp-masked',' H$\alpha$ $\lambda$-integrated & unsharp-masked',' Ca II 8542 $\angstrom$ core', ' LP at photosphere',' continuum 4000 $\angstrom$',' CP at photosphere',' Ca II K core',' H$\alpha$ core']


  ;caption color
  ccolor = [replicate('white',5),'black',replicate('white',2)]
  
  ;axis formats
  ax0 = {xticklen:-0.02,yticklen:-0.02,xtitle:'x [arcsec]',xtickformat:"(I0)",ytickformat:"(I0)",ytitle:'y [arcsec]'}
  ax1 = {xticklen:-0.02,yticklen:-0.02,xtitle:'x [arcsec]',xtickformat:"(I0)",ytickformat:"(A1)",ytitle: ''}
  ax2 = {xticklen:-0.02,yticklen:-0.02,xtitle:'',xtickformat:"(A1)",ytickformat:"(I0)",ytitle: 'y [arcsec]'}
  ax3 = {xticklen:-0.02,yticklen:-0.02,xtitle:'',xtickformat:"(A1)",ytickformat:"(A1)",ytitle: ''}
  axis = [ax0,ax1,ax2,ax3,ax2,ax3,ax2,ax3]
  
  ;positions
  pos = [[pos1(*,6)],[pos1(*,7)],[pos1(*,4)],[pos1(*,5)],[pos1(*,2)],[pos1(*,3)],[pos1(*,0)],[pos1(*,1)]]

;PLOTTING THE BACKGROUND IMAGE
;=============================
for p = 0L, n_elements(map(*,0))/nxx -1 do begin
  ;plotting background maps
   cgimage,map(p*nxx:(p+1)*nxx-1,*),/scale,/keep,/axes,/addcmd,axkeywords = axis(p),xrange = [0,nxx*res],yrange = [0,nyy*res],position = pos[*,p]

   ;caption
   cgtext,pos(0,p)+1, pos(1,p)+1, caption(p), /add, color = ccolor[p]
;stop   
   if p lt 2 then begin
      
  ;overplotting fibrils
  ;=====================
      for n = 0L, n_elements(file_ca)/2-1 do begin

         ;fibrils
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

  
         ;fibril backgrounds
         ;Ca K
         restore, file_ca(2*n+1)
         B_ca = loop_slab
         B_ca_x = x_coords
         B_ca_y = y_coords
         B_ca_l = loop_size

         ;overplot fibril
         cgplot, F_ca_x*res,F_ca_y*res,/addcmd,/noerase,color = 'red',thick = 2.,position = pos1[*,0],/overplot,aspect=2./3 ;, linestyle = 2
         ;overplot fibril background
         cgplot, B_ca_x*res,B_ca_y*res,/addcmd,/noerase,color = 'dodger blue',thick = 0.9,position = pos1[*,0],/overplot

      endfor
      
   endif

endfor

  ;save the maps
cgcontrol,output = savedir+'path_'+strtrim(fn,2)+'_maps.eps',/ps_encap, ps_char = 0.8,/adjustsize, resize = [2*nxx-100,4*nyy]
 ;cgcontrol,output = savedir+'path_'+strtrim(fn,2)+'_maps_present.eps',/ps_encap, ps_char = 0.8,/adjustsize, resize = [2*nxx-100,4*nyy]


  ;stop
  cgdelete
  
  stop
end
