pro fibril_map

  cgdelete,/all
  ;TO BE DEFINED
  w_h = 7                       ;h alpha line core position
  w_ca = 7
  fn = 01                       ;desired frame

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
  lp_header, datadir+cak_cube, nx = nxx, ny = nyy ;image dimensions
  lp_header, datadir+caksp_cube, nx = nw_ca, ny = ntt ;frame number and Ca K wavelength
  lp_header,  datadir+hsp_cube, nx = nw_h ;H alpha wavelength

  img_ca = lp_get(datadir+cakintun_cube,fn)
  img_h = lp_get(datadir+hintun_cube,fn)
  ;img_ca_bri = lp_get(datadir+cak_cube,fn*nw_ca+0)
  ;img_ca_sing = lp_get(datadir+cak_cube,fn*nw_ca+8)
  img_cont = lp_get(datadir+cak_cube,fn*nw_ca+nw_ca-1)
  ;img_ca_int = lp_get(datadir+cakint_cube,fn)

  ;test = img_ca_int - 2*img_ca_bri
  ;less = where(test lt min(img_ca_int))
  ;test(less) = min(img_ca_int)

  ;stop

  ;display fibril paths on the data
  cgwindow, wxsize =2*nxx,wysize = nyy;,wmulti = [0,2,1] ;,woxmargin = [1,0],woymargin = [4,0.5]
  pos1 = cgLayout([2,1], OXMargin=[2,0], OYMargin=[8,2], XGap=0)
  ;maps
  map = [img_ca,img_h]
  ;axis formats
  ax0 = {xticklen:-0.02,yticklen:-0.02,xtitle:'[arcsec]',ytickformat:"(I0)",ytitle:'[arcsec]'}
  ax1 = {xticklen:-0.02,yticklen:-0.02,xtitle:'[arcsec]',ytickformat:"(A1)",ytitle: ''}
  axis = [ax0,ax1]

for p = 0L, n_elements(map(*,0))/nxx -1 do begin
  ;plotting background maps
  cgimage,map(p*nxx:(p+1)*nxx-1,*),/scale,/keep_,/axes,/addcmd,axkeywords = axis(p),xrange = [0,nxx*res],yrange = [0,nyy*res],position = pos1[*,p]


endfor

stop

  ;save the path map
  cgcontrol,output = savedir+'path_'+strtrim(fn,2)+'_un_unpdated.eps',/ps_encap, ps_char = 0.7,/adjustsize
  cgdelete
  
  stop
end
