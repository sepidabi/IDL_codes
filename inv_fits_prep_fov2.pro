; inv_fits_prep.pro

;; purpose = Prepare fits cubes to be used in
;; 'prepare_data.py' for inversions

pro inv_fits_prep_fov

cgdelete,/all

  ;TO BE DEFINED
;  w_h = 7                       ;h alpha line core position
;  w_h_cont = 0
;  w_ca = 7
;  w_ca_cont = 42
fn = 28                       ;desired frame
fbn = 5
pxn = 1 ; test pixel index along its path

datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
fibdir = 'fr'+strtrim(long(fn),2)+'/'
savedir = '/scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/'
invdir = '/home/seki2695/INV/stic/example_me/'
  
pref = ['6302','8542','3950']
pref_name = ['fe','ca8','cak']
cube =file_search(datadir+['crispex.*.6302.*_CHROMIS.fcube' , 'crispex.*.8542.*_CHROMIS.fcube' , 'crispex_3950_*_time-corrected.fcube'])
sp_cube = file_search(datadir+['crispex.*.6302.*_CHROMIS_sp.fcube' , 'crispex.*.8542.*_CHROMIS_sp.fcube' , 'crispex_3950_*_time-corrected_sp.fcube'])
calib_fits = ['calib.6302.fits','calib.8542.fits','calib.3934.fits']
cw = [6302.4935d0, 8542.091d0, 3933.6640d]
cak_int_un = lp_get(datadir+'cak_int_un.fcube',fn)

for i = 0L, n_elements(pref)-1 do begin

   ;; Image Dimensions
   lp_header, cube[i], nx = nxx, ny = nyy      ;image dimensions
   lp_header, sp_cube[i], nx = nw, ny = ntt    ;frame number and wavelength
   map = lp_read(cube[i])

   ;; Calibration info
   calib = readfits(datadir+calib_fits[i],/sile) ; [I_cal, w_cal, fx = 1.]

   ;; Ca II K Calibration
   ;; ===========
   if pref[i] eq '3950' then begin ; in case of no pol data(Cak)
      w_up = 20                    ; to take only Ca ii K line into account
      spect=f0(file_search(datadir+'wav*'+pref[i]+'*f0'))
      spect_pos = (spect-spect[10])*1.e10 ; ->also converted to angstrom
      nw_cak = n_elements(spect_pos)
      w = spect_pos(0:w_up)     ; spectral positions
      wcont = spect_pos(-1)     ; continuum point
      map_calib = dblarr(nxx,nyy,w_up+2,5) ; final calibrated cak cube (including cont. point)

      ;; No polarimetry
      map_I = map[*,*, fn*nw+0:fn*nw+w_up]
      map_cont = map[*,*, fn*nw+nw_cak - 1]
      
      for xx = 0L, nxx-1 do begin
         for yy = 0L, nyy-1 do begin
            ;; filling the spectral line
            map_calib[xx,yy,0:w_up,0] = w+calib[1]+cw[i]
            map_calib[xx,yy,0:w_up,1] = (map_I[xx,yy,*]/calib[0])*1000. ;/calib[0] ; filling stokes I
            map_calib[xx,yy,0:w_up,2:4] = 0. ; filling stokes Q, U, V as zero
            ;; filling the continuum point 4000
            map_calib[xx,yy,w_up+1,0] = wcont+calib[1]+cw[i]
            map_calib[xx,yy,w_up+1,1] = (map_cont[xx,yy]/calib[0])*1000. ;; filling stokes I -> continuum point
            map_calib[xx,yy,w_up+1,2:4] = 0. ; filling stokes Q, U, V as zero -> cont
            ;stop
         endfor
      endfor
      writefits, savedir+'map_'+pref[i]+'_fr'+strtrim(fn,2)+'.fits', map_calib
      print, 'Created Cube:' + savedir+'map_'+pref[i]+'_fr'+strtrim(fn,2)+'.fits'
      ;stop

   ;; Polarimetric Lines Calibration
   ;; =================
   endif else begin
      w_up = nw - 1
      spect_pos=f0(file_search(datadir+'wav*'+pref[i]+'*f0'))
      w = spect_pos(0:w_up)     ;spectral positions
      map_calib = dblarr(nxx,nyy,nw,5) ; final calibrated cube
      
      ;; in case of polarimetric data (Fe and Ca8)
      map_I = map[*,*,fn*nw*4+nw*0+0:fn*nw*4+nw*0+w_up]
      map_Q = map[*,*,fn*nw*4+nw*1+0:fn*nw*4+nw*1+w_up]
      map_U = map[*,*,fn*nw*4+nw*2+0:fn*nw*4+nw*2+w_up]
      map_V = map[*,*,fn*nw*4+nw*3+0:fn*nw*4+nw*3+w_up]

      for xx = 0L, nxx-1 do begin
         for yy = 0L, nyy-1 do begin
            map_calib[xx,yy,*,0] = (w+calib[1])*calib[2] +cw[i]
            map_calib[xx,yy,*,1] = map_I[xx,yy,*]/calib[0] ; filling stokes I
            map_calib[xx,yy,*,2] = map_Q[xx,yy,*]/calib[0] ; filling stokes Q
            map_calib[xx,yy,*,3] = map_U[xx,yy,*]/calib[0] ; filling stokes U
            map_calib[xx,yy,*,4] = map_V[xx,yy,*]/calib[0] ; filling stokes V
            ;stop
         endfor
      endfor
      
      writefits, savedir+'map_'+pref[i]+'_fr'+strtrim(fn,2)+'.fits', map_calib
      print, 'Created Cube:' + savedir+'map_'+pref[i]+'_fr'+strtrim(fn,2)+'.fits'

   endelse
   ;stop
endfor
print, 'To create the test individual pixel, first modify the script and then continue!! XD'
stop
  ;plotting the test inversion pixel on the FOV
  cgwindow
  c_map = lp_get(cube[2], fn*nw_cak+nw_cak-1)
  cgimage, c_map, /keep, /axis, /add, /scal
  cgimage, cak_int_un, /keep, /axis, /add, /scal
  cgplot,/add,/over,f_x_pts,f_y_pts,color = 'white'
  cgplot,/add,/over,f_x_pts[pxn],f_y_pts[pxn],color = 'red',psym = 1, symsize = 2

  cgplot,/add,/over,f_x_pts,f_y_pts,color = 'white'
  cgplot,/add,/over,f_x_pts[pxn],f_y_pts[pxn],color = 'red',psym = 1, symsize = 2

  cgplot,/add,/over,b_x_pts,b_y_pts,color = 'black'
  cgplot,/add,/over,b_x_pts[pxn],b_y_pts[pxn],color = 'blue',psym = 1, symsize = 2


  save, filename = invdir + 'inv_in'+strtrim(long(fbn),2)+'_'+strtrim(long(pxn),2)+'.sav', c_map, cak_int_un, f_x_pts, f_y_pts, b_x_pts, b_y_pts, pxn, fbn, invdir
  ;cgcontrol,output = invdir+'inv_map_'+strtrim(long(fbn),2)+'_'+strtrim(long(pxn),2)+'.eps'

stop
end
