; inv_fits_prep2.pro

;; purpose = Prepare fits cubes to be used in
;; 'prepare_data.py' for inversions of Project II

pro inv_fits_prep2

cgdelete,/all

;; TO BE DEFINED

fr = 29                         ;desired frame
fbn = 5
pxn = 1 ; test pixel index along its path

datadir = '/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/'
fibdir = 'fr'+strtrim(long(fr),2)+'/'
calibdir = '/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/w_calibrated/'
invdir = '/home/seki2695/INV/stic/example_me/'
  
pref = ['6173','8542','3950']
pref_name = ['fe','ca8','cak']

cube =file_search(datadir+['crispex.*.6173.*_CHROMIS.fcube' ,$
                           'crispex.*.8542.*_CHROMIS.fcube' ,$
                           'crispex_3950_*_time-corrected.fcube'])
   
sp_cube = file_search(datadir+['crispex.*.6173.*_CHROMIS_sp.fcube' ,$
                               'crispex.*.8542.*_CHROMIS_sp.fcube' ,$
                               'crispex_3950_*_time-corrected_sp.fcube'])

calib_fits = [(file_search(calibdir+'*6173*fits'))[fr],$
              (file_search(calibdir+'*8542*fits'))[fr],$
              (file_search(calibdir+'*3934*fits'))[fr]]

cw = [6173.3355d0, 8542.091d0, 3933.6640d]

for i = 0L, n_elements(pref)-1 do begin

  ;Image Dimensions
  lp_header, cube[i], nx = nxx, ny = nyy ;image dimensions
  lp_header, sp_cube[i], nx = nw, ny = ntt ;frame number and wavelength

  ;; slabs file
  file = file_search(datadir+fibdir+'*'+pref[i]+'*.csav')
  calib = readfits(calib_fits[i],/sile) ; [I_cal, w_cal, fx = 1.]
    
  for n = 0L, n_elements(file)-1 do begin     

     ;; CRISPEX slab
     restore, file[n]
     
     if pref[i] eq '3950' then begin       ; in case of no pol data(Cak)
        w_up = 26               ; to take only Ca ii K line into account
        w_c = 13
        spect=f0(file_search(datadir+'wav*'+pref[i]+'*f0'))
        spect_pos = (spect-spect[w_c])*1.e10 ; ->also converted to angstrom
        nw_cak = n_elements(spect_pos)
        w = spect_pos(0:w_up)   ;spectral positions
        wcont = spect_pos(-1)   ; continuum point
        slab = dblarr(loop_size, w_up+2, 5)
        
        for ff = 0L, loop_size-1 do begin
           ;; filling the spectral line
           slab[ff,0:w_up,0] = w+cw[i]+calib[1] ;(w+calib[1])*calib[2]+cw[i]
           slab[ff,0:w_up,1] = (loop_slab[ff, fr, 0:w_up])*1000. ;/calib[0] ; filling stokes I
           slab[ff,0:w_up,2:4] = 0. ; filling stokes Q, U, V as zero
           
           ;; filling the continuum point 4000
           slab[ff,w_up+1,0] = wcont+cw[i]+calib[1] ;(wcont+calib[1])*calib[2]+cw[i]
           slab[ff,w_up+1,1] = (loop_slab[ff, fr, -1])*1000. ;/calib[0] ; filling stokes I -> continuum point
           slab[ff,w_up+1,2:4] = 0. ; filling stokes Q, U, V as zero -> cont

        endfor
        filename = calibdir + 'slab_obs'+pref[i]+'_'+strmid(file[n],121,9)+'.fits'
        filename_test = calibdir + 'test_obs'+pref[i]+'_'+strmid(file[n],121,9)+'.fits'

     ;; in case of polarimetric data (Fe and Ca8)
     endif else begin
        w_up = nw - 1
        spect_pos=f0(file_search(datadir+'wav*'+pref[i]+'*f0'))
        w = spect_pos(0:w_up)   ;spectral positions
        slab = dblarr(loop_size, w_up+1, 5)

        for px = 0L, loop_size-1 do begin
           slab[px,*,0] = (w+calib[1])*calib[2] +cw[i]
           slab[px,*,1] = loop_slab[px, fr, *,0]/calib[0] ; filling stokes I
           slab[px,*,2] = loop_slab[px, fr, *,1]/calib[0] ; filling stokes Q
           slab[px,*,3] = loop_slab[px, fr, *,2]/calib[0] ; filling stokes U
           slab[px,*,4] = loop_slab[px, fr, *,3]/calib[0] ; filling stokes V
        ;stop
        endfor

        filename = calibdir + 'slab_obs'+pref[i]+'_'+strmid(file[n],113,9)+'.fits'
        filename_test = calibdir + 'test_obs'+pref[i]+'_'+strmid(file[n],113,9)+'.fits'

     endelse
     ;stop

     writefits, filename, slab
     print, 'Created Cube: ' + filename
     writefits, filename_test, slab[0, *,*]
     print, 'Created test pixel: '+ filename_test
     ;stop
  endfor

;stop
endfor

stop
;; create 1-px test files
  ;; fibril:
  if pref[i] eq '3950' then begin
     test_f = dblarr(w_up+2,5)
     test_f[*,*] = slab[pxn,*,*]
  endif else begin
     ;;fibril:
     test_f = dblarr(w_up+1,5)
     test_f[*,*] = slab[pxn,*,*]
 endelse

  ;; Save the 1-pixel test to .fits
  writefits, invdir+'testf_'+pref_name[i]+'obs'+strtrim(long(fbn),2)+'_'+strtrim(long(pxn),2)+'.fits', test_f  
  writefits, invdir+'testb_'+pref_name[i]+'obs'+strtrim(long(fbn),2)+'_'+strtrim(long(pxn),2)+'.fits', test_b

  ;writefits, invdir+'testf_'+pref_name[i]+'obs.fits', test_f  
  ;writefits, invdir+'testb_'+pref_name[i]+'obs.fits', test_b
  
  ;stop

  ;plotting the test inversion pixel on the FOV
  cgwindow
  c_map = lp_get(cube[2], fr*nw_cak+nw_cak-1)
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

