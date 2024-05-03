; inv_fits_prep.pro

;; purpose = Prepare fits cubes to be used in
;; 'prepare_data.py' for inversions

pro inv_fits_prep

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
;art_shift = [0.5,0.,-0.3] ; defined to adopt with example shifts
cube =file_search(datadir+['crispex.*.6302.*_CHROMIS.fcube' , 'crispex.*.8542.*_CHROMIS.fcube' , 'crispex_3950_*_time-corrected.fcube'])
sp_cube = file_search(datadir+['crispex.*.6302.*_CHROMIS_sp.fcube' , 'crispex.*.8542.*_CHROMIS_sp.fcube' , 'crispex_3950_*_time-corrected_sp.fcube'])
calib_fits = ['calib.6302.fits','calib.8542.fits','calib.3934.fits']
cw = [6302.4935d0, 8542.091d0, 3933.6640d]
cak_int_un = lp_get(datadir+'cak_int_un.fcube',fn)

for i = 0L, n_elements(pref)-1 do begin

  ;Image Dimensions
  lp_header, cube[i], nx = nxx, ny = nyy ;image dimensions
  lp_header, sp_cube[i], nx = nw, ny = ntt ;frame number and wavelength

  ;;fibrils file
  file = file_search(datadir+fibdir+'*'+pref[i]+'*.csav')
  path_n = n_elements(file)  ; number of all the defined paths
  fib_n = path_n/2           ; number of fibrils
  calib = readfits(datadir+calib_fits[i],/sile) ; [I_cal, w_cal, fx = 1.]

  fbn = path_n/2-1
  for n = 0L, fbn do begin     

     ;fibril
     restore, file(2*n)
     F = loop_slab  ; intensity values in the desired frame
     F_x = x_coords
     F_x_pts = x_loop_pts
     F_y = y_coords
     F_l = loop_size
     F_y_pts = y_loop_pts
     
     ;fibril background
     restore, file(2*n+1)
     B = loop_slab ; intensity values in the desired frame
     B_x = x_coords
     B_x_pts = x_loop_pts
     B_y_pts = y_loop_pts
     B_y = y_coords
     B_l = loop_size

     if pref[i] eq '3950' then begin       ; in case of no pol data(Cak)
        w_up = 20 ; to take only Ca ii K line into account
        ;restore,datadir+'spectfile.'+pref[i]+'.idlsave'
        spect=f0(file_search(datadir+'wav*'+pref[i]+'*f0'))
        spect_pos = (spect-spect[10])*1.e10   ; ->also converted to angstrom
        nw_cak = n_elements(spect_pos)
        w = spect_pos(0:w_up)   ;spectral positions
        wcont = spect_pos(-1)   ; continuum point
        obs_f = dblarr(F_l, w_up+2, 5)
        obs_B = dblarr(B_l, w_up+2, 5)
        for ff = 0L, F_l-1 do begin
           ; filling the spectral line
           obs_f[ff,0:w_up,0] = w+cw[i]+calib[1] ;(w+calib[1])*calib[2]+cw[i]
           obs_f[ff,0:w_up,1] = (F[ff, fn, 0:w_up])*1000. ;/calib[0] ; filling stokes I
           obs_f[ff,0:w_up,2:4] = 0.                       ; filling stokes Q, U, V as zero
           ;filling the continuum point 4000
           obs_f[ff,w_up+1,0] = wcont+cw[i]+calib[1] ;(wcont+calib[1])*calib[2]+cw[i]
           obs_f[ff,w_up+1,1] = (F[ff, fn, -1])*1000. ;/calib[0] ; filling stokes I -> continuum point
           obs_f[ff,w_up+1,2:4] = 0.                       ; filling stokes Q, U, V as zero -> cont
         ;stop
        endfor
        for bb = 0L, B_l-1 do begin
           obs_b[bb,0:w_up,0] = w+cw[i]+calib[1] ; (w+calib[1])*calib[2]+cw[i]
           obs_b[bb,0:w_up,1] = B[bb, fn, 0:w_up]*1000. ;/calib[0] ; filling stokes I
           obs_b[bb,0:w_up,2:4] = 0.                       ; filling stokes Q, U, V as zero
           ;filling the continuum point 4000
           obs_b[bb,w_up+1,0] = wcont+cw[i]+calib[1] ; (wcont+calib[1])*calib[2]+cw[i]
           obs_b[bb,w_up+1,1] = B[bb, fn, -1]*1000. ;/calib[0] ; filling stokes I -> continuum point
           obs_b[bb,w_up+1,2:4] = 0.                       ; filling stokes Q, U, V as zero -> cont
        ;stop
        endfor
        writefits, savedir+'f_obs'+pref[i]+'_'+strmid(file[2*n],112,9)+'.fits', obs_f
        writefits, savedir+'b_obs'+pref[i]+'_'+strmid(file[2*n],112,9)+'.fits', obs_b
        print, 'Created Cube:' + savedir+'..._obs'+pref[i]+'_'+strmid(file[2*n],112,9)+'.fits'
     ;stop
     endif else begin           ; in case of polarimetric data (Fe and Ca8)
        w_up = nw - 1
        ;restore,datadir+'spectfile.'+pref[i]+'.idlsave'
        spect_pos=f0(file_search(datadir+'wav*'+pref[i]+'*f0'))
        w = spect_pos(0:w_up)   ;spectral positions
        obs_f = dblarr(F_l, w_up+1, 5)
        obs_B = dblarr(B_l, w_up+1, 5)

        for ff = 0L, F_l-1 do begin
           obs_f[ff,*,0] = (w+calib[1])*calib[2] +cw[i]
           obs_f[ff,*,1] = F[ff, fn, *,0]/calib[0] ; filling stokes I
           ;obs_f[ff,*,2] = F[ff, fn, *,1]/calib[0] ; filling stokes Q
           ;obs_f[ff,*,3] = F[ff, fn, *,2]/calib[0] ; filling stokes U
           ;obs_f[ff,*,4] = F[ff, fn, *,3]/calib[0] ; filling stokes V
        ;stop
        endfor
        for bb = 0L, B_l-1 do begin
           obs_b[bb,*,0] = (w+calib[1])*calib[2] +cw[i]
           obs_b[bb,*,1] = B[bb, fn, *,0]/calib[0] ; filling stokes I
           ;obs_b[bb,*,2] = B[bb, fn, *,1]/calib[0] ; filling stokes Q
           ;obs_b[bb,*,3] = B[bb, fn, *,2]/calib[0] ; filling stokes U
           ;obs_b[bb,*,4] = B[bb, fn, *,3]/calib[0] ; filling stokes V
        ;stop
        endfor
        writefits, savedir+'f_obs'+pref[i]+'_'+strmid(file[2*n],104,9)+'.fits', obs_f
        writefits, savedir+'b_obs'+pref[i]+'_'+strmid(file[2*n],104,9)+'.fits', obs_b
        print, 'Created Cube:' + savedir+'..._obs'+pref[i]+'_'+strmid(file[2*n],104,9)+'.fits'

     endelse
     ;stop
  endfor

  ;;create 1-px test files
  ;;fibril:
  if pref[i] eq '3950' then begin
     test_f = dblarr(w_up+2,5)
     test_f[*,*] = obs_f[pxn,*,*]
     ;;BG:
     test_b = dblarr(w_up+2,5)
     test_b[*,*] = obs_b[pxn,*,*]
  endif else begin
     ;;fibril:
     test_f = dblarr(w_up+1,5)
     test_f[*,*] = obs_f[pxn,*,*]
     ;;BG:
     test_b = dblarr(w_up+1,5)
     test_b[*,*] = obs_b[pxn,*,*]
  endelse

  ;Save the 1-pixel test to .fits
  writefits, invdir+'testf_'+pref_name[i]+'obs'+strtrim(long(fbn),2)+'_'+strtrim(long(pxn),2)+'.fits', test_f  
  writefits, invdir+'testb_'+pref_name[i]+'obs'+strtrim(long(fbn),2)+'_'+strtrim(long(pxn),2)+'.fits', test_b

  ;writefits, invdir+'testf_'+pref_name[i]+'obs.fits', test_f  
  ;writefits, invdir+'testb_'+pref_name[i]+'obs.fits', test_b
  
  ;stop
endfor
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
