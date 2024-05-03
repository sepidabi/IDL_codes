
cgdelete, /all

  dir='/scratch/sepid/DATA/AR/plage/2016.09.15/'
  savedir = '/scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/'
  fn = 28                       ; frame number
  ;;Gregal's satlas
  satlas = '~gviss/bin/idl/collibs/tlibb/fts_disk_center_SI.idlsave'

;; This should be the coords of
;; a quiet region IN THE DATA
;; to fit the solar atlas to:
x0=911
y0=68
x1=1262
y1=260

if(0) then begin
   ;; calibrated cube
   cube = readfits('/scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/map_8542_fr28.fits')
   w = cube[0,0,*,0]
   cw = 8542.091d0              ; central wavelength
   w = w - cw      ; calibrated data wavelength positions
   ;; actual cubes
   fil0 = dir+'crispex.stokes.8542.*.time_corrected_CHROMIS.fcube'
   fil1 = dir+'crispex.stokes.8542.*.time_corrected_CHROMIS_sp.fcube'
   mapdatastokes,fil0,fil1,d,nx,ny,nt,nw
   pref = '8542'
   ptot = 1.
endif

if(1) then begin
   ;; calibrated cube
   cube = readfits('/scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/map_6302_fr28.fits')
   w = cube[0,0,*,0]
   cw = 6302.4935d0
   w = w - cw      ; calibrated data wavelength positions
   ;; actual cubes
   fil0 = dir+'crispex.stokes.6302.*.time_corrected_CHROMIS.fcube'
   fil1 = dir+'crispex.stokes.6302.*.time_corrected_CHROMIS_sp.fcube'
   mapdatastokes,fil0,fil1,d,nx,ny,nt,nw
   pref = '6302'
   ptot = 1.
endif


;;calculates the average value
;;over the defined QS area
;;in each wavelength
imean = dblarr(nw)
for ii=0,nw-1 do imean[ii] = mean(cube[x0:x1, y0:y1,ii,1])


;; Extract the range of wavelength and corresponding intensity
;; and put it in x , y respectively "FROM SOLAR ATLAS" + to cgs
red_satlas, w[0]+cw-1.0,  w[-1]+cw+1.0, x, y, /cgs
x-= cw          ; substituting the cw to get the difference from the line center position
;stop
fpi = cfpi(cw) ; object related to instrument profile
dw = x[1]-x[0] ; wave length range
np = 51 ; ???
tw = (dindgen(np)-np/2)*dw ; adjusting the number positions from satlas to 51
tr = fpi->dual_fpi(tw+cw, erh=-0.015) ;;it returns a structure with the nominal reflectivities
                                                             ;;and cavity separations of the CRISP spectropolarimeter
                                                             ;;USING THE NEW POSITIONS FROM TW
y1 = fftconvol(y*ptot, tr/total(tr)) ; convolves the satlas intensity to a broader profile

cgwindow
cgplot, x+cw, y1, linestyle = 1,/add,xrange = [min(x+cw),max(x+cw)]
cgplot, w+cw, imean,/add,/over
cglegend,titles = ['S-Atlas (cgs units)','Fit'],linestyles = ['1','0'],/add

;cgcontrol,create_pdf = savedir+'test.calib.'+pref+'.pdf'
;writefits, dir+'calib.'+pref+'.fits', par

stop
end
