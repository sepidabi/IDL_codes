;pro calib_data

  function fitme, p, x=x, y=y, imean=imean, w=w, sig=sig

  y2 = interpol(y, x-p[1], w*p[2])*p[0]
  plot, w, imean, /line
  oplot, w, y2
  ;stop
  return, (y2-imean)*sig
end


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

;x0=561
;y0=744
;x1=872
;y1=1087


if(0) then begin

   w=f0(dir+'wav.8542.f0') ; gives the spectral steps

fil0 = dir+'crispex.stokes.8542.*.time_corrected_CHROMIS.fcube'
fil1 = dir+'crispex.stokes.8542.*.time_corrected_CHROMIS_sp.fcube'
mapdatastokes,fil0,fil1,d,nx,ny,nt,nw
   pref = '8542'
   cw = 8542.091d0 ; central wavelength
   ;;defines the weight for the
   ;;wavelength steps for the fitting process
   sig = dblarr(nw)+1.0
   sig[4:12] = 0.4
   ;ptot=0.72           ; DON'T CHANGE IT defined only to equilize the Icont for all lines
                                ; according to cont 4000
   ptot = 1.
 endif

if(1) then begin

   ;restore,'slit_6302.csav'
   w=f0(dir+'wav.6302.f0')
  fil0 = dir+'crispex.stokes.6302.*.time_corrected_CHROMIS.fcube'
  fil1 = dir+'crispex.stokes.6302.*.time_corrected_CHROMIS_sp.fcube'
  mapdatastokes,fil0,fil1,d,nx,ny,nt,nw
   pref = '6302'
   cw = 6302.4935d0
    sig = dblarr(nw)+1.0
    sig[4:5] = 0.3
    sig[10:14] = 0.3
    sig[9] = 5
    ;ptot=0.64           ; DON'T CHANGE IT defined only to equilize the Icont for all lines
                                ; according to cont 4000
    ptot = 1.
 endif


;; Loads the desired time step
cube=d[fn]

;;calculates the average value
;;over the defined QS area
;;in each wavelength
imean = dblarr(nw)
for ii=0,nw-1 do imean[ii] = mean(cube[x0:x1, y0:y1,ii,0])

;------------------------------------------
;;Gregal's test
;restore,satlas ; restoring the solar atlas
;ftsint = ftsint_si/(1e-7 / (1e-2)^2)     ; get to cgs
;ftsint_sel = fltarr(nw) ; intensity corrosponding to the selected wl range
;ftswav_sel = fltarr(nw) ; selected wl range
;for ii=0,nw-1 do begin
;   dum = min(abs(ftswav-(w[ii]+cw)), wherewav)
;   ftsint_sel[ii] = ftsint[wherewav]
;   ftswav_sel[ii] = ftswav[wherewav]
;endfor
;offset_factor = (imean / ftsint_sel) ; intensity offset factor == calib[0]
;;------------------------------------------

;; Extract the range of wavelength and corresponding intensity
;; and put it in x , y respectively "FROM SOLAR ATLAS" + to cgs
red_satlas, w[0]+cw-1.0,  w[-1]+cw+1.0, x, y, /cgs
x-= cw ; substituting the cw to get the difference from the line center position
;stop
fpi = cfpi(cw) ; object related to instrument profile
dw = x[1]-x[0] ; wave length range
np = 51 ; ???
tw = (dindgen(np)-np/2)*dw ; adjusting the number positions from satlas to 51
tr = fpi->dual_fpi(tw+cw, erh=-0.015) ;;it returns a structure with the nominal reflectivities
                                                             ;;and cavity separations of the CRISP spectropolarimeter
                                                             ;;USING THE NEW POSITIONS FROM TW
y1 = fftconvol(y*ptot, tr/total(tr)) ; convolves the satlas intensity to a broader profile
;y1_test = convol()
;stop

fpar = replicate({limits:[0d0,0d0], limited:[1,1], fixed:0}, 3)
fpar[0].limits[*] = [0d0, mean(imean)/mean(y1)*10d0]
fpar[1].limits[*] = [-1d0,1.0d0]
fpar[2].limits[*] = [0.98d0,1.02d0]
fpar[2].fixed=1
farg = {x:x, y:y1, w:w, sig:sig, imean:imean}

par = [double(mean(imean)/mean(y1)), 0.001d0, 1.0d0]
print, par
par = mpfit('fitme', par, parinfo=fpar, fun=farg)

;------------------------------------
;;from Gregal's
;par[0] = mean(offset_factor)
;stop
;------------------------------------

cgwindow
;plot, x, y1, /line, xrange=[-2,2], /xstyle, ysty=3
;oplot, (w+par[1])*par[2], imean/par[0]
cgplot, x+cw, y1, linestyle = 1,/add,xrange = [min(x+cw),max(x+cw)]
cgplot, (w+par[1]+cw)*par[2], imean/par[0],/add,/over, psym = 16, symsize = 0.75
cglegend,titles = ['S-Atlas (cgs units)','Fit'],linestyles = ['1','0'],/add

cgcontrol,create_pdf = savedir+'test.calib.'+pref+'.pdf'
writefits, dir+'calib.'+pref+'.fits', par

stop
end
