;pro calib_data_2

  function fitme, p, x=x, y=y, imean=imean, w=w, sig=sig

  y2 = interpol(y, x-p[1], w*p[2])*p[0]
  ;plot, w, imean, /line
  ;oplot, w, y2
  ;stop
  return, (y2-imean)*sig
end
cgdelete,/all


dir='/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/'
savedir = '/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/w_calibrated/'
;fr = 212                         ; frame number
;;Gregal's satlas
satlas = '~gviss/bin/idl/collibs/tlibb/fts_disk_center_SI.idlsave'

;; This should be the coords of
;; a quiet region IN THE DATA
;; to fit the solar atlas to:
;; qs region 1 (2020-06-25):
x0=800 ;923
y0=226 ;575
x1=1119 ;1138
y1=549 ;869


if(0) then begin

   ;restore,'slit_6173.csav'
   w=f0(dir+'wav.6173.f0')
  fil0 = dir+'crispex.stokes.6173.*.time_corrected_CHROMIS.fcube'
  fil1 = dir+'crispex.stokes.6173.*.time_corrected_CHROMIS_sp.fcube'
  mapdatastokes,fil0,fil1,d,nx,ny,nt,nw
   pref = '6173'
   cw = 6173.3355d0
    sig = dblarr(nw)+1.0
;    sig[4:5] = 0.3
    sig[4:5] = 0.3
    sig[13] = 5
    ;sig[12]=5
    ;sig [0]=5
    ;sig[0:1] = 2
    ;sig[12:13] = 10
    ;ptot=0.64           ; DON'T CHANGE IT defined only to equilize the Icont for all lines
                                ; according to cont 4000
    ptot = 1.
 endif


if(1) then begin

   w=f0(dir+'wav.8542.f0') ; gives the spectral steps

   ;stop
   
   fil0 = dir+'crispex.stokes.8542.*.time_corrected_CHROMIS.fcube'
   fil1 = dir+'crispex.stokes.8542.*.time_corrected_CHROMIS_sp.fcube'
   mapdatastokes,fil0,fil1,d,nx,ny,nt,nw
   pref = '8542'
   cw = 8542.091d0 ; central wavelength
   ;;defines the weight for the
   ;;wavelength steps for the fitting process
   sig = dblarr(nw)+1.0
   sig[4:12] = 0.4
   sig[11]=5
   sig[0] = 5
   sig[-1] = 5
   
   ;ptot=0.72           ; DON'T CHANGE IT defined only to equilize the Icont for all lines
                                ; according to cont 4000
   ptot = 1.
 endif

calib = fltarr(3,nt)

for fr = 0L, nt-1 do begin

   ;; Loads the desired time step
   cube=d[4*fr]


   ;;calculates the average value
   ;;over the defined QS area
   ;;in each wavelength
   imean = dblarr(nw)
   for ii=0,nw-1 do imean[ii] = mean(cube[x0:x1, y0:y1,ii,0])
      
   ;; Extract the range of wavelength and corresponding intensity
   ;; and put it in x , y respectively "FROM SOLAR ATLAS" + to cgs
   red_satlas, w[0]+cw-1.0,  w[-1]+cw+1.0, x, y, /cgs
   x-= cw       ; substituting the cw to get the difference from the line center position

   fpi = cfpi(cw)                     ; object related to instrument profile
   dw = x[1]-x[0]                     ; wave length range
   np = 51                            ; ???
   tw = (dindgen(np)-np/2)*dw         ; adjusting the number positions from satlas to 51
   tr = fpi->dual_fpi(tw+cw, erh=-0.015) ;;it returns a structure with the nominal reflectivities
   ;;and cavity separations of the CRISP spectropolarimeter
   ;;USING THE NEW POSITIONS FROM TW
   ycon = fftconvol(y*ptot, tr/total(tr)) ; convolves the satlas intensity to a broader profile
   
   fpar = replicate({limits:[0d0,0d0], limited:[1,1], fixed:0}, 3)
   fpar[0].limits[*] = [0d0, mean(imean)/mean(ycon)*10d0]
   fpar[1].limits[*] = [-1d0,1.0d0]
   fpar[2].limits[*] = [0.98d0,1.02d0]
   fpar[2].fixed=1
   farg = {x:x, y:ycon, w:w, sig:sig, imean:imean}
   
   par = [double(mean(imean)/mean(ycon)), 0.001d0, 1.0d0]
   print, par
   par = mpfit('fitme', par, parinfo=fpar, fun=farg)
   calib[*,fr] = par

   ; plot the fitting results
   ;cgwindow
   ;plot, x, ycon, /line, xrange=[-2,2], /xstyle, ysty=3
   ;oplot, (w+par[1])*par[2], imean/par[0]
   ;cgplot, x+cw, ycon, linestyle = 1,/add,xrange = [min(x+cw),max(x+cw)]
   ;cgplot, (w+par[1]+cw)*par[2], imean/par[0],/add,/over ;, psym = 16, symsize = 0.75
   ;cglegend,titles = ['S-Atlas (cgs units)','Fit'],linestyles = ['1','0'],/add


   ; save the plot and fitting parameters
   ;cgcontrol,create_pdf = savedir+'calib_'+pref+'_'+string(fr, format='(I3.3)')+'.pdf'
   
   ;stop

endfor
writefits, savedir+'calib_'+pref+'_all.fits', calib

end
