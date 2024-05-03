
function fitme, p, x=x, y=y, imean=imean, w=w, sig=sig

  y2 = interpol(y, x-p[1], w*p[2])*p[0]
  ;plot, w, imean, /line
  ;oplot, w, y2
  
  return, (y2-imean)*sig
end

cgdelete, /all

fn = 28
dir='/scratch/sepid/DATA/AR/plage/2018.07.22/08:23:57/'
savedir = '/nadir-scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/'

if(1) then begin
   w=(f0(dir+'wav_3950.f0'))[0:20]*1e10 ; extract wavelength of the data and convert the unit
   w= w-w[10]                   ; the wavelength positions from the wavelength center
   fil0 = dir+'crispex_3950_*corrected.fcube' ; data cube
   fil1 = dir+'crispex_3950_*corrected_sp.fcube' ; spectral cube
   mapdatastokes,fil0,fil1,d,nx,ny,nt,nw
   pref = '3934'
   cw = 3933.6640d              ; labatory value for the wavelength center

   ; quiet sun patch coords
   x0=900
   y0=100
   xx=1100
   yy=300

   ; wavelength position weight to fit -> sigma
   sig = dblarr(21)+1.0 ; container for sigma
   sig[3:6] = 0.4
   sig[8:11]=0.2
   sig[13:15] = 0.4
   ;sig[5:14] = 0.4              ; wings
   ;sig[8:10]=0.2                ; around line core
   ptot=1
    ;ptot = 0.428727272817914
    ;ptot = 0.1
   cube=d[fn] ; the desired frame of the data
   imean = dblarr(21)           ; container for the mean intensity profile
   
   ; visualizing the selected quiet patch on the FOV
   cgwindow
   cgimage,cube[*,*,-1],/scal,/add,/keep,/axis
   tvbox, [xx-x0, yy-y0], (x0+xx)/2., (y0+yy)/2., color = 'red', /add, /data,thick = 2.5

   ; mean intensity profile over the quiet patch
   for ii=0,20 do imean[ii] = mean(cube[x0:xx, y0:yy,ii])
    
endif

; extracting the profile from Solar Atlas
red_satlas, w[0]+cw-5.0,  w[-1]+cw+5.0, x, y, /cgs ; x -> wave-length, y -> I in [cgs]
x-= cw ; change the atlas wavelength value into distance from the line-core
dw = x[1]-x[0]                  ; S-atlas delta-wavelength

; Convolution of the S-atlas with Chromis profile
np = round((0.080 * 8) / dw)
tw = (dindgen(np)-np/2)*dw
tr = chromis_profile(tw+cw, erh=-0.07d0)
tr /= total(tr)
y1 = fftconvol(y*ptot, tr) ; convolved S-atlas intensity

; fitting preparation
fpar = replicate({limits:[0d0,0d0], limited:[1,1], fixed:0}, 3)
fpar[0].limits[*] = [0d0, mean(imean)/mean(y1)*10d0]
fpar[1].limits[*] = [-1d0,1.0d0]
fpar[2].limits[*] = [0.95d0,1.08d0]
fpar[2].fixed=1
farg = {x:x, y:y1, w:w, sig:sig, imean:imean}

par = [double(mean(imean)/mean(y1)), 0.001d0, 1.0d0]
print, par
par = mpfit('fitme', par, parinfo=fpar, fun=farg) ; fitting

; visual result of the fitting
cgwindow
cgplot, x, y1, linestyle = 1,/add, xrange=[-1,1] ; S-Atlas
cgplot, (w+par[1])*par[2], imean/par[0],/add,/over ; Calibrated data
cglegend,titles = ['S-Atlas (cgs units)','Fit'],linestyles = ['1','0'],/add

;cgcontrol,create_pdf = savedir+'calib.'+pref+'.pdf'
;writefits, dir+'calib.'+pref+'.fits', par
stop
end
