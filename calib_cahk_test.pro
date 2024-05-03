cgdelete,/all
fn = 28
dir='/nadir-scratch/sepid/DATA/AR/plage/2016.09.15/'
savedir = '/nadir-scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/'

; test
ck = readfits('/scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/map_3950_fr28.fits')

if(1) then begin
   pref = '3934'
   cw = 3933.6640d              ; labatory value for the wavelength center
   ;w=(f0(dir+'wav_3950.f0'))[0:20]*1e10 ; extract wavelength of the data and convert the unit
   w=ck[0,0,0:20,0] ; extract wavelength of the data and convert the unit
   w= w-cw                   ; the wavelength positions from the wavelength center
   fil0 = dir+'crispex_3950_*.fcube' ; data cube
   fil1 = dir+'crispex_3950_*_sp.fcube' ; spectral cube
   mapdata,fil0,fil1,d,nx,ny,nt,nw

   ; quiet sun patch coords
   x0=10
   y0=10
   xx=1500
   yy=1000

   imean = dblarr(21)           ; container for the mean intensity profile
   
   ; visualizing the selected quiet patch on the FOV
   cgwindow
   cgimage,ck[*,*,-1, 1],/scal,/add,/keep,/axis
   tvbox, [xx-x0, yy-y0], (x0+xx)/2., (y0+yy)/2., color = 'red', /add, /data,thick = 2.5

   ; mean intensity profile over the quiet patch
   for ii=0,20 do imean[ii] = mean(ck[x0:xx, y0:yy,ii, 1]/1000.) ; test
    
endif

; extracting the profile from Solar Atlas
red_satlas, w[0]+cw-5.0,  w[-1]+cw+5.0, x, y, /cgs ; x -> wave-length, y -> I in [cgs]
x-= cw ; change the atlas wavelength value into distance from the line-core
dw = x[1]-x[0]                  ; S-atlas delta-wavelength

; Convolution of the S-atlas with Chromis profile
ptot=1.
np = round((0.080 * 8) / dw)
tw = (dindgen(np)-np/2)*dw
tr = chromis_profile(tw+cw, erh=-0.07d0)
tr /= total(tr)
y1 = fftconvol(y*ptot, tr) ; convolved S-atlas intensity

; visual result of the fitting
cgwindow
cgplot, x, y1, linestyle = 1,/add, xrange=[-1,1] ; S-Atlas
cgplot, w, imean,/add,/over ; test
cglegend,titles = ['S-Atlas (cgs units)','Fit'],linestyles = ['1','0'],/add

stop
end
