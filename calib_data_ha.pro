function fitme, p, x=x, y=y, imean=imean, w=w, sig=sig

  y2 = interpol(y, x-p[1], w*p[2])*p[0]
  plot, w, imean, /line
  oplot, w, y2
  
  return, (y2-imean)*sig
end
cgdelete,/all
fn = 28
dir='/nadir-scratch/sepid/DATA/AR/plage/2016.09.15/'
savedir = '/nadir-scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/'

;restore,'slit_3950.csav'
w=f0(dir+'wav.6563.f0') ;*1e10
;w= w-w[7]
fil0 = dir+'crispex.6563*.fcube'
fil1 = dir+'crispex.6563*sp.fcube'
mapdatastokes,fil0,fil1,d,nx,ny,nt,nw
pref = '6563'
cw = 6562.7963d
;coords of a quiet patch
x0=900
y0=68
xx=1100
yy=260
;specifing weights for fitting the profile
sig = dblarr(16)+1.0
sig[1:3] = 0.4
sig[4:10] = 0.2
sig[11:13] = 0.4
ptot=1
;ptot = 0.428727272817914
;ptot = 0.1
cube=d[fn]
imean = dblarr(16)
cgwindow
cgimage,cube[*,*,-1],/scal,/add,/keep,/axis
tvbox, [xx-x0, yy-y0], (x0+xx)/2., (y0+yy)/2., color = 'red', /add, /data,thick = 2.5

;stop
for ii=0,15 do imean[ii] = mean(cube[x0:xx, y0:yy,ii])

;extracting from the solar atlas 
red_satlas, w[0]+cw-5.0,  w[-1]+cw+5.0, x, y, /cgs
x-= cw

dw = x[1]-x[0]
np = round((0.080 * 8) / dw)
tw = (dindgen(np)-np/2)*dw
;tr = chromis_profile(tw+cw, erh=-0.07d0)
fpi = cfpi(cw) ; object related to instrument profile
tr = fpi->dual_fpi(tw+cw, erh=-0.015) ;;it returns a structure with the nominal reflectivities
                                                             ;;and cavity separations of the CRISP spectropolarimeter
tr /= total(tr)
y1 = fftconvol(y*ptot, tr)


fpar = replicate({limits:[0d0,0d0], limited:[1,1], fixed:0}, 3)
fpar[0].limits[*] = [0d0, mean(imean)/mean(y1)*10d0]
fpar[1].limits[*] = [-1d0,1.0d0]
fpar[2].limits[*] = [0.95d0,1.08d0]
fpar[2].fixed=1
farg = {x:x, y:y1, w:w, sig:sig, imean:imean}

par = [double(mean(imean)/mean(y1)), 0.001d0, 1.0d0]
print, par
par = mpfit('fitme', par, parinfo=fpar, fun=farg)

cgwindow
;plot, x, y1, /line, /xstyle, ysty=3, xrange=[-1,1]
;oplot, (w+par[1])*par[2], imean/par[0]
cgplot, x, y1, linestyle = 1,/add, xrange=[-1,1]
cgplot, (w+par[1])*par[2], imean/par[0],/add,/over
cglegend,titles = ['S-Atlas (cgs units)','Fit'],linestyles = ['1','0'],/add

cgcontrol,create_pdf = savedir+'test.calib.'+pref+'.pdf'
writefits, dir+'test.calib.'+pref+'.fits', par
stop
end
