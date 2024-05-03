; data_3934.pro
; purpose:       (1) Converts the CHROMIS/crisp fcubes from QS to
; polarization maps
; (Here only applied to Full-Stokes data of chromosphere -> applicable
; for any other full-stokes parameters)
;                specified wavelenths + I continuum!
;                (4) Smooth with a box-car of 3*3 pixels
; Input:         fcube data
; Output:        polarization maps of each frame

close, /all

;variable to be defined
;read,sm_fact,prompt = 'smoothing factor = '
lcore = 3934
fcb = 12                         ; for full-stokes data of Ca II 8542

;Data and Output directory  
nadir = '/nadir-scratch/jaime/CHROMIS/'
QSdir = '2016.10.12_quietsun_10_46_37_cakfine_crisp/'
datadir = '/nadir-scratch/sepid/DATA/'

cd, nadir+QSdir
print,fcb

fcube,f

;data dimensions
lp_header,f[fcb],nx=nxx,ny=nyy   ;spatial dimensions
lp_header,f[fcb+1],nx=nww,ny=ntt ;number of spectral lines and frames

def = fltarr(nxx,nyy,ntt)
Ik = def
Ik_0 = def
Ik_3 = def
Ik_7 = def
Ik_11 = def
Ik_13 = def
Ik_14 = def
Ik_16 = def
Ik_19 = def
Ik_21 = def
Ik_23 = def
Ik_27 = def
Ik_32 = def
Ik_34 = def
Icont400 = def

;Reading the data
for i=0L , ntt-1 do begin

Ik(*,*,i) = lp_get(f[fcb], i*nww + 17)
Ik_0(*,*,i) = lp_get(f[fcb], i*nww + 0)
Ik_3(*,*,i) = lp_get(f[fcb], i*nww + 3)
Ik_7(*,*,i) = lp_get(f[fcb], i*nww + 7)
Ik_11(*,*,i) = lp_get(f[fcb], i*nww + 11)
Ik_13(*,*,i) = lp_get(f[fcb], i*nww + 13)
Ik_14(*,*,i) = lp_get(f[fcb], i*nww + 14)
Ik_16(*,*,i) = lp_get(f[fcb], i*nww + 16)
Ik_19(*,*,i) = lp_get(f[fcb], i*nww + 19)
Ik_21(*,*,i) = lp_get(f[fcb], i*nww + 21)
Ik_23(*,*,i) = lp_get(f[fcb], i*nww + 23)
Ik_27(*,*,i) = lp_get(f[fcb], i*nww + 27)
Ik_32(*,*,i) = lp_get(f[fcb], i*nww + 32)
Ik_34(*,*,i) = lp_get(f[fcb], i*nww + 34)
Icont400(*,*,i) = lp_get(f[fcb], i*nww + 35)

print,'...frame no. '+ strtrim(long(i),2) + ' processed'
;stop
endfor

print,'STORING MAPS IN:'
print,datadir
print,' '

save,Ik,Ik_0,Ik_3,Ik_7,Ik_11,Ik_13,Ik_14,Ik_16,Ik_19,Ik_21,Ik_23,Ik_27,Ik_32,Ik_34,Icont400,$
filename = datadir+strtrim(long(lcore),2)+'_intensity_maps.sav'

stop

end
