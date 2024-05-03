pro stokes_6302
; data_6302.pro
; purpose:       (1) Converts the CHROMIS/crisp fcubes from QS to
; polarization maps
; (Here only applied to Full-Stokes data of photosphere -> applicable
; for any other full-stokes parameters)
;                (2) Calculates the total LP and CP averaged over the
;                specified wavelenths + I continuum!
;                (3) Normalizes the maps to the I continuum
;                (4) Smooth with a box-car of 3*3 pixels
; Input:         fcube data
; Output:        polarization maps of each frame
; noise value of the maps

close, /all

;variable to be defined
sm_fact = 1
lcore = 6302
fcb = 8                         ; for full-stokes data of Fe I 6302

;Data and Output directory  
datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
data = [datadir+ '*6302*.fcube']
stop
;data dimensions
lp_header,f[fcb],nx=nxx,ny=nyy   ;spatial dimensions
lp_header,f[fcb+1],nx=nww,ny=ntt ;number of spectral lines and frames
nss = 4                          ;number Stokes parameters

;area with ACTUAL DATA (visually selected) => ONLY FOR CALCULATING THE NOISE LEVEL
x1 = 210 & x2 = 1630
y1 = 0   & y2 = 1085

default = fltarr(nxx,nyy,ntt)
CP=default & LP=default & CPcont=default & LPcont=default & Icont=default

sigmaLP = fltarr(ntt)
sigmaCP = fltarr(ntt)

;Reading the data
for i=0L , ntt-1 do begin

   Ic = lp_get(f[fcb], i*nww*nss + 0*nww + 9)
   Qc = lp_get(f[fcb], i*nww*nss + 1*nww + 9)
   Uc = lp_get(f[fcb], i*nww*nss + 2*nww + 9)
   Vc = lp_get(f[fcb], i*nww*nss + 3*nww + 9)
;stop
   ;Reform and Smoothing Continuum maps
   Icont(*,*,i) = smooth(reform(Ic), [sm_fact,sm_fact], /edge_truncate)
   Qcont = smooth(reform(Qc), [sm_fact,sm_fact], /edge_truncate)
   Ucont = smooth(reform(Uc), [sm_fact,sm_fact], /edge_truncate)
   Vcont = smooth(reform(Vc), [sm_fact,sm_fact], /edge_truncate)
   
;==============================================================================
;        PHOTOSPHERE        PHOTOSPHERE        PHOTOSPHERE        PHOTOSPHERE
;==============================================================================

   ;Initial Value of LP and CP
   LPt_l = fltarr(nxx,nyy)
   LPt_r = fltarr(nxx,nyy)
   CPt_l = fltarr(nxx,nyy)
   CPt_r = fltarr(nxx,nyy)
   StokesI = fltarr(nxx,nyy)

   ;IN CASE OF FE I:
   ;Total CP over the left-side of "selected" spectral positions
   for l = 0L, 3 do begin
      del = 1.
      StV_l = lp_get(f[fcb], i*nww*nss + 3*nww + l)
      CPt_l = CPt_l + del*smooth(reform(StV_l), [sm_fact,sm_fact], /edge_truncate)

      StQ = lp_get(f[fcb], i*nww*nss + 1*nww + l)
      StU = lp_get(f[fcb], i*nww*nss + 2*nww + l)
      LPt_l = LPt_l + sqrt(smooth(reform(StQ), [sm_fact,sm_fact], /edge_truncate)^2 + smooth(reform(StU), [sm_fact,sm_fact], /edge_truncate)^2)

   endfor

   ;Total CP over the right-side of "selected" spectral positions
   for r = 0L, 3 do begin
      del = -1.
      StV_r = lp_get(f[fcb], i*nww*nss + 3*nww + r+5)
      CPt_r = CPt_r + del*smooth(reform(StV_r), [sm_fact,sm_fact], /edge_truncate)

      StQ = lp_get(f[fcb], i*nww*nss + 1*nww + l)
      StU = lp_get(f[fcb], i*nww*nss + 2*nww + l)
      LPt_r = LPt_r + sqrt(smooth(reform(StQ), [sm_fact,sm_fact], /edge_truncate)^2 + smooth(reform(StU), [sm_fact,sm_fact], /edge_truncate)^2)

   endfor

      ;StI = lp_get(f[fcb], i*nww*nss + 0*nww + l)
      ;StokesI = StokesI + smooth(reform(StI), [sm_fact,sm_fact], /edge_truncate)
      ;Total LP over the "selected" spectral positions
      
      

   ;Averaging over the 8 wavelengths around the 1st line core (630.1 nm)
   CP(*,*,i) = float(CPt_l+CPt_r)/float(8.*Icont)
   LP(*,*,i) = float(LPt_l+LPt_r)/float(8.*Icont)
   CPcont(*,*,i) = float(Vcont)/float(Icont)
   LPcont(*,*,i) = float(sqrt(Qcont^2+Ucont^2))/float(Icont)
   sigmaCP(i) = stddev(reform(CPcont[x1:x2,y1:y2,i]))/sqrt(8.)
   sigmaLP(i) = stddev(reform(LPcont[x1:x2,y1:y2,i]))/sqrt(8.)

   print,'...frame no. '+ strtrim(long(i),2) + ' processed'

;stop

endfor

print,'STORING MAPS IN:'
print,datadir
print,' '

if sm_fact eq 1 then begin
	save,LP,CP,LPcont,sigmaLP,sigmaCP,CPcont,Icont,sm_fact,filename = datadir+strtrim(long(lcore),2)+'_polarization_maps.sav'
endif else begin
	save,LP,CP,LPcont,sigmaLP,sigmaCP,CPcont,Icont,sm_fact,filename = datadir+strtrim(long(lcore),2)+'_polarization_maps__smoothed.sav'
endelse

stop

end
