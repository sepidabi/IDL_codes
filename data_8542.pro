; data_8542.pro
; purpose:       (1) Converts the CHROMIS/crisp fcubes from QS to
; polarization maps
; (Here only applied to Full-Stokes data of chromosphere -> applicable
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
read,sm_fact,prompt = 'smoothing factor = '
lcore = 8542
fcb = 10                         ; for full-stokes data of Ca II 8542

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
nss = 4                          ;number Stokes parameters

;Area with ACTUAL DATA (visually selected) => ONLY FOR SPECIFING NOISE LEVEL
x1 = 295 & x2 = 1635
y1 = 0   & y2 = 1040

def = fltarr(nxx,nyy,ntt)
CP=def
LP=def
CPcont=def
LPcont=def
Icont=def
I854 = def
I854_13 = def
I854_16 = def
I854_19 = def
I854_7 = def
I854_4 = def 
I854_1 = def

sigmaLP = fltarr(ntt)
sigmaCP = fltarr(ntt)

;Reading the data
for i=0L , ntt-1 do begin

   Ic = lp_get(f[fcb], i*nww*nss + 0*nww + 20)
I854(*,*,i) 	= lp_get(f[fcb], i*nww*nss + 0*nww + 10)
I854_13(*,*,i) 	= lp_get(f[fcb], i*nww*nss + 0*nww + 13)
I854_16(*,*,i) 	= lp_get(f[fcb], i*nww*nss + 0*nww + 16)
I854_19(*,*,i) 	= lp_get(f[fcb], i*nww*nss + 0*nww + 19)
I854_7(*,*,i) 	= lp_get(f[fcb], i*nww*nss + 0*nww + 7)
I854_4(*,*,i) 	= lp_get(f[fcb], i*nww*nss + 0*nww + 4)
I854_1(*,*,i) 	= lp_get(f[fcb], i*nww*nss + 0*nww + 1)
   Qc = lp_get(f[fcb], i*nww*nss + 1*nww + 20)
   Uc = lp_get(f[fcb], i*nww*nss + 2*nww + 20)
   Vc = lp_get(f[fcb], i*nww*nss + 3*nww + 20)
;stop
   ;Initial Value of LP and CP
   LPt_l = fltarr(nxx,nyy)
   LPt_r = fltarr(nxx,nyy)
   CPt_l = fltarr(nxx,nyy)
   CPt_r = fltarr(nxx,nyy)
   StokesI = fltarr(nxx,nyy)

   ;Reform and Smoothing Continuum maps
   Icont(*,*,i) = smooth(reform(Ic), [sm_fact,sm_fact], /edge_truncate)
   Qcont = smooth(reform(Qc), [sm_fact,sm_fact], /edge_truncate)
   Ucont = smooth(reform(Uc), [sm_fact,sm_fact], /edge_truncate)
   Vcont = smooth(reform(Vc), [sm_fact,sm_fact], /edge_truncate)
   
;=================================================================================
;     chromosphere        chromosphere        chromosphere        chromosphere
;=================================================================================

   ;Total CP over the left-side of "selected" spectral positions
   for l = 0L, 8 do begin
      del = 1.
      StV_l = lp_get(f[fcb], i*nww*nss + 3*nww + l+1)
      CPt_l = CPt_l + del*smooth(reform(StV_l), [sm_fact,sm_fact], /edge_truncate)

      StQ = lp_get(f[fcb], i*nww*nss + 1*nww + l)
      StU = lp_get(f[fcb], i*nww*nss + 2*nww + l)
      LPt_l = LPt_l + sqrt(smooth(reform(StQ), [sm_fact,sm_fact], /edge_truncate)^2 + smooth(reform(StU), [sm_fact,sm_fact], /edge_truncate)^2)

   endfor

   ;Total CP over the right-side of "selected" spectral positions
   for r = 0L, 8 do begin
      del = -1.
      StV_r = lp_get(f[fcb], i*nww*nss + 3*nww + r+11)
      CPt_r = CPt_r + del*smooth(reform(StV_r), [sm_fact,sm_fact], /edge_truncate)
      StQ = lp_get(f[fcb], i*nww*nss + 1*nww + l)
      StU = lp_get(f[fcb], i*nww*nss + 2*nww + l)
      LPt_r = LPt_r + sqrt(smooth(reform(StQ), [sm_fact,sm_fact], /edge_truncate)^2 + smooth(reform(StU), [sm_fact,sm_fact], /edge_truncate)^2)

   endfor

   ;Averaging over the 18 wavelengths around line core (8542)
   LP(*,*,i) = float(LPt_l+LPt_r)/float(18.*Icont)
   CP(*,*,i) = float(CPt_l+CPt_r)/float(18.*Icont)
   LPcont(*,*,i) = float(sqrt(Qcont^2+Ucont^2))/float(Icont)
   CPcont(*,*,i) = float(Vcont)/float(Icont)
   sigmaCP(i) = stddev(reform(CPcont[x1:x2,y1:y2,i]))/sqrt(18.)
   sigmaLP(i) = stddev(reform(LPcont[x1:x2,y1:y2,i]))/sqrt(18.)

print,'...frame no. '+ strtrim(long(i),2) + ' processed'
;stop
endfor

print,'STORING MAPS IN:'
print,datadir
print,' '
LPch = LP
CPch = CP
LPcontch = LPcont
CPcontch = CPcont
Icontch = Icont
sigmaCPch = sigmaCP

if sm_fact eq 1 then begin
	save,LPch,CPch,LPcontch,sigmaCPch,CPcontch,Icontch,sm_fact,$
I854,I854_13,I854_16,I854_19,I854_7,I854_4,I854_1,$
filename = datadir+strtrim(long(lcore),2)+'_polarization_maps.sav'
endif else begin
	save,LPch,CPch,LPcontch,sigmaCPch,CPcontch,Icontch,sm_fact,$
I854,I854_13,I854_16,I854_19,I854_7,I854_4,I854_1,$
filename = datadir+strtrim(long(lcore),2)+'_polarization_maps__smoothed.sav'
endelse

stop

end
