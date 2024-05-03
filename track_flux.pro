;To track the detected features from "detect_flux.pro"
;keywords same as detection code

pro track_flux,photo=photo,chrom=chrom,nth1=nth1,nth2=nth2,lp=lp,cp=cp

close, /all

sm_fact = 3.

c = 1.
size_threshold = 20.
  
nadir = '/nadir-scratch/jaime/CHROMIS/'
QSdir = '2016.10.12_quietsun_10_46_37_cakfine_crisp/'
AnaDir = '/home/seki2695/PROJECTS/Flux_Em/Analyze/'
TrDir = '/home/seki2695/PROJECTS/Flux_Em/Track/'
outdir = '/home/seki2695/PROJECTS/Flux_Em/Outputs/'

cd, nadir+QSdir

fcube,f

;--------------------------------------------------------------------------
;  0 Bazi_ME_6301_6302_CHROMIS.fcube
;  1 Bhor_ME_6301_6302_CHROMIS.fcube
;  2 Bincl_ME_6301_6302_CHROMIS.fcube
;  3 Bstrength_ME_6301_6302_CHROMIS.fcube
;  4 Bz_ME_6301_6302_CHROMIS.fcube
;  5 Vlos_ME_6301_6302_CHROMIS.fcube
;  6 crispex.6563.10:46:35.time_corrected_CHROMIS.fcube
;  7 crispex.6563.10:46:35.time_corrected_CHROMIS_sp.fcube
;  8 crispex.stokes.6302.10:46:35.time_corrected_CHROMIS.fcube
; 9 crispex.stokes.6302.10:46:35.time_corrected_CHROMIS_sp.fcube
; 10 crispex.stokes.8542.10:46:35.time_corrected_CHROMIS.fcube
; 11 crispex.stokes.8542.10:46:35.time_corrected_CHROMIS_sp.fcube
; 12 crispex_3950_2016-10-12T10:46:37_scans=0-44_time-corrected.fcube
; 13 crispex_3950_2016-10-12T10:46:37_scans=0-44_time-corrected_sp.fcube
; 14 wb_3950_2016-10-12T10:46:37_scans=0-44_corrected.icube
;--------------------------------------------------------------------------

if keyword_set(photo) then fcb = 8 else fcb = 10
if keyword_set(lp) then st = 1 else st = 2

;data dimensions
lp_header,f[fcb],nx=nxx,ny=nyy   ;spatial dimensions
lp_header,f[fcb+1],nx=nww,ny=ntt ;number of spectral lines and frames
nss = 4                          ;number Stokes parameters

;selecting the area with ACTUAL DATA (visually selected)
x1 = 195 & x2 = 1615
y1 = 0   & y2 = 1085

;Defining Cubes for the time-series data
default = fltarr(x2-x1+1,y2-y1+1,ntt)

Icont = default  &  LP = default  &  CP = default  &  CPf_l = default  &  CPf_h = default
im = lonarr(x2-x1+1,y2-y1+1,ntt)

;Reading the data
for i=0L , ntt-1 do begin


   print,'  Reading... CHROMIS DATA - frame: '+strtrim(long(i),2)
   print,'====================================='
   
   ;Getting initial frame of the maps
   StI = lp_get(f[fcb], i*nww*nss + 0*nww + 0)
   StQ = lp_get(f[fcb], i*nww*nss + 1*nww + 0)
   StU = lp_get(f[fcb], i*nww*nss + 2*nww + 0)
   StV = lp_get(f[fcb], i*nww*nss + 3*nww + 0)

   Ic = lp_get(f[fcb], i*nww*nss + 0*nww + 9)
   Qc = lp_get(f[fcb], i*nww*nss + 1*nww + 9)
   Uc = lp_get(f[fcb], i*nww*nss + 2*nww + 9)
   Vc = lp_get(f[fcb], i*nww*nss + 3*nww + 9)

   ;cropping the edges and smoothing the data
   StokesI = smooth(reform(StI[x1:x2,y1:y2]), [sm_fact,sm_fact], /edge_truncate)
   StokesQ = smooth(reform(StQ[x1:x2,y1:y2]), [sm_fact,sm_fact], /edge_truncate)
   StokesU = smooth(reform(StU[x1:x2,y1:y2]), [sm_fact,sm_fact], /edge_truncate)
   StokesV = smooth(reform(StV[x1:x2,y1:y2]), [sm_fact,sm_fact], /edge_truncate)

   ;Initial Value of LP and CP
   LPt = sqrt(StokesQ^2 + StokesU^2)
   CPt = StokesV

   ;Continuum maps
   Icont(*,*,i) = smooth(reform(Ic[x1:x2,y1:y2]), [sm_fact,sm_fact], /edge_truncate)
   Qcont = smooth(reform(Qc[x1:x2,y1:y2]), [sm_fact,sm_fact], /edge_truncate)
   Ucont = smooth(reform(Uc[x1:x2,y1:y2]), [sm_fact,sm_fact], /edge_truncate)
   Vcont = smooth(reform(Vc[x1:x2,y1:y2]), [sm_fact,sm_fact], /edge_truncate)

   ;Noise
   LPcont = float(sqrt(Ucont^2+Vcont^2))/float(Icont)
   CPcont = float(Vcont)/float(Icont)
   
;==============================================================================
;        PHOTOSPHERE        PHOTOSPHERE        PHOTOSPHERE        PHOTOSPHERE
;==============================================================================
if fcb eq 8 then begin

   ;IN CASE OF FE I:
   ;Averaging over the 8 wavelengths around the 1st line core (6301 nm)
   for l = 0L, nww-7 do begin

      StI = lp_get(f[fcb], i*nww*nss + 0*nww + l)
      StokesI = StokesI + smooth(reform(StI[x1:x2,y1:y2]), [sm_fact,sm_fact], /edge_truncate)

      if l le 4 then begin
         del = 1.
      endif else begin
         del = -1.
      endelse

      ;Total CP over the "selected" spectral positions
      StV = lp_get(f[fcb], i*nww*nss + 3*nww + l)
      CPt = CPt + del*smooth(reform(StV[x1:x2,y1:y2]), [sm_fact,sm_fact], /edge_truncate)

      ;Total CP over the "selected" spectral positions
      StQ = lp_get(f[fcb], i*nww*nss + 1*nww + l)
      StU = lp_get(f[fcb], i*nww*nss + 2*nww + l)
      LPt = LPt + sqrt(smooth(reform(StQ[x1:x2,y1:y2]), [sm_fact,sm_fact], /edge_truncate)^2 + smooth(reform(StU[x1:x2,y1:y2]), [sm_fact,sm_fact], /edge_truncate)^2)
      
   endfor

   ;Averaging LP and CP over spectral positions
   LP(*,*,i) = LPt/float(8.*Icont)
   CP(*,*,i) = CPt/float(8.*Icont)

   ;Noise Level -> based on the number of spectral positions
   sigmaLP = stddev(LPcont)/sqrt(8.)
   sigmaCP = stddev(CPcont)/sqrt(8.)

   ;dimension of cropped image
   im_size = size(CP)
   nnx = im_size[1]  &  nny = im_size[2]
   
   print,'  Reading... DETECTION RESULTS - frame: '+strtrim(long(i),2)
   print,'=========================================='   

   ;LP tracking
   if st eq 1 then begin
      
print,'Tracking CPfs...'
print,'================'

   ;CP tracking
   endif else begin

      CPf_l(*,*,i) = readfits(AnaDir+'PHOTOSPHERE/CP/CP__nth_'+strtrim(long(nth1),2)+'__frame_'+strtrim(long(i),2)+'__features.fits')
      CPf0_l = CPf_l(*,*,0)

      CPf_h(*,*,i) = readfits(AnaDir+'PHOTOSPHERE/CP/CP__nth_'+strtrim(long(nth2),2)+'__frame_'+strtrim(long(i),2)+'__features.fits')
      CPf0_h = CPf_h(*,*,0)
      


   endelse
      
;==============================================================================
;      CHROMOSPHERE       CHROMOSPHERE       CHROMOSPHERE       CHROMOSPHERE
;==============================================================================

endif else begin

;Ca I 8542 feature detection
   
endelse
;stop

endfor

;====================================================
;		START TRACKING PHASE
;====================================================

print,'Tracking CPfs...'
print,'================'

	; parameters needed for tracking
	displ= 2 ; maximum displacement allowed -- maybe need to be modified!
	CPsize = 21 ; an odd number, roughly two times of typical size of patches in pixel -- maybe need to be modified!
	intensityLim = 2 ; a pixel value threshold -- maybe need to be modified!

	; Track and Measure the following parameters:
	;(1) x-centroid (2) y-centroid (3) total brightness (4) radius of gyration and (5) eccentricity
	f_l = feature(CPf0_l,CPsize)
	help, f_l
        prefix_l = AnaDir+'pretrack__CPfs__nth_'+strtrim(long(nth1),2)+'.gdf'

        epretrack_me, CPf_l, bplo=1, bphi=BPsize, dia=BPsize, mass=intensityLim, prefix=prefix_l, /nobpass

	pt_l = read_gdf(prefix_l)
	help, pt_l

	t_l = track_mbe(pt_l,displ,goodenough=1,memory=2)
	help, t_l

        nFOI_l = n_elements(t_l(0,*))
	id = t_l(6,*)
	frame_n = t_l(5,*)
	xxi_l = t_l(0,*)
	yyi_l = t_l(1,*)
	bri_l = t_l(2,*)
	ssize_l = t_l(3,*)
	lastid_l = id_l(nFOI_l-1)
        
stop
end
