;The program is to detect the magnetic features on CHROMIS QS area.
;
;keywords:
; /photo -> investigate the Photosphere layer (Fe I 6302)
; /chrom -> investigate the low Chromosphere layer (Ca I 8542)
; cp -> for cpf detection
; lp -> for lp detection
; nth1 -> lower signal threshold factor
; nth2 -> higher signal threshold factor
; 

pro detect_flux, photo = photo, chrom = chrom, lp = lp, cp = cp, nth1 = nth1, nth2 = nth2

  close, /all

  sm_fact = 3.

  cc = 1.
  size_threshold = 20.
  
nadir = '/nadir-scratch/jaime/CHROMIS/'
QSdir = '2016.10.12_quietsun_10_46_37_cakfine_crisp/'
AnaDir = '/home/seki2695/PROJECTS/Flux_Em/Analyze/'
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

;Reading the data
for i=0L , ntt-1 do begin

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
   Icont = smooth(reform(Ic[x1:x2,y1:y2]), [sm_fact,sm_fact], /edge_truncate)
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
   LP = LPt/float(8.*Icont)
   CP = CPt/float(8.*Icont)

   ;Noise Level -> based on the number of spectral positions
   sigmaLP = stddev(LPcont)/sqrt(8.)
   sigmaCP = stddev(CPcont)/sqrt(8.)

   ;dimension of cropped image
   im_size = size(CP)
   nnx = im_size[1]  &  nny = im_size[2]
   image = fltarr(nnx,nny )     ;to be used either as LP or CP

;If LP selected
if st eq 1 then begin

   file_mkdir, AnaDir+'PHOTOSPHERE/LP/frame_'+strtrim(long(i),2)
   FDir = AnaDir+'PHOTOSPHERE/LP/frame_'+strtrim(long(i),2)+'/'   

   aLP = ABS(LP)                ; absolut values of LP >>> for event detection

   ;Loop to have both threshold levels      
   iLP = where(aLP ge nth1*sigmaLP)

print,'===================================================='	
print,'+++++++++++++++++++ LP Detection +++++++++++++++++++'
print,'===================================================='

   image(iLP) = aLP(iLP)
	
   image=bytscl(image)

   ; Define a structuring kernel for an opening operation on the image.
   radius = 1
   kernel = SHIFT(DIST(2*radius), radius, radius) LE radius

   ; Apply the opening operator to the image.
   openImage = MORPH_OPEN(image, kernel, /GRAY)

; Threshold the image 
   mask = openImage GE 0

; Do the analysis (LPF = LP Features)
   LPF = Obj_New('blob_analyzer', openImage, MASK=mask, SCALE=[1.0, 1.0])

   count = LPF -> NumberOfBlobs()

	print, ' '
	print, '  -----------------------------------------------'
	print, ' '
	print, ' >>> Number of initially detected LP features: ',float(count),format='(a,f6.0,$)' 
	print, ' '
	print, '  -----------------------------------------------'
	print, ' '


	; Report and save stats
	LPF -> ReportStats
	
	filestats = AnaDir+'stat__frame_'+strtrim(long(i),2)+'.txt'
        LPF -> ReportStats, file=filestats, tofile=1, /noscale
        
	; create and save a new frame that only includes the detected events: to be used for, e.g., tracking procedure
	blank = fltarr(nxx,nyy)

	; Display the CPF candidates
	s = Size(image, /DIMENSIONS)
	Window, 0, XSIZE=s[0]/cc, YSIZE=s[1]/cc, xpos=1925, ypos=1000
	nx=s[0]   &   ny=s[1]
        tvscl, histo_opt(congrid(LP, nx/cc, ny/cc, cubic=-0.5))
        
	; Display the LPFs we located with LABEL_REGION
	counto = 0

        file_mkdir, AnaDir+'PHOTOSPHERE/LP/frame_'+strtrim(long(i),2)
        FDir = AnaDir+'PHOTOSPHERE/LP/frame_'+strtrim(long(i),2)+'/'
        
        FOR j=0,count-1 DO BEGIN
           
           stats = LPF -> GetStats(j, /NoScale, XYINDICES=XYINDICES, INDICES=INDICES)
           
           if stats.count ge size_threshold then begin
              
		PLOTS, stats.perimeter_pts[0,*]/cc, stats.perimeter_pts[1,*]/cc, /Device, COLOR=cgColor('Red'), thick=2
		XYOUTS, (stats.center[0]/cc)+5., (stats.center[1]/cc)-5, /Device, StrTrim(j,2), $
		    COLOR=cgColor('Navy'), ALIGNMENT=0., CHARSIZE=1.3, CHARTHICK=2
			
		; save coordinates and indices of individual detected features
		outsc = fltarr(3,n_elements(INDICES))
		outnamecoords = FDir+'frame_'+strtrim(long(i),2)+'__'+'_feature_'+strtrim(long(j),2)+'.txt'
                outsc(0,*)=XYINDICES(0,*)  &  outsc(1,*)=XYINDICES(1,*)  &  outsc(2,*)=INDICES
                write_table, outsc, outnamecoords
			
		blank[INDICES] = image[INDICES]
             endif
        ENDFOR
        
        saveimage, AnaDir+'PHOTOSPHERE/LP/LP__frame_'+strtrim(long(i),2)+'__features.png',quality =100,png=png
        writefits, AnaDir+'PHOTOSPHERE/LP/LP__frame_'+strtrim(long(i),2)+'__features.fits', blank

print, ' '
print, ' >>> frame '+strtrim(i,2)+' is analysed :)'
print, ' '
print, ' >>> type .c to continue with the next frame'
print, ' '

endif

if st eq 2 then begin

   file_mkdir, AnaDir+'PHOTOSPHERE/CP/frame_'+strtrim(long(i),2)
   FDir = AnaDir+'PHOTOSPHERE/CP/frame_'+strtrim(long(i),2)+'/'

   aCP = ABS(CP)                ; absolut values of LP >>> for event detection

   ;To have both threshold levels
   iCP_L = where(aCP ge nth1*sigmaCP)
   iCP_H = where(aCP ge nth2*sigmaCP)
        
   print,'===================================================='	
   print,'+++++++++++++++++++ CP Detection +++++++++++++++++++'
   print,'===================================================='

   image_l = fltarr(nnx,nny)
   image_h = fltarr(nnx,nny)
   image_l(iCP_l) = aCP(iCP_l)
   image_h(iCP_h) = aCP(iCP_h)
	
   image_l=bytscl(image_l)
   image_h=bytscl(image_h)

   ; Define a structuring kernel for an opening operation on the image.
   radius = 1
   kernel = SHIFT(DIST(2*radius), radius, radius) LE radius

   ; Apply the opening operator to the image.
   openImage_l = MORPH_OPEN(image_l, kernel, /GRAY)
   openImage_h = MORPH_OPEN(image_h, kernel, /GRAY)

; Threshold the image 
   mask_l = openImage_l GE 0
   mask_h = openImage_h GE 0
; Do the analysis (LPF = LP Features)
   CPF_l = Obj_New('blob_analyzer', openImage_l, MASK=mask, SCALE=[1.0, 1.0])
   CPF_h = Obj_New('blob_analyzer', openImage_h, MASK=mask, SCALE=[1.0, 1.0])


   ;LOW THRESHOLD CASE
   ;==================
   count_l = CPF_l -> NumberOfBlobs()

	print, ' '
	print, '  -----------------------------------------------'
	print, ' '
	print, ' >>> Detected CP features above LOW threshold: ',float(count_l),format='(a,f6.0,$)' 
	print, ' '
	print, '  -----------------------------------------------'
	print, ' '


	; Report and save stats
	CPF_l -> ReportStats
	
	filestats = AnaDir+'PHOTOSPHERE/CP/stat__frame_'+strtrim(long(i),2)+'__nth_'+strtrim(long(nth1),2)+'.txt'
        CPF_l -> ReportStats, file=filestats, tofile=1, /noscale


   ;HIGH THRESHOLD CASE
   ;===================

count_h = CPF_h -> NumberOfBlobs()

	print, ' '
	print, '  -----------------------------------------------'
	print, ' '
	print, ' >>> Detected CP features above LOW threshold: ',float(count_h),format='(a,f6.0,$)' 
	print, ' '
	print, '  -----------------------------------------------'
	print, ' '


	; Report and save stats
        CPF_h -> ReportStats

	filestats = AnaDir+'PHOTOSPHERE/CP/stat__frame_'+strtrim(long(i),2)+'__nth_'+strtrim(long(nth2),2)+'.txt'
        CPF_h -> ReportStats, file=filestats, tofile=1, /noscale

    ;DISPLAYING the detected features
        
	; create and save a new frame that only includes the detected events: to be used for, e.g., tracking procedure
	blank_l = fltarr(nnx,nny)
        blank_h = fltarr(nnx,nny)
        
	s = Size(image_l, /DIMENSIONS)
	Window, 0, XSIZE=s[0]/cc, YSIZE=s[1]/cc, xpos=1925, ypos=1000
	nx=s[0]   &   ny=s[1]
	loadct, 0
        tvscl, histo_opt(congrid(CP, nx/cc, ny/cc, cubic=-0.5))

	; Display the CPF_ls we located with LABEL_REGION
        counto = 0
        
        ;OPLOTing the LOW threshold
        FOR j=0, count_l-1 DO BEGIN
           stats = CPF_l -> GetStats(j, /NoScale, XYINDICES=XYINDICE_l, INDICES=INDICE_l)
           if stats.count ge size_threshold then begin
              PLOTS, stats.perimeter_pts[0,*]/cc, stats.perimeter_pts[1,*]/cc, /Device, COLOR=cgColor('cyan'), thick=2
              XYOUTS, (stats.center[0]/cc)+5., (stats.center[1]/cc)-5, /Device, StrTrim(j,2), $
              COLOR=cgColor('navy'), ALIGNMENT=0., CHARSIZE=1.3, CHARTHICK=2
			
              ; save coordinates and indices of
              ; individual detected features
              outsc = fltarr(3,n_elements(INDICE_l))
              outnamecoords = FDir+'frame_'+strtrim(long(i),2)+'__nth_'+strtrim(long(nth1),2)+'_feature_'+strtrim(long(j),2)+'.txt'
              outsc(0,*)=XYINDICE_l(0,*)  &  outsc(1,*)=XYINDICE_l(1,*)  &  outsc(2,*)=INDICE_l
              write_table, outsc, outnamecoords
			
              blank_l[INDICE_l] = image_l[INDICE_l]
           endif
        ENDFOR

        ;OPLOTing the HIGH threshold
        FOR j=0, count_h-1 DO BEGIN
           stats = CPF_h -> GetStats(j, /NoScale, XYINDICES=XYINDICE_h, INDICES=INDICE_h)
           if stats.count ge size_threshold then begin
              PLOTS, stats.perimeter_pts[0,*]/cc, stats.perimeter_pts[1,*]/cc, /Device, COLOR=cgColor('red'), thick=2
              XYOUTS, (stats.center[0]/cc)+5., (stats.center[1]/cc)-5, /Device, StrTrim(j,2), $
              COLOR=cgColor('red7'), ALIGNMENT=0., CHARSIZE=1.3, CHARTHICK=2
			
              ; save coordinates and indices of
              ; individual detected features
              outsc = fltarr(3,n_elements(INDICE_h))
              outnamecoords = FDir+'frame_'+strtrim(long(i),2)+'__nth_'+strtrim(long(nth2),2)+'_feature_'+strtrim(long(j),2)+'.txt'
              outsc(0,*)=XYINDICE_h(0,*)  &  outsc(1,*)=XYINDICE_h(1,*)  &  outsc(2,*)=INDICE_h
              write_table, outsc, outnamecoords
			
              blank_h[INDICE_h] = image_h[INDICE_h]
           endif
        ENDFOR
   
saveimage, AnaDir+'PHOTOSPHERE/CP/CP__frame_'+strtrim(long(i),2)+'__features.png',quality =100,png=png
writefits, AnaDir+'PHOTOSPHERE/CP/CP__nth_'+strtrim(long(nth1),2)+'__frame_'+strtrim(long(i),2)+'__features.fits', blank_l
writefits, AnaDir+'PHOTOSPHERE/CP/CP__nth_'+strtrim(long(nth2),2)+'__frame_'+strtrim(long(i),2)+'__features.fits', blank_h


print, ' '
print, ' >>> frame '+strtrim(i,2)+' is analysed :)'
print, ' '
print, ' >>> type .c to continue with the next frame'
print, ' '

endif

;==============================================================================
;      CHROMOSPHERE       CHROMOSPHERE       CHROMOSPHERE       CHROMOSPHERE
;==============================================================================

endif

if fcb eq 10 then begin

print,'Ca I 8542 feature detection'
   
endif

wdelete
stop   
endfor

stop
end
