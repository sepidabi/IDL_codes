pro cube_make

close, /all

; -------- To be modified:
datadir = '/home/sepid/SunriseData/NR/'
savedir = '/home/sepid/SUNRISE_DATA_OUTPUTS/DATA_CUBES/'

sm_fact = 3.		; smooth box car dimension (odd number; sm_fact=1 means no smoothing)
nthre = 3.	 	; signal threshold in noise level (not necessarily, but usually larger than 4 -- to be decided)
size_threshold = 10 	; size threshold in pixel (if you want to have it; I do recommend that -- 10 makes you sure that spatially 				unresolved features are not included)
cadnc = 33 		; cadence of observations in seconds
res = 45./936.

cc = 1 			; a compresion factor for window's size
;---------------------------------------------------------------

; Cube making
CPube_name = savedir+'CP_NR_smooth'+strtrim(long(sm_fact),2)+'.fcube'		;  name of the output cube (in assoc float format)
LPube_name = savedir+'LP_NR_smooth'+strtrim(long(sm_fact),2)+'.fcube'		;  name of the output cube (in assoc float format)
ICube_name = savedir+'IC_NR_smooth'+strtrim(long(sm_fact),2)+'.fcube'		;  name of the output cube (in assoc float format)

nx= 936. & ny= 936. ; size of images
nt = 58 ; number of frames in the cube

header=make_lp_header(fltarr(nx,ny),nt=nt)
openw, luwCP, CPube_name, /get_lun
openw, luwLP, LPube_name, /get_lun
openw, luwIC, ICube_name, /get_lun

;defining the cube 'CPube'
CPube=assoc(luwCP, bytarr(512)) & CPube[0]=byte(header)
CPube=assoc(luwCP,fltarr(nx,ny),512)

;defining the cube 'LPube'
LPube=assoc(luwLP, bytarr(512)) & LPube[0]=byte(header)
LPube=assoc(luwLP,fltarr(nx,ny),512)

;defining the cube 'ICube'
ICube=assoc(luwIC, bytarr(512)) & ICube[0]=byte(header)
ICube=assoc(luwIC,fltarr(nx,ny),512)

;-----------------------------------------------------------------

; read data:
	files = file_search(datadir+'imax_nr_163_*.fits')
	nim = n_elements(files) ; number of files/images

	print, ' '
	print, ' ===================================================='
	print, '  begin to (1)read NR data & (2)set the data in cube'
	print, ' ===================================================='
	print, ' '


    for t=0L, nim-1 do begin

	; Stokes
		
	imcube = readfits(files[t],/silent)
	
	StokesI = reform(imcube[*,*,*,0])
	StokesQ = reform(imcube[*,*,*,1])
	StokesU = reform(imcube[*,*,*,2])
	StokesV = reform(imcube[*,*,*,3])
	
	nx = (size(StokesI))[1]  &  ny = (size(StokesI))[2] ; dimension of images
	
	; Making Smoothed Stokes dataset

	I1 = smooth(reform(StokesI[*,*,0]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	I2 = smooth(reform(StokesI[*,*,1]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	I3 = smooth(reform(StokesI[*,*,2]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	I4 = smooth(reform(StokesI[*,*,3]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
     Icont = smooth(reform(StokesI[*,*,4]),[sm_fact,sm_fact],/EDGE_TRUNCATE)	; Stokes I continuum

	Q1 = smooth(reform(StokesQ[*,*,0]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	Q2 = smooth(reform(StokesQ[*,*,1]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	Q3 = smooth(reform(StokesQ[*,*,2]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	Q4 = smooth(reform(StokesQ[*,*,3]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
     Qcont = smooth(reform(StokesQ[*,*,4]),[sm_fact,sm_fact],/EDGE_TRUNCATE)	; Stokes Q continuum

	U1 = smooth(reform(StokesU[*,*,0]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	U2 = smooth(reform(StokesU[*,*,1]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	U3 = smooth(reform(StokesU[*,*,2]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	U4 = smooth(reform(StokesU[*,*,3]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
     Ucont = smooth(reform(StokesU[*,*,4]),[sm_fact,sm_fact],/EDGE_TRUNCATE)	; Stokes U continuum

	V1 = smooth(reform(StokesV[*,*,0]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	V2 = smooth(reform(StokesV[*,*,1]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	V3 = smooth(reform(StokesV[*,*,2]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	V4 = smooth(reform(StokesV[*,*,3]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
     Vcont = smooth(reform(StokesV[*,*,4]),[sm_fact,sm_fact],/EDGE_TRUNCATE)	; Stokes V continuum

	; Averaged Stokes over 4 wavelength positions
	Imean = float(I1+I2+I3+I4)/4.

	Qmean = float(Q1+Q2+Q3+Q4)/4.
	
	Umean = float(U1+U2+U3+U4)/4.	

	; total circular polarisation:
	CP = float(V1+V2-V3-V4)/(4.*float(Icont))
    sigmaV = stddev(float(Vcont)/float(Icont)) ; Stokes V noise level, normalised to local Stokes I continuum
   CPsigma = sigmaV/2.

	; total linear polarisation:
	LP = sqrt(float(Qmean)^2+float(Umean)^2)/float(Icont)	;low-noise
	
    LPcont = sqrt(float(Qcont)^2+float(Ucont)^2)/float(Icont)	; LP continuum
   LPsigma = stddev(float(LPcont))/2.

;=================================================================
	CPube[t] = CP
	LPube[t] = LP
	ICube[t] = Icont

	print,'Frame No. '+strtrim(t+1,2)+' finished.'

endfor
print, '================================================'
print, 'Finished making the NR cubes'
print, 'cube saved in "'+ savedir+'"'
print, '================================================'
print, ' '

stop
free_lun, luwCP
free_lun, luwLP
free_lun, luwIC


;=========================================================================================================================================
;=========================================================================================================================================

; -------- To be modified:
datadir = '/home/sepid/SunriseData/R/'
;---------------------------------------------------------------

; Cube making
CPube_name = savedir+'CP_smooth'+strtrim(long(sm_fact),2)+'.fcube'		;  name of the output cube (in assoc float format)
LPube_name = savedir+'LP_smooth'+strtrim(long(sm_fact),2)+'.fcube'		;  name of the output cube (in assoc float format)
ICube_name = savedir+'IC_smooth'+strtrim(long(sm_fact),2)+'.fcube'		;  name of the output cube (in assoc float format)

nx= 936. & ny= 936. ; size of images
nt = 58 ; number of frames in the cube

header=make_lp_header(fltarr(nx,ny),nt=nt)
openw, luwCP, CPube_name, /get_lun
openw, luwLP, LPube_name, /get_lun
openw, luwIC, ICube_name, /get_lun

;defining the cube 'CPube'
CPube=assoc(luwCP, bytarr(512)) & CPube[0]=byte(header)
CPube=assoc(luwCP,fltarr(nx,ny),512)

;defining the cube 'LPube'
LPube=assoc(luwLP, bytarr(512)) & LPube[0]=byte(header)
LPube=assoc(luwLP,fltarr(nx,ny),512)

;defining the cube 'ICube'
ICube=assoc(luwIC, bytarr(512)) & ICube[0]=byte(header)
ICube=assoc(luwIC,fltarr(nx,ny),512)

;-----------------------------------------------------------------

; read data:
	files = file_search(datadir+'imax_163_*.fits')
	nim = n_elements(files) ; number of files/images

	print, ' '
	print, ' ================================================================='
	print, '  begin to (1)read NR data & (2)set the data in reconsructed cubes'
	print, ' ================================================================='
	print, ' '


    for t=0L, nim-1 do begin

	; Stokes
		
	imcube = readfits(files[t],/silent)
	
	StokesI = reform(imcube[*,*,*,0])
	StokesQ = reform(imcube[*,*,*,1])
	StokesU = reform(imcube[*,*,*,2])
	StokesV = reform(imcube[*,*,*,3])
	
	nx = (size(StokesI))[1]  &  ny = (size(StokesI))[2] ; dimension of images
	
	; Making Smoothed Stokes dataset

	I1 = smooth(reform(StokesI[*,*,0]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	I2 = smooth(reform(StokesI[*,*,1]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	I3 = smooth(reform(StokesI[*,*,2]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	I4 = smooth(reform(StokesI[*,*,3]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
     Icont = smooth(reform(StokesI[*,*,4]),[sm_fact,sm_fact],/EDGE_TRUNCATE)	; Stokes I continuum

	Q1 = smooth(reform(StokesQ[*,*,0]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	Q2 = smooth(reform(StokesQ[*,*,1]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	Q3 = smooth(reform(StokesQ[*,*,2]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	Q4 = smooth(reform(StokesQ[*,*,3]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
     Qcont = smooth(reform(StokesQ[*,*,4]),[sm_fact,sm_fact],/EDGE_TRUNCATE)	; Stokes Q continuum

	U1 = smooth(reform(StokesU[*,*,0]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	U2 = smooth(reform(StokesU[*,*,1]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	U3 = smooth(reform(StokesU[*,*,2]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	U4 = smooth(reform(StokesU[*,*,3]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
     Ucont = smooth(reform(StokesU[*,*,4]),[sm_fact,sm_fact],/EDGE_TRUNCATE)	; Stokes U continuum

	V1 = smooth(reform(StokesV[*,*,0]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	V2 = smooth(reform(StokesV[*,*,1]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	V3 = smooth(reform(StokesV[*,*,2]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
	V4 = smooth(reform(StokesV[*,*,3]),[sm_fact,sm_fact],/EDGE_TRUNCATE)
     Vcont = smooth(reform(StokesV[*,*,4]),[sm_fact,sm_fact],/EDGE_TRUNCATE)	; Stokes V continuum

	; Averaged Stokes over 4 wavelength positions
	Imean = float(I1+I2+I3+I4)/4.

	Qmean = float(Q1+Q2+Q3+Q4)/4.
	
	Umean = float(U1+U2+U3+U4)/4.	

	; total circular polarisation:
	CP = float(V1+V2-V3-V4)/(4.*float(Icont))
    sigmaV = stddev(float(Vcont)/float(Icont)) ; Stokes V noise level, normalised to local Stokes I continuum
   CPsigma = sigmaV/2.

	; total linear polarisation:
	LP = sqrt(float(Qmean)^2+float(Umean)^2)/float(Icont)	;low-noise
	
    LPcont = sqrt(float(Qcont)^2+float(Ucont)^2)/float(Icont)	; LP continuum
   LPsigma = stddev(float(LPcont))/2.

;=================================================================
	CPube[t] = CP
	LPube[t] = LP
	ICube[t] = Icont

	print,'Frame No. '+strtrim(t+1,2)+' finished.'

endfor
print, '================================================'
print, 'Finished making the reconstructed cubes'
print, 'cube saved in "'+ savedir+'"'
print, '================================================'
print, ' '

stop
free_lun, luwCP
free_lun, luwLP
free_lun, luwIC

end
