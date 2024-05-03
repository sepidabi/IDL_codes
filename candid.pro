;candid.pro
;	purpose = check the bipolars to find loop candids

; set some initial parameters
;====================================
;read, lcore, prompt = 'Line Core Wavelength = '
lcore = 6302
;read,sm_fact,prompt = 'Smoothing Factor = '
sm_fact = 3
;read, nthre, prompt = 'Signal threshold factor = '
nthre = 4
;read, min_size, prompt = 'Min size allowed (pixels) = '
min_size = 10
res = 0.0375

nadir = '/nadir-scratch/jaime/CHROMIS/'
QSdir = '2016.10.12_quietsun_10_46_37_cakfine_crisp/'
datadir = '/nadir-scratch/sepid/DATA/'
outputdir = '/nadir-scratch/sepid/YAFTA/'
figuredir = '/home/seki2695/OUTPUT/figures/'

print,'=========================='
print,'  RESTORING PARAMETERS...'
PRINT,'=========================='

;if sm_fact eq 1 then begin
;	restore, outputdir+strtrim(long(lcore),2)+'_bipolar__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'.sav'
;	restore, datadir+strtrim(long(lcore),2)+'_polarization_maps.sav'
;	restore,datadir+strtrim(long(8542),2)+'_polarization_maps.sav'
;endif else begin
	restore, datadir+strtrim(long(lcore),2)+'_polarization_maps__smoothed.sav'
	restore,datadir+strtrim(long(8542),2)+'_polarization_maps__smoothed.sav'
	restore,datadir+strtrim(long(3934),2)+'_intensity_maps.sav'
	restore,outputdir+strtrim(long(8542),2)+'_YAFTA_output__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'__smoothed.sav'
	restore,outputdir+strtrim(long(lcore),2)+'_bipolar__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'__smoothed.sav'
	restore, outputdir+'YAFTA_output_Bh__4sigma__size12__smoothed.sav'
	restore, outputdir+'6302_YAFTA_output_LP__7sigma__size12__smoothed.sav'
;endelse

n_bipolar = n_elements(bipolar)/2		;number of bipolar elements
nf = n_elements(features)			;total # of features

totLP = fltarr(nx,ny)
def = fltarr(nx,ny,nt)
inc = def
stren = def
vlos = def
LPmask = def
Bhmask = def
CPmask = def
CPchmask = def
Bh = def
maskI = def
maskF = def

candid = intarr(2)
tag = intarr(2)

cd, nadir+QSdir
fcube,f

;=============================
; Gefring all the maps needed
;=============================
for fr = 0L, 44 do begin
	totLP = totLP + LP(*,*,fr)
	Bh(*,*,fr) = lp_get(f[2], fr)
	inc(*,*,fr) = lp_get(f[2], fr)
	stren(*,*,fr) = lp_get(f[3], fr)
	vlos(*,*,fr) = lp_get(f[5], fr)
;stop
	;masking LP by [0,1]
	tmask = masksLP(*,*,fr)
	tmask(*,*) = 0
	LPup = where(masksLP(*,*,fr) gt 0)	
	tmask(LPup) = 1
	LPmask(*,*,fr) = tmask

	;masking Bh by [0,1]
	bmask = masksBh(*,*,fr)
	bmask(*,*) = 0
	Bhup = where(masksBh(*,*,fr) gt 0)	
	bmask(Bhup) = 1
	Bhmask(*,*,fr) = Bmask

	;masking CP-6302 by [-1,0,1]
	cmask = masks(*,*,fr)
	cmask(*,*) = 0
	cup = where(masks(*,*,fr) gt 0)	
	cmask(cup) = 1
	clow = where(masks(*,*,fr) lt 0)	
	cmask(clow) = -1
	CPmask(*,*,fr) = cmask

	;masking CP-8542 by [-1,0,1]
	chmask = masksch(*,*,fr)
	chmask(*,*) = 0
	chup = where(masksch(*,*,fr) gt 0)	
	chmask(chup) = 1
	chlow = where(masksch(*,*,fr) lt 0)	
	chmask(chlow) = -1
	CPchmask(*,*,fr) = chmask
endfor
	LPm = totLP/45.

;stop
;SAVE GRAPHICS
;      thisDevice = !D.Name
;      Set_Plot, "PS"
;      Device, filename = figuredir+'CP'+strtrim(long(fr),2)+'.eps'
;      cgimage,CP(*,*,0),xrange = [0,nx-1]*res,yrange = [0,ny-1]*res,/axis,/scale,/keep_aspect_ratio,$
;axkeywords = {xticklen:-0.02,yticklen:-0.02},xtitle = '[arcsec]',ytitle = '[arcsec]', charsize = 1.
;      ;cgcontour,Bhmask(*,*,0), /onimage,/noerase,color = 'cyan'
;      Device, /Close_File
;      Set_Plot, thisDevice

;========================================================
;    Extracting information from each pair of bipolars
;========================================================
for l=0L,n_bipolar-1 do begin

	bimaski = masks(*,*,0)
	bimaskf = masks(*,*,0)
	bimaski(*,*) = 0
	bimaskf(*,*) = 0

	poslabel = bipolar(0,l)
	neglabel = bipolar(1,l)
	bipos1 = features(where(features.label eq poslabel,npos))	;positive side
	bineg1 = features(where(features.label eq neglabel,nneg))	;negative side

  if npos gt 0 then begin
	
		print,''
		print,'STOP!'
		print,''
		print,'candid index = '+strtrim(long(l),2)
		print,''
		print,'Signal Label (+) = '+strtrim(long(poslabel),2)
		print,'Signal Label (-) = '+strtrim(long(neglabel),2)

		xnegi = (bineg1.x)(0) & xnegf = (bineg1.x)(nneg-1)
		ynegi = (bineg1.y)(0) & ynegf = (bineg1.y)(nneg-1)
		xposi = (bipos1.x)(0) & xposf = (bipos1.x)(npos-1)
		yposi = (bipos1.y)(0) & yposf = (bipos1.y)(npos-1)

		if n_elements(tag) eq 0 then tag = [(xnegi+xposi)/2.,(ynegi+yposi)/2.] $
		else tag = [[tag],[(xnegi+xposi)/2.,(ynegi+yposi)/2.]]

		;extracting features info in 1ST FRAME
		;=====================================
		;masking the (+) polarity feature
		str_posi = (features(where(features.label eq poslabel)))(0).mask_str
		addposi = long(strsplit(str_posi,/extract))
		;masking the (-) polarity feature
		str_negi = (features(where(features.label eq neglabel)))(0).mask_str
		addnegi = long(strsplit(str_negi,/extract))

		;BI-POLAR element map:
		bimaski(addposi) = strtrim(long(bipos1(0).label),2)
		bimaski(addnegi) = strtrim(long(bineg1(0).label),2)

		;extracting features info in last FRAME
		;=====================================
		;masking the (+) polarity feature
		str_posf = (features(where(features.label eq poslabel)))(npos-1).mask_str
		addposf = long(strsplit(str_posf,/extract))

		;masking the (-) polarity feature
		str_negf = (features(where(features.label eq neglabel)))(nneg-1).mask_str
		addnegf = long(strsplit(str_negf,/extract))
			
		;BI-POLAR element map:
		bimaskf(addposf) = strtrim(long(bipos1(npos-1).label),2)
		bimaskf(addnegf) = strtrim(long(bineg1(nneg-1).label),2)

	;===============================================================
	; After Checking the movies decide to keep which pair to candid
	;===============================================================
	window,0,xsize = nx,ysize = ny
	bgimg = CPmask(*,*,0)
	plot_image,bgimg,/keep_aspect_ratio ,charsize = 1.5;, title = '(+) = '+strtrim(long(poslabel),2)+'   ';+'(-) = '+strtrim(long(neglabel),2)
	;plot candids total mask
	plot_edges, bimaski, setcolor = 210,thick = 0.25
	plot_edges, bimaskf, setcolor = 130,thick = 0.25

	;decide whether to keep the bipolar
	decide = dialog_message('Candid?',/question)

	if decide eq 'Yes' then begin
		if n_elements(candid) eq 0 then candid = [(bipos1.label)(0),(bineg1.label)(0)] $
		else candid = [[candid],[(bipos1.label)(0),(bineg1.label)(0)]]
	endif

  endif

endfor
;++++++++++++++++++++++++++++++++++++++++++
;	Movies must be prepared here
;++++++++++++++++++++++++++++++++++++++++++
print,'CREATE THE MOVIES!!!'
stop

print,''
print,'Saving all parameters...'

save,candid,addresses,nearby,features,masks,CP,sigmaCP,LP,Icont,CPch,sigmaCPch,$
	min_size,dx,nx,ny,nt,nthre,filename=outputdir+strtrim(long(lcore),2)+'_candid__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'__smoothed.sav'

stop
end

