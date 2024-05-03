;mask_candid.pro
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
min_size = 12

nadir = '/nadir-scratch/jaime/CHROMIS/'
QSdir = '2016.10.12_quietsun_10_46_37_cakfine_crisp/'
datadir = '/home/seki2695/DATA/'
outputdir = '/home/seki2695/OUTPUT/YAFTA/'

print,'=========================='
print,'  RESTORING PARAMETERS...'
PRINT,'=========================='

if sm_fact eq 1 then begin
	restore, outputdir+strtrim(long(lcore),2)+'_candid__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'.sav'
	restore, datadir+strtrim(long(lcore),2)+'_polarization_maps.sav'
	restore, datadir+strtrim(long(8542),2)+'_polarization_maps.sav'
endif else begin
	restore, outputdir+strtrim(long(lcore),2)+'_candid__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'__smoothed.sav'
	restore, datadir+strtrim(long(lcore),2)+'_polarization_maps__smoothed.sav'
	restore, datadir+strtrim(long(8542),2)+'_polarization_maps__smoothed.sav'
	restore, outputdir+'YAFTA_output_Bh__4sigma__size12__smoothed.sav'
	restore, outputdir+'6302_YAFTA_output_LP__7sigma__size12__smoothed.sav'
endelse

n_candid = n_elements(candid)/2		;number of bipolar elements
nf = n_elements(features)			;total # of features

totLP = fltarr(nx,ny)
def = fltarr(nx,ny,nt)
inc = def
stren = def
vlos = def
Bh = def
cont400 = def
;LPmask = def
maskI = def
maskF = def

cd, nadir+QSdir
fcube,f

for fr = 0L, 44 do begin
	totLP = totLP + LP(*,*,fr)
	Bh(*,*,fr) = lp_get(f[1], fr)
	inc(*,*,fr) = lp_get(f[2], fr)
	stren(*,*,fr) = lp_get(f[3], fr)
	vlos(*,*,fr) = lp_get(f[5], fr)
	cont400(*,*,fr) = congrid(lp_get(f[12], fr*36+35),nx,ny)

	;masking LP maps
	;nthreLP = 7.
	;below = LP(*,*,fr)
	;below_add = where(LP(*,*,fr) lt nthreLP*sigmaLP(fr))
	;below(below_add) = 0
	;LPmask(*,*,fr) = below

	;candid mask (initial state)
	bimaski = masks(*,*,fr)
	bimaski(*,*) = 0

	;candid mask (final state)
	bimaskf = masks(*,*,0)
	bimaskf(*,*) = 0
stop
   for l=0L,n_candid-1 do begin

	poslabel = candid(0,l)
	neglabel = candid(1,l)
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

	endif
   endfor
;stop
  ;total masks BI-CPFs candids:
  maskI(*,*,fr) = bimaski
  maskF(*,*,fr) = bimaskf

  map = ''
  mp = 10  
  window,0,xsize = nx,ysize = ny
  
  ;for mp = 0L , 9 do begin
;	if mp eq 0 then map = CP
;	if mp eq 1 then map = CPch
;	if mp eq 2 then map = LP
;	if mp eq 3 then map = LPmask
;	if mp eq 4 then map = Icont
;	if mp eq 5 then map = stren
;	if mp eq 6 then map = inc
;	if mp eq 7 then map = vlos
;	if mp eq 8 then map = Icontch
	if mp eq 9 then map = Bh
	if mp eq 10 then map = cont400

	;plot the background map
	bgimg = map(*,*,fr)
	plot_image,bgimg,/keep_aspect_ratio ,charsize = 1.5;, title = '(+) = '+strtrim(long(poslabel),2)+'   ';+'(-) = '+strtrim(long(neglabel),2)
	;plot candids total mask
	plot_edges, maskI(*,*,fr), setcolor = 210,thick = 0.25
	plot_edges, maskF(*,*,fr), setcolor = 130,thick = 0.25
	stop
	WRITE_PNG, outputdir+'figures/frames/'+'fr_'+string(fr,format = '(i2.2)')+'__map_'+strtrim(long(mp),2)+'.png', TVRD(True=1)
  ;endfor
stop		
endfor
wdelete, 0, 1, 2
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;						EVOLUTION PLOTS
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		;===============
		;TRAJECTORY plot					
		window,1, title = 'TRAJECTORY'			
		xmin = min([bipos1.x,bineg1.x])	&	xmax = max([bipos1.x,bineg1.x])
		ymin = min([bipos1.y,bineg1.y])	&	ymax = max([bipos1.y,bineg1.y])
		cgplot,bipos1.x,bipos1.y, ytitle = '[px]', xtitle = '[px]', $
		xrange = [xmin,xmax], yrange = [ymin,ymax], color = 'red', title = '(+) = '+strtrim(long(poslabel),2)+'   '+'(-) = '+strtrim(long(neglabel),2)
		cgplot,bineg1.x,bineg1.y, ytitle = '[px]', xtitle = '[px]',color = 'blue',/overplot
		cgarrow, xnegi,ynegi,xnegf,ynegf,/data,/solid,/overplot, color = 'blue', linestyle = 2
		cgarrow, xposi,yposi,xposf,yposf,/data,/solid,/overplot, color = 'red', linestyle = 2
WRITE_PNG, outputdir+'figures/frames/'+'candid_'+strtrim(long(l),2)+'_trajectory.png', TVRD(True=1)

		;SIZE & SIGNAL
		window,2, title = 'EVOLUTION'
		cgplot,bipos1.step,(bipos1.phi/bipos1.size)*100, ytitle = 'Avg. Abs. CP (% I$\downc$)', xtitle = 'Steps', yrange = [0,max([bipos1.phi/bipos1.size,bineg1.phi/bineg1.size])*100], xrange = [0,nt],color = 'red', YStyle=8, thick = 1.5, title = '(+) = '+strtrim(long(poslabel),2)+'   '+'(-) = '+strtrim(long(neglabel),2)
		cgplot,bineg1.step,(bineg1.phi/bineg1.size)*100,/overplot,color = 'blue', thick = 1.5

		cgAxis, YAxis=1, YRange=[0,max([bipos1.size,bineg1.size])], /Save, ytitle = 'Size (px)'
		cgOplot,bipos1.step,bipos1.size,color = 'red', linestyle = 2, thick = 1.5
		cgOplot,bineg1.step,bineg1.size,/overplot,color = 'blue', linestyle = 2, thick = 1.5
		cgLegend, Title=['CP', 'size'], LineStyle=[0,2], Color=['black'], Location=[0.8,0.96]
		WRITE_PNG, outputdir+'figures/frames/'+'candid_'+strtrim(long(l),2)+'_evolution.png', TVRD(True=1)

;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


STOP

print,''
print,'Saving all parameters...'

save,candid,addresses,footpoint,nearby,features,masks,CP,sigmaCP,LP,Icont,CPch,sigmaCPch,vlos,inc,stren,$
	min_size,dx,nx,ny,nt,nthre,filename=outputdir+strtrim(long(lcore),2)+'_candid__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'__smoothed.sav'

stop
end

