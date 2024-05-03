;check_candid.pro
;	purpose = check the bipolars to find loop candids

; set some initial parameters
;====================================
read,lcore, prompt = 'Line Core Wavelength = '
read,sm_fact,prompt = 'Smoothing Factor = '
read, nthre, prompt = 'Signal threshold factor = '
read, min_size, prompt = 'Min size allowed (pixels) = '

nadir = '/nadir-scratch/jaime/CHROMIS/'
QSdir = '2016.10.12_quietsun_10_46_37_cakfine_crisp/'
datadir = '/home/seki2695/DATA/'
outputdir = '/home/seki2695/OUTPUT/YAFTA/'

if sm_fact eq 1 then begin
	restore, outputdir+strtrim(long(lcore),2)+'_footpoint__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'.sav'
	restore, datadir+strtrim(long(lcore),2)+'_polarization_maps.sav'
	restore,datadir+strtrim(long(8542),2)+'_polarization_maps.sav'
endif else begin
	restore, outputdir+strtrim(long(lcore),2)+'_footpoint__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'__smoothed.sav'
	restore, datadir+strtrim(long(lcore),2)+'_polarization_maps__smoothed.sav'
	restore,datadir+strtrim(long(8542),2)+'_polarization_maps__smoothed.sav'
endelse

n_footpoint = n_elements(footpoint)/2		;number of bipolar elements
nf = n_elements(features)			;total # of features

candid = intarr(2)

for l=0L,n_footpoint-1 do begin

	poslabel = footpoint(0,l)
	neglabel = footpoint(1,l)
	bipos1 = features(where(features.label eq poslabel,npos))	;positive side
	bineg1 = features(where(features.label eq neglabel,nneg))	;negative side

	if npos gt 0 then begin
	
		print,''
		print,'STOP!'
		print,''
		print,'footpoint index = '+strtrim(long(l),2)
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
		bimaski = masks(*,*,0)
		bimaski(*,*) = 0
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
		bimaskf = masks(*,*,0)
		bimaskf(*,*) = 0
		bimaskf(addposf) = strtrim(long(bipos1(npos-1).label),2)
		bimaskf(addnegf) = strtrim(long(bineg1(nneg-1).label),2)

		;==========================
		;bg -> desired frame number
		;========================== 
		;window,0,xsize = nx,ysize = ny;, title = 'CP MAP'
		;plot the background map with specified FOI
		;plot_image,bgimg,/keep_aspect_ratio,title = 'step 1',charsize = 1.5,/silent

		;plot_edges,bimaski, setcolor = 210,thick = 0.25

		;plot_edges,bimaskf, setcolor = 130,thick = 0.25

		;window,0,xsize = nx,ysize = ny;, title = 'CP MAP'
		;plot the background map with specified FOI
		;plot_image,bgimg,/keep_aspect_ratio,title = 'step 1',charsize = 1.5,/silent
		;masking the (+) polarity feature
		;str_pos = (features(where(features.label eq poslabel)))(bg).mask_str
		;addpos = long(strsplit(str_pos,/extract))
		;masking the (-) polarity feature
		;str_neg = (features(where(features.label eq neglabel)))(bg).mask_str
		;addneg = long(strsplit(str_neg,/extract))
		
		;BI-POLAR element map:
		;bimask = masks(*,*,0)
		;bimask(*,*) = 0
		;bimask(addpos) = strtrim(long(bipos1(0).label),2)
		;bimask(addneg) = strtrim(long(bineg1(0).label),2)
		;plot_edges,bimask, setcolor = 300,thick = 0.25
			
		;feature movement direction		
		cgarrow, xnegi,ynegi,xnegf,ynegf,/data,/solid,/overplot, color = 'dodger blue', hsize = 10.
		cgarrow, xposi,yposi,xposf,yposf,/data,/solid,/overplot, color = 'red', hsize = 10.

;WRITE_PNG, outputdir+'figures/frames/'+'indx_'+strtrim(long(l),2)+'.png', TVRD(True=1)

		;EVOLUTION PLOTS
		;===============
		;TRAJECTORY plot					
		window,1, title = 'TRAJECTORY'			
		xmin = min([bipos1.x,bineg1.x])	&	xmax = max([bipos1.x,bineg1.x])
		ymin = min([bipos1.y,bineg1.y])	&	ymax = max([bipos1.y,bineg1.y])
		cgplot,bipos1.x,bipos1.y, ytitle = '[px]', xtitle = '[px]', $
		xrange = [xmin,xmax], yrange = [ymin,ymax], color = 'red'
		cgplot,bineg1.x,bineg1.y, ytitle = '[px]', xtitle = '[px]',color = 'blue',/overplot
		cgarrow, xnegi,ynegi,xnegf,ynegf,/data,/solid,/overplot, color = 'blue', linestyle = 2
		cgarrow, xposi,yposi,xposf,yposf,/data,/solid,/overplot, color = 'red', linestyle = 2
WRITE_PNG, outputdir+'figures/frames/'+'indx_'+strtrim(long(l),2)+'_trajectory.png', TVRD(True=1)

		;SIZE & SIGNAL
		window,2, title = 'EVOLUTION'
		cgplot,bipos1.step,(bipos1.phi/bipos1.size)*100, ytitle = 'Avg. Abs. CP (% I$\downc$)', xtitle = 'Steps', yrange = [0,max([bipos1.phi/bipos1.size,bineg1.phi/bineg1.size])*100], xrange = [0,nt],color = 'red', YStyle=8, thick = 1.5
		cgplot,bineg1.step,(bineg1.phi/bineg1.size)*100,/overplot,color = 'blue', thick = 1.5

		cgAxis, YAxis=1, YRange=[0,max([bipos1.size,bineg1.size])], /Save, ytitle = 'Size (px)'
		cgOplot,bipos1.step,bipos1.size,color = 'red', linestyle = 2, thick = 1.5
		cgOplot,bineg1.step,bineg1.size,/overplot,color = 'blue', linestyle = 2, thick = 1.5
		cgLegend, Title=['CP', 'size'], LineStyle=[0,2], Color=['black'], Location=[0.8,0.96]
		WRITE_PNG, outputdir+'figures/frames/'+'indx_'+strtrim(long(l),2)+'_evolution.png', TVRD(True=1)

		;storing the accepted footpoint candidates
		;=========================================
		dec = ''		
		read,dec, prompt = 'Accept the candid (y/n)? '

		if dec eq 'y' then begin

			if n_elements(candid) eq 0 then candid = [(bipos1.label)(0),(bineg1.label)(0)] $
			else candid = [[candid],[(bipos1.label)(0),(bineg1.label)(0)]]

		endif
		wdelete, 0, 1, 2
	endif
endfor
stop
print,''
print,'Saving all parameters...'

save,candid,addresses,footpoint,nearby,features,masks,CP,sigmaCP,LP,Icont,CPch,sigmaCPch,$
	min_size,dx,nx,ny,nt,nthre,min_size,filename=outputdir+strtrim(long(lcore),2)+'_candid__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'__smoothed.sav'

stop
end

