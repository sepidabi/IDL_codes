;save_frame.pro
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

if sm_fact eq 1 then begin
	restore, outputdir+strtrim(long(lcore),2)+'_candid__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'.sav'
	restore, datadir+strtrim(long(lcore),2)+'_polarization_maps.sav'
	restore,datadir+strtrim(long(8542),2)+'_polarization_maps.sav'
endif else begin
	restore, outputdir+strtrim(long(lcore),2)+'_candid__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'__smoothed.sav'
	restore, datadir+strtrim(long(lcore),2)+'_polarization_maps__smoothed.sav'
	restore,datadir+strtrim(long(8542),2)+'_polarization_maps__smoothed.sav'
endelse

n_candid = n_elements(candid)/2		;number of bipolar elements
nf = n_elements(features)			;total # of features

totLP = fltarr(nx,ny)	&	def = fltarr(nx,ny,nt)
inc = def	&	stren = def	&	vlos = def

cd, nadir+QSdir
fcube,f
for tt = 0L, 44 do begin
	totLP = totLP + LP(*,*,tt)
	inc(*,*,tt) = lp_get(f[2], tt)
	stren(*,*,tt) = lp_get(f[3], tt)
	vlos(*,*,tt) = lp_get(f[5], tt)
endfor
LPm = totLP/45.

for n_map = 0L, 5 do begin

		;plot the background map
		map = ''
		print, ''
		print, 'Set map to desired one ;) and press .c to continue...'
		print, '========================='
stop
		mapsav = ''
		print, 'Enter the maps name'		
		read, mapsav

	for fr = 0L, nt-1 do begin
		bgimg = map(*,*,fr)
		window,0,xsize = nx,ysize = ny
		plot_image,bgimg,/keep_aspect_ratio ,charsize = 1.5, title = 'FRAME '+string(fr)
stop
		WRITE_PNG, outputdir+'figures/frames/'+mapsav+'__fr_'+string(fr,format = '(i2.2)')+'.png', TVRD(True=1)
	endfor
endfor
stop
print,''
print,'Saving all parameters...'

save,candid,addresses,footpoint,nearby,features,masks,CP,sigmaCP,LP,Icont,CPch,sigmaCPch,vlos,inc,stren,$
	min_size,dx,nx,ny,nt,nthre,filename=outputdir+strtrim(long(lcore),2)+'_candid__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'__smoothed.sav'

stop
end

