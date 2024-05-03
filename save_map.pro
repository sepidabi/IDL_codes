;save_map.pro
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
figuredir = '/home/seki2695/OUTPUT/figures/sample_seq/'

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
	restore, outputdir+strtrim(long(8542),2)+'_YAFTA_output__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'__smoothed.sav'
	restore, outputdir+'6302_YAFTA_output_LP__7sigma__size12__smoothed.sav'
	restore, outputdir+'YAFTA_output_Bh__4sigma__size12__smoothed.sav'
	restore, outputdir+'6302_YAFTA_output__4sigma__size12__smoothed.sav'
;endelse

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

cd, nadir+QSdir
fcube,f

xpar1 = 600
xpar2 = 700
ypar1 = 120
ypar2 = 200

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

;SAVE GRAPHICS

  ;for mp = 0L,27 do begin
mp = 25

      ;SET THE BACKGROUND MAP
	MAPS = ''
	if mp eq 0 then begin
	 map = 'CPmask' & title = 'Fe I 6301 CP mask' & MAPS = CPmask
	endif
	if mp eq 1 then begin
	 map = 'CPch' & title = 'Ca II 8542 CP' & MAPS = CPch
	endif
	if mp eq 2 then begin
	 map = 'Bhmask' & title = 'Fe I B$\downh$' & MAPS = Bhmask
	endif
	if mp eq 3 then begin
	 map = 'LP' & title = 'Fe I 6301 LP' & MAPS = LP
	endif
	if mp eq 4 then begin
	 map = 'Icont' & title = 'Continuum (630.1 + 0.07 nm)' & MAPS = Icont
	endif
	if mp eq 5 then begin
	 map = 'Icontch' & title = 'Continuum (854.2 + 0.17)' & MAPS = Icontch
	endif
	if mp eq 6 then begin
	 map = 'I854' & title = 'Ca II Core (854.20 nm)' & MAPS = I854
	endif
	if mp eq 7 then begin
	 map = 'I854_1' & title = 'Ca II Wing (854.20 - 0.08 nm)' & MAPS = I854_1
	endif
	if mp eq 8 then begin
	 map = 'I854_4' & title = 'Ca II Wing (854.20 - 0.05 nm)' & MAPS = I854_4
	endif
	if mp eq 9 then begin
	 map = 'I854_7' & title = 'Ca II Core (854.20 - 0.02 nm)' & MAPS = I854_7
	endif
	if mp eq 10 then begin
	 map ='I854_13' & title = 'Ca II Core (854.20 + 0.02 nm)' & MAPS = I854_13
	endif
	if mp eq 11 then begin
	 map = 'I854_16' & title = 'Ca II Wing (854.20 + 0.05 nm)' & MAPS = I854_16
	endif
	if mp eq 12 then begin
	 map = 'I854_19' & title = 'Ca II Wing (854.20 - 0.08 nm)' & MAPS = I854_19
	endif
	if mp eq 13 then begin
	 map = 'Ik' & title = 'Ca II K3 (393.40 nm)' & MAPS = Ik
	endif
	if mp eq 14 then begin
	 map = 'Ik_0' & title = 'Ca II K Wing (393.40 - 0.14 nm)' & MAPS = Ik_0
	endif
	if mp eq 15 then begin
	 map = 'Ik_3' & title = 'Ca II K Wing (393.40 - 0.10 nm)' & MAPS = Ik_3
	endif
	if mp eq 16 then begin
	 map = 'Ik_7' & title = 'Ca II K Wing (393.40 - 0.06 nm)' & MAPS = Ik_7
	endif
	if mp eq 17 then begin
	 map = 'Ik_11' & title = 'Ca II K1V (393.40 - 0.04 nm)' & MAPS = Ik_11
	endif
	if mp eq 18 then begin
	 map = 'Ik_13' & title = 'Ca II K2V (393.40 - 0.02 nm)' & MAPS = Ik_13
	endif
	if mp eq 19 then begin
	 map = 'Ik_14' & title = 'Ca II K2V (393.40 - 0.02 nm)' & MAPS = Ik_14
	endif
	if mp eq 20 then begin
	 map = 'Ik_16' & title = 'Ca II Core (393.40 - 0.00 nm)' & MAPS = Ik_16
	endif
	if mp eq 21 then begin
	 map = 'Ik_19' & title = 'Ca II K2R (393.40 + 0.01 nm)' & MAPS = Ik_19
	endif
	if mp eq 22 then begin
	 map = 'Ik_21' & title = 'Ca II K1R (393.40 + 0.02 nm)' & MAPS = Ik_21
	endif
	if mp eq 23 then begin
	 map = 'Ik_23' & title = 'Ca II K Wing (393.40 + 0.03 nm)' & MAPS = Ik_23
	endif
	if mp eq 24 then begin
	 map = 'Ik_27' & title = 'Ca II K Wing (393.40 + 0.06 nm)' & MAPS = Ik_27
	endif
	if mp eq 25 then begin
	 map = 'Ik_32' & title = 'Ca II K Wing (393.40 + 0.12 nm)' & MAPS = Ik_32
	endif
	if mp eq 26 then begin
	 map = 'Ik_34' & title = 'Ca II K Wing (393.40 + 0.14 nm)' & MAPS = Ik_34
	endif
	if mp eq 27 then begin
	 map = 'Icont400' & title = 'Continuum (400.00 nm)' & MAPS = Icont400
	endif

	cgwindow,wxsize = 400,wysize = 650, wmulti = [0, 5, 9], woxmargin = [0,5], woymargin = [0,10]

		for fr=0L,44 DO BEGIN
			bgimg = MAPS(xpar1:xpar2,ypar1:ypar2,fr)
			cimg1 = CPmask(xpar1:xpar2,ypar1:ypar2,fr)
			cimg2 = Bhmask(xpar1:xpar2,ypar1:ypar2,fr)

			axformat = {XTICKFORMAT:"(A1)",YTICKFORMAT:"(A1)",xticklen:-0.025,yticklen:-0.025}
			cgimage,/addcmd,bgimg,/scale,/keep_aspect_ratio, MultiMargin = 0.25,/axis,axkeywords = axformat
			cgcontour,cimg1, /onimage,/noerase,color = 'blue',levels = [0], /addcmd, c_charsize = 0.01, c_labels = ['']
			cgcontour,cimg1, /onimage,/noerase,color = 'red',levels = [1.], /addcmd, c_charsize = 0.01, c_labels = ['']
			cgcontour,cimg2, /onimage,/noerase,color = 'cyan', /addcmd, c_charsize = 0.01, c_labels = ['']
		endfor

	;!P.Multi = 0
	cgColorbar, Range=[min(MAPS),max(MAPS)], /addcmd,/top, title = title,position = [0.1,0.93,0.9,0.94],tcharsize = 1.
	;cgText, 'Ca II K Wing (393.4 - 0.14 nm)', /Normal, Alignment=0.5, charsize = 1.,/addcmd, /place
	;!X.OMargin=0 & !Y.OMargin=0
stop
	cgcontrol, output = figuredir+map+'_partial.eps', /ps_encapsulated, /adjustsize
	cgdelete
;stop
;endfor

stop
end

