; detect_Bhmass.pro
;+++++++++++++++++++++++++++++++++++++++++

; Tracks both convex & contiguous features, and matches labels.
;===============================================================
device,decompose=0  	; gets color table correct
device,retain=2		; refreshes windows when made foreground	

; set some initial parameters
;====================================
read,sm_fact,prompt = 'Smoothing Factor = '
read,nthre,prompt = 'Signal Threshold Factor = '
read, min_size, prompt = 'Min size allowed (pixels) = '

; get data --- restore only if not done already
;==============================================
datadir = '/nadir-scratch/sepid/DATA/'
outputdir = '/nadir-scratch/sepid/YAFTA/'
nadir = '/nadir-scratch/jaime/CHROMIS/'
QSdir = '2016.10.12_quietsun_10_46_37_cakfine_crisp/'


if sm_fact eq 1 then begin
	restore, datadir+strtrim(long(6302),2)+'_polarization_maps.sav'
endif else begin
	restore, datadir+strtrim(long(6302),2)+'_polarization_maps__smoothed.sav'
endelse

;Getting image dimensions
nn = size(LP)
nx = nn[1] & ny = nn[2] & nt = nn[3]

;dx = 1				; pixel size 
jpg = 1                         ; set to 1 to generate .jpg at each step
ps = 0                          ; set to 1 to generat .ps at each step

;Directing to cubes directory
cd, nadir+QSdir
fcube,f
Bh = fltarr(nx,ny,nt)

for i = 0, nt-1 do begin
    
    print, "Tracking step: "+string(i+1) 

    Bh(*,*,i) = lp_get(f[1], i)
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++
;			Bh DETECTION
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ; get current data array (suffix "2" is from
    ; current step, "1" is from prev. step)
    ;===========================================
     
     img2 = Bh(*,*,i)
     threshold = nthre*stddev(img2)
     
    ; do same for contiguous pixel groupings & their structures
    ;===========================================================

    print, threshold
    contiguous_mask2, float(img2), mask2, threshold=threshold, pad=1, /unipolar

    create_features,img2, mask2, features2, min_size=min_size, $
      vx=vx, vy=vy, dx=dx, peakthreshold=peakthreshold
    features2.step=i+1 ; insert current time step
; HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ; after at least two steps, begin matching 
    ;==========================================
    if (i gt 0) then begin 
 
        ; now match contiguous features
        ;================================================
;        match_features_v01, features1, features2, $
;		mask1, mask2, img1, img2, $
;                old_max_label=max_label

        ; concantenate masks from current step with prev. steps' masks
        ;==========================================================
        masks = [[[masks]],[[mask2]]] ; contiguous masks

        ; concantenate features from current step with prev. steps'
        ;==========================================================
        if (n_elements(features) eq 0) then  features = features1 $
        else features = [features, features1] ; & contiguous

        ; find highest existing label, for unique labels on next step
        ;=============================================================
        max_label = max([features.label, features2.label]) 

    endif else begin
	masks = mask2  	; store current contiguous mask
    endelse

    ; raname current data array, mask, and
    ; feature as same from previous step
    ;=========================================
    img1 = temporary(img2)

    mask1 = temporary(mask2)
    features1 = temporary(features2)

    ; graphical output of features, labels
    ;=========================================

;window, 1, xsize = nx, ysize = ny

;cgimage,img1, title = 'Massifs at step'+string(i+1), /axes, xtitle = '[pixel]', ytitle = '[pixel]',/scale
;plot_edges,mask1                                                            
;plot_labels,features1

;stop
;wdelete,1
endfor

; Include most recent step's features.
;=====================================
features = [features, features1]
featuresBh = features
masksBh = masks
nthreBh = nthre

; Write data to .sav
;=====================================

if sm_fact eq 1. then begin
save,CP,LP,featuresBh,masksBh,$
	min_size,dx,sigmaCP,sigmaLP,nx,ny,nt,nthreBh,filename=outputdir+'YAFTA_output_Bh__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'.sav'
endif else begin
save,CP,LP,featuresBh,masksBh,$
	min_size,dx,sigmaCP,sigmaLP,nx,ny,nt,nthreBh,filename=outputdir+'YAFTA_output_Bh__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'__smoothed.sav'
endelse
stop
end
