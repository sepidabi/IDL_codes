;fix_detect.pro
; -> YAFTA_dualmask.pro edited for CP maps
;+++++++++++++++++++++++++++++++++++++++++

; Tracks both convex & contiguous features, and matches labels.
;===============================================================
device,decompose=0  	; gets color table correct
device,retain=2		; refreshes windows when made foreground	

; set some parameters of tracking run
;====================================
read,sm_fact,prompt = 'smoothing factor = '
read, nthre, prompt = 'Signal threshold factor = '
read, min_size, prompt = 'Min size allowed (pixels) = '

cc = 1.                         ;to define window size

amp_sig = 300.

; get data --- restore only if not done already
;==============================================
datadir = '/home/seki2695/DATA/photosphere/'
outputdir = '/home/seki2695/PROJECTS/YAFTA/'

if sm_fact eq 1 then begin
	restore, datadir+'6302_polarization_maps.sav'
endif else begin
	restore, datadir+'6302_polarization_maps__smoothed.sav'
endelse

IN = 1                          ;set to 1 to restrict to IN area otherwise only edges are cropeed

;Dimensions to crop the image
if IN eq 1 then begin ; -> cropping edges + restricting to IN region
   x1 = 195 & y1 = 0
   x2 = 195+1200 & y2 = 0+800
   nx = x2-x1+1 & ny = y2-y1+1
endif else begin ; -> cropping edges
   x1 = 195 & y1 = 0
   x2 = 1615 & y2 = 1085
endelse

;Getting image dimensions
nn = size(CP(x1:x2,y1:y2,*))
nx = nn[1] & ny = nn[2] & nt = nn[3]

sigmaCP = fltarr(nt)

;threshold = 30                 ; ignore pixels below this
; only track features w/min_size pixels
;dx = 1			; pixel size 
jpg = 1                         ; set to 1 to generate .jpg at each step
ps = 0                          ; set to 1 to generat .ps at each step

for i = 0, nt-1 do begin

    ; for generating graphical output
    ;================================
    IF (keyword_set(JPG) or keyword_set(ps)) THEN BEGIN 
        tmpnum = STRING(1000 + fix(i), FORMAT = '(I4)')
        filenum = STRMID(tmpnum , 1)
        fileid = STRCOMPRESS('file-' + filenum, /REMOVE_ALL)

        if keyword_set(ps) then begin
            square_plot,'ps' 
            device,bits=8,/color,filename=STRCOMPRESS(fileid+'.ps',/re)
        endif	    
    endif

    print, "Tracking step:",string(i+1) 

    ; get current data array (suffix "2" is from
    ; current step, "1" is from prev. step)
    ;===========================================

     noise = CPcont(x1:x2,y1:y2,i)
    
     sigmaCP(i) = stddev(abs(noise))/sqrt(8.)
     threshold = nthre*sigmaCP(i)
 
     img2 = CP(x1:x2,y1:y2,i)

    ; do same for contiguous pixel groupings & their structures
    ;===========================================================

    print,threshold
    contiguous_mask2, float(img2), mask2, threshold=threshold, pad=1

    create_features,img2, mask2, features2, min_size=min_size, $
      vx=vx, vy=vy, dx=dx, peakthreshold=peakthreshold
    features2.step=i+1 ; insert current time step










; HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












    ; after at least two steps, begin matching 
    ;==========================================
    if (i gt 0) then begin 
 
        ; now match contiguous features
        ;================================================
        match_features_v01, features1, features2, $
		mask1, mask2, img1, img2, $
                old_max_label=max_label

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

window, 1, xsize = nx, ysize = ny

cgimage,img1, title = 'Massifs at step'+string(i+1), /axes, xtitle = '[pixel]', ytitle = '[pixel]'
plot_edges,mask1                                                            
plot_labels,features1

    ; Output to .jpg file
    ;=====================
    IF keyword_set(JPG) THEN image = $
      tvread(FILENAME = outputdir+'frames/CPF__step_'+strtrim(long(i+1),2)+'__'+strtrim(long(nthre),2)+'sigma.jpeg', /jpeg, QUALITY = 100, /NODIALOG)

    ; Close output to .ps file
    ;=========================
    if keyword_set(ps) then begin
        device,/cl
        set_plot,'x'
    endif
;stop
endfor

; Include most recent step's features.
;=====================================
features = [features, features1]

; Write data to .sav
;=====================================

if sm_fact eq 1. then begin
save,CP,features,masks,$
	min_size,dx,sigmaCP,nx,ny,nt,nthre,min_size,x1,x2,y1,y2,filename=outputdir+'YAFTA_output__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'.sav'
endif else begin
save,CP,features,masks,$
	min_size,dx,sigmaCP,nx,ny,nt,nthre,min_size,x1,x2,y1,y2,filename=outputdir+'YAFTA_output__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'__smoothed.sav'
endelse
stop
end
