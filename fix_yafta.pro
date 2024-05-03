; fix_yafta.rpo
; pupose = 	to avoid duplicating of YAFTA detections in
; 		disappearing features and merging/splitting case
; input =	detected features from detect_mass.pro
; output =	same as input but with reduced number of features
;		...to be processed in

datadir = '/home/seki2695/DATA/photosphere/'
outputdir = '/home/seki2695/PROJECTS/YAFTA/'

gap = 1  ;definition of nearby distance

print,' '
read,sm_fact,prompt = 'smoothing factor = '
read, nthre, prompt = 'Signal threshold factor = '
read, min_size, prompt = 'Min size allowed (pixels) = '

if sm_fact eq 1 then begin
restore, outputdir + 'YAFTA_output__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'.sav'
endif else begin
restore, outputdir + 'YAFTA_output__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'__smoothed.sav'
endelse

dis = features(where(features.trm lt 0 and features.label ne features.trm))

stop
end
