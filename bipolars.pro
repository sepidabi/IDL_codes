;bipolars.pro
;	purpose = relates the labels of bipolars detected in 'nearby_me.pro'
;		to their information stored in 'all_/alt_features'

; set some initial parameters
;====================================
read,lcore, prompt = 'Line Core Wavelength = '
read,sm_fact,prompt = 'Smoothing Factor = '
read, nthre, prompt = 'Signal threshold factor = '
read, min_size, prompt = 'Min size allowed (pixels) = '

datadir = '/nadir-scratch/sepid/DATA/'
outputdir = '/nadir-scratch/sepid/YAFTA/'

if sm_fact eq 1 then begin
restore, outputdir+strtrim(long(lcore),2)+'_nearby__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'.sav'
endif else begin
restore, outputdir+strtrim(long(lcore),2)+'_nearby__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'__smoothed.sav'
endelse

n_nearby = n_elements(nearby)/2		;number of bipolar elements
nf = n_elements(features)	;total # of all_/alt_features

bipolar = intarr(2)

for l=0L,n_nearby-1 do begin

	poslabel = nearby(0,l)
	neglabel = nearby(1,l)
	bipos1 = features(where(features.label eq poslabel,npos))	;positive side
	bineg1 = features(where(features.label eq neglabel,nneg))	;negative side

	if npos gt 0 then begin

	    ;if their lifetimes are long enough
            if n_elements(bipos1) gt 4 and n_elements(bineg1) gt 4 then begin

	    ;if bipolars appear at once
	    if abs((bipos1.step)(0)-(bineg1.step)(0)) lt 5 then begin

	    ;distance between the features
	    xnegi = (bineg1.x)(0) & xnegf = (bineg1.x)(nneg-1)
	    ynegi = (bineg1.y)(0) & ynegf = (bineg1.y)(nneg-1)
	    xposi = (bipos1.x)(0) & xposf = (bipos1.x)(npos-1)
	    yposi = (bipos1.y)(0) & yposf = (bipos1.y)(npos-1)

	    Ri = sqrt((xposi-xnegi)^2+(yposi-ynegi)^2)

	    Rf = sqrt((xposf-xnegf)^2+(yposF-ynegf)^2)

	    ;if foopoints are getting apart
	    if Rf gt Ri then begin

		print,''
		print,'STOP!'
		print,''
		print,'nearby index = '+strtrim(long(l),2)
		print,''
		print,'Signal Label (+) = '+strtrim(long(poslabel),2)
		print,'Signal Label (-) = '+strtrim(long(neglabel),2)

;place evolution plot here!!!!!!!!!!!!!!!!!!!!!!!!!

		;dec = ''		
		;read,dec, prompt = 'Accept the bipolar (y/n)? '

		;if dec eq 'y' then begin

			if n_elements(bipolar) eq 0 then bipolar = [(bipos1.label)(0),(bineg1.label)(0)] $
			else bipolar = [[bipolar],[(bipos1.label)(0),(bineg1.label)(0)]]

		;endif
		
		;stop
            endif
	    endif
	    endif
	endif

endfor

print,''
print,'Saving all parameters...'

save,bipolar,addresses,nearby,CP,features,masks,$
	min_size,dx,sigmaCP,nx,ny,nt,nthre,min_size,filename=outputdir+strtrim(long(lcore),2)+'_bipolar__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'__smoothed.sav'

stop
end

