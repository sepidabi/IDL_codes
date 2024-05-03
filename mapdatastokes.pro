;; Cube dimensions, obtained with lp_heade
PRO mapdatastokes, file,file1,d,nx,ny,nt,nw, intc=intc

;wavnumber=StrMid(file,15,4)
;wavfile='wav.'+wavnumber+'.f0'

;; Polarized data: 6302,8542
ns=4

lp_header, file, nx = nx, ny = ny
lp_header, file1,nx = nw, ny = nt

;; Load wavelength
;wav = f0(wavfile)

;; Open file
openr, lun, file, /get_lun
;; Map file into an integer array
if KEYWORD_SET(intc) then begin
   d = assoc(lun, intarr(nx, ny, nw,ns, /nozero), 512)
endif else begin
   d = assoc(lun, fltarr(nx, ny, nw,ns,/nozero), 512) ;; This reads all images within a time step
endelse
return
end


PRO mapdatastokesspectra, file,file1,x,y,cc

;wavnumber=StrMid(file,15,4)
;wavfile='wav.'+wavnumber+'.f0'

ns=1
;; UNpolarized data: Ca K

lp_header, file, nx = nx, ny = ny
lp_header, file1,nx = nw, ny = nt

;; Load wavelength
;wav = f0(wavfile)

;; Open file
openr, lun, file, /get_lun
;; Map file into an integer array

dims=n_elements(x)

cc=fltarr(dims,nw,nt)
for i=0, nt-1 do begin
d = assoc(lun, fltarr(nx, ny, nw,ns,/nozero), 512) ;; This reads all images within a time step
cube=d[i]
for j=0, dims-1 do cc[j,*,i]=cube[x[j],y[j],*,0]
print,i
endfor

;save,cc, filename='slit8542.csav'

return
end
