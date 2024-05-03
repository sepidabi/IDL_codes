;; Cube dimensions, obtained with lp_header
PRO mapdata, file,file1,d,nx,ny,nt,nw, intc=intc

;wavnumber=StrMid(file,15,4)
;wavfile='wav.'+wavnumber+'.f0'

;; UNpolarized data: Ca K
ns = 1

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

