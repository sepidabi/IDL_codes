pro photo_data

   close,/all

  datadir = '/home/seki2695/DATA/photosphere/'

  nthre = 4.5
  size_thre = 12.

  nt = 45
  nx = 1615 - 195 + 1
  ny = 1085 - 0 + 1

  CP = fltarr(nx,ny,nt)
  ;mask = fltarr(nx,ny,nt)
  sigmaCP = fltarr(nt)

  for t = 0L, nt-1 do begin

     CP(*,*,t) = readfits(datadir+'CHROMIS_PHOTOSPHERE_CP_'+strtrim(long(t),2)+'.fits', /silent)
     sigmaCP(t) = stddev(abs(CP[*,*,t]))/sqrt(8.)
print, ' '
     print, 'Frame '+strtrim(long(t),2)+ ' finished... '
     print, 'Noise Level ('+strtrim(nthre,2)+'sigma) = '+strtrim(sigmaCP(t))
     print, ' '

    
     ;CgCONTOUR,/noerase,/onimage,CP(*,*,t), color = 'red',levels = [-nthre*sigmaCP(t),nthre*sigmaCP(t)], thick = 1.5
     ;display_labels, CP(*,*,t), mask
     
     
  endfor

  SAVE, /ALL, FILENAME = datadir + 'photo_data.sav'
  
  stop
end
