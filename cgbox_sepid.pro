; script to draw box
; suitable for overplotting
; on coyote graphical plots/maps
; ©sepideh kianfar

pro cgbox_sepid, xmin, xmax, ymin, ymax, linewidth = linewidth, linestyle = linestyle, color = color

  ; regtangle dimensions
  dx = xmax-xmin+1
  dy = ymax-ymin+1
  
  ; array definitions
  bottom = intarr(2,2)
  top = intarr(2,2)
  left = intarr(2,2)
  right = intarr(2,2)

  ; bottom line
  bottom[0,0] = xmin
  bottom[0,1] = xmax
  bottom[1,0] = ymin
  bottom [1,1] = ymin
  cgplot, bottom[0,*], bottom[1,*], /addcmd, /noerase, /axes, /over, thick = linewidth, lineStyle = linestyle, /data, color = color

  ; top line
  top[0,0] =xmin
  top[0,1] = xmax
  top[1,0] = ymax
  top [1,1] = ymax
  cgplot, top[0,*], top[1,*], /addcmd, /noerase, /axes, /over, thick = linewidth, lineStyle = linestyle, /data, color = color

  ; left line
  left[1,0] =  ymin
  left[1,1] = ymax
  left[0,0] = xmin
  left[0,1] = xmin
  cgplot, left[0,*], left[1,*], /addcmd, /noerase, /axes, /over, thick = linewidth, lineStyle = linestyle, /data, color = color

  ; right line
  right[1,0] = ymin
  right[1,1] = ymax
  right[0,0] = xmax
  right[0,1] = xmax
  cgplot, right[0,*], right[1,*], /addcmd, /noerase, /axes, /over, thick = linewidth, lineStyle = linestyle, /data, color = color

  print, xmin, xmax
  print, ymin, ymax
  ;stop
end
