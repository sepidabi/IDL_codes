pro fibril_4stic

  cgdelete,/all
  ;TO BE DEFINED
  ;w_h = 7                       ;h alpha line core position
  ;w_ca = 7
  fn = 28                       ;desired frame
  
  datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
  savedir = '/scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/'
  fibdir = 'fr'+strtrim(long(fn),2)+'/'
  file_ca = file_search(datadir+fibdir+'*3950*.csav')
  file_h = file_search(datadir+fibdir+'*6563*.csav')

  ;Extracting Fibrils and BGs paths
  for n = 0L, n_elements(file_ca)/2-1 do begin

  ;fibril

  ;Ca K  
  restore, file_ca(2*n)
  F_ca = loop_slab
  F_ca_x = x_coords
  F_ca_y = y_coords
  F_ca_l = loop_size

  ;H alpha
  restore, file_h(2*n)
  F_h = loop_slab
  F_h_x = x_coords
  F_h_y = y_coords
  F_h_l = loop_size

  
  ;fibril background

  ;Ca K
  restore, file_ca(2*n+1)
  B_ca = loop_slab
  B_ca_x = x_coords
  B_ca_y = y_coords
  B_ca_l = loop_size

  stop

endfor
  
  stop
end
