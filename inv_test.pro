pro inv_test

cgdelete,/all
fn = 28                       ;desired frame

datadir = '/scratch/sepid/DATA/AR/plage/2016.09.15/'
fibdir = 'fr'+strtrim(long(fn),2)+'/'
savedir = '/nadir-scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/'

pref = ['6302','8542','3950']
cube = ['crispex.*.6302.*_CHROMIS.fcube' , 'crispex.*.8542.*_CHROMIS.fcube' , 'crispex_3950_*_time-corrected.fcube']
sp_cube = ['crispex.*.6302.*_CHROMIS_sp.fcube' , 'crispex.*.8542.*_CHROMIS_sp.fcube' , 'crispex_3950_*_time-corrected_sp.fcube']

for i = 0L, n_elements(pref)-1 do begin

  ;Image Dimensions
  lp_header, datadir+cube[i], nx = nxx, ny = nyy ;image dimensions
  lp_header, datadir+sp_cube[i], nx = nw, ny = ntt ;frame number and wavelength

  ;;fibrils file
  file = file_search(datadir+fibdir+'*'+pref[i]+'*.csav')
  path_n = n_elements(file)  ; number of all the defined paths
  fib_n = path_n/2              ; number of fibrils


  stop
endfor

stop
end
