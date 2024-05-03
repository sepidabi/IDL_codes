;pro obs_time

  ;;Cak
chrom_data = '/scratch/sepid/2016.09.15/CHROMIS/data/08:49:53/Chromis-N/'
crisp_data = '/scratch/sepid/2016.09.15/CRISP/data/08:49:51/Crisp-R/'
scan = [1,2,3,4,5,8,9,11,12,13,14,15,17,18,19,20,21,22,24,25,26,28,29,30,32] ;good-seeing scans
nscan = n_elements(scan)
cn = [0:28]                     ;# of crisp frames
right_h = fltarr(nscan)         ;defined for later
right_ca = fltarr(nscan)         ;defined for later
test_cak = strarr(nscan)
dif_sec_h = fltarr(nscan)
dif_sec_ca = fltarr(nscan)

for n = 0L, nscan-1 do begin
   chrom = file_search(chrom_data+'*000'+strtrim(scan[n]+24,2)+'*3934_+0*') ;looks for the LC scan
   info = red_readhead(chrom)
   chrom_avg = date_conv(strmid(info[22],11,26),'j') ;average time in chromis LC scan
   test_cak[n] = date_conv(chrom_avg,'f')
   chv = date_conv(chrom_avg,'v')
   prob_h = fltarr(n_elements(cn))                     ;defined for later
   prob_ca = fltarr(n_elements(cn))                    ;defined for later
   test_h = strarr(n_elements(cn))
   test_ca = strarr(n_elements(cn))
   for c = 0L, n_elements(cn)-1 do begin
      if c lt 10 then thingy = '0000' else thingy = '000'

      ;H alpha
      crisp_h = file_search(crisp_data+'*'+thingy+strtrim(cn[c],2)+'*6563_+0000*') ;looks for the LC scanS
      lcn_h = n_elements(crisp_h)
      cr_beg_h = strmid((red_readhead(crisp_h[0]))[16],11,26) ;Begining of the LC scan
      cr_end_h = strmid((red_readhead(crisp_h[lcn_h-1]))[17],11,26) ;End of the LC scan
      cr_avg_h = mean([date_conv(cr_beg_h,'j'), date_conv(cr_end_h,'j')]) ;average time in Ha LC scan
      prob_h[c] = abs(cr_avg_h-chrom_avg) ;array of the time difference
      test_h[c] = date_conv(cr_avg_h,'f')
      
      ;Ca 8542
      crisp_ca = file_search(crisp_data+'*'+thingy+strtrim(cn[c],2)+'*8542_+0000*') ;looks for the LC scanS
      lcn_ca = n_elements(crisp_ca)
      cr_beg_ca = strmid((red_readhead(crisp_ca[0]))[16],11,26) ;Begining of the LC scan
      cr_end_ca = strmid((red_readhead(crisp_ca[lcn_ca-1]))[17],11,26) ;End of the LC scan
      cr_avg_ca = mean([date_conv(cr_beg_ca,'j'), date_conv(cr_end_ca,'j')]) ;average time in Ha LC scan
      prob_ca[c] = abs(cr_avg_ca-chrom_avg) ;array of the time difference
      test_ca[c] = date_conv(cr_avg_ca,'f')

;stop
   endfor

   time_cak = date_conv(test_cak[n],'v')
   sec_cak = time_cak[2]*3600.+time_cak[3]*60.+time_cak[4]
   
   dif_h = min(prob_h, can_h)
   right_h[n] = dif_h ;Ha
   time_h = date_conv(test_h[can_h],'v')
   sec_h = time_h[2]*3600.+time_h[3]*60.+time_h[4]
   dif_sec_h[n] = abs(sec_h - sec_cak)
   
   dif_ca = min(prob_ca, can_ca)
   right_ca[n] = dif_ca ;Ca8
   time_ca = date_conv(test_ca[can_ca],'v')
   sec_ca = time_ca[2]*3600.+time_ca[3]*60.+time_ca[4]
   dif_sec_ca[n] = abs(sec_ca - sec_cak)

   print, 'Ca K: i='+strtrim(scan[n],2)+'  '+test_cak[n],'     Ha: i='+strtrim(can_h,2)+'  '+test_h[can_h] ,'     Ca 8542: i='+strtrim(can_ca,2)+'  '+test_ca[can_ca]
   
;   stop
endfor
;ha = min(right_h,fr_h)
;print, 'Best chrom frame for H = ', scan[fr_h]

;ha = min(right_ca,fr_ca)
;print, 'Best chrom frame for Ca 8542 = ', scan[fr_ca]
;stop
nam = strarr(25)
for na = 0L,24 do nam[na] = strtrim(scan[na],2)
cgwindow
cgbarplot,dif_sec_h,barnames = nam,colors = 'red',barspace = 0.5,/add, xtitle = 'Ca II K good-seeing scans (#)', ytitle = 'time difference (s)', title = 'Closest Frame in Time'
cgbarplot,dif_sec_ca,barnames = replicate(' ',25),/overplot,baroffset = 2,barspace = 0.5,color = 'blue',/add
cglegend, colors = ['red','blue'],psym = 15,/addcmd, titles = ['H$\alpha$','Ca 8542'], symsize = 2.,/cent, loc = [0.8, 0.85],length = 0

cgcontrol, output = '/scratch/sepid/DATA/AR/plage/2016.09.15/obs_time.pdf',resize = [900,350],ps_charsize = 1.
stop
end



;  fr = '00025'
;  line = '3934'
;  bw = '-1331'
;  lc = '+0'
;  rw= '+1331'

;  bw_file = file_search(chrom_data+'*'+fr+'*'+line+'*'+bw+'*')
;  rw_file = file_search(chrom_data+'*'+fr+'*'+line+'*'+rw+'*')
;  lc_file = file_search(chrom_data+'*'+fr+'*'+line+'*'+lc+'*')

;  bw_info = red_readhead(bw_file)
;  rw_info = red_readhead(rw_file)
;  lc_info = red_readhead(lc_file)

;  print, '==============='
;  print, '<< Ca II K time >>'
;  print, '--------------------------'
;  print, bw_info[25] ;xposure
;  ;print, '--------------------------'
;  ;print, bw_info[26] ;cadence
;  print, '--------------------------'
;  print, 'SCAN DURATION'
;  print, rw_info[20],bw_info[21]
;  print, '--------------------------'
;  print, 'LINE-CENTER AVG. TIME'
;  print, lc_info[22]

  ;;Crisp
;  crisp_data = '/scratch/sepid/2016.09.15/CRISP/data/08:49:51/Crisp-T/'

;  fr = '00009'

  ;;Ha
;  line = '6563'
;  bw = '-1550'
;  lc = '+0'
;  rw= '+1550'

  ;bw_file = file_search(crisp_data+'*'+fr+'*'+line+'*'+bw+'*lc4*') ;fr 9
  ;rw_file = file_search(crisp_data+'*'+fr+'*'+line+'*'+rw+'*lc4*') ;fr 9
  ;bw_file = file_search(crisp_data+'*'+fr+'*'+line+'*'+bw+'*lc4*') ;fr10
  ;rw_file = file_search(crisp_data+'*'+fr+'*'+line+'*'+rw+'*lc4*') ;fr 10
;  bw_file = file_search(crisp_data+'*'+fr+'*'+line+'*'+bw+'*lc4*') ;fr 11
;  rw_file = file_search(crisp_data+'*'+fr+'*'+line+'*'+rw+'*lc4*') ;fr 11

;  began = bw_file(0)
;  ended = rw_file(n_elements(rw_file)-1)

;  bw_info = red_readhead(began)
;  rw_info = red_readhead(ended)
  ;lc_info = red_readhead(lc_file)

;  print, '==============='
;  print, '<< H alpha time >>'
;  print, '--------------------------'
;  print, bw_info[18] ;xposure
;  print, '--------------------------'
;  ;print, bw_info[26] ;cadence
;  ;print, '--------------------------'
;  print, 'SCAN DURATION'
;  print, bw_info[16],rw_info[17]
;  ;print, '--------------------------'
;  ;print, 'LINE-CENTER AVG. TIME'
;  ;print, lc_info[22]
  
  ;;Ca8
;  line = '8542'
;  bw = '-1700'
;  lc = '+0'
;  rw= '+1700'

  ;bw_file = file_search(crisp_data+'*'+fr+'*'+line+'*'+bw+'*lc0*') ;fr 9
  ;rw_file = file_search(crisp_data+'*'+fr+'*'+line+'*'+rw+'*lc3*') ;fr 9
  ;bw_file = file_search(crisp_data+'*'+fr+'*'+line+'*'+bw+'*lc0*') ;fr 10
  ;rw_file = file_search(crisp_data+'*'+fr+'*'+line+'*'+rw+'*lc3*') ;fr 10
;  bw_file = file_search(crisp_data+'*'+fr+'*'+line+'*'+bw+'*lc0*') ;fr 11
;  rw_file = file_search(crisp_data+'*'+fr+'*'+line+'*'+rw+'*lc3*') ;fr 11

;  began = bw_file(0)
;  ended = rw_file(n_elements(rw_file)-1)

;  bw_info = red_readhead(began)
;  rw_info = red_readhead(ended)
  ;lc_info = red_readhead(lc_file)

;  print, '==============='
;  print, '<< Ca 8542 time >>'
;  print, '--------------------------'
;  print, bw_info[18] ;xposure
;  print, '--------------------------'
;  ;print, bw_info[26] ;cadence
;  ;print, '--------------------------'
;  print, 'SCAN DURATION'
;  print, bw_info[16],rw_info[17]
;  ;print, '--------------------------'
;  ;print, 'LINE-CENTER AVG. TIME'
;  ;print, lc_info[22]
;  print, '============'

  ;cak_time =
;  stop
;end
