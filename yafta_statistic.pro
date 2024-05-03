  close, /all

  datadir = '/home/seki2695/DATA/photosphere/'
  outputdir = '/home/seki2695/PROJECTS/YAFTA/'

  print,' '
  read, nthre, prompt = 'Signal threshold factor = '
  read, min_size, prompt = 'Min size allowed (pixels) = '

  restore, outputdir+'YAFTA_output__'+strtrim(long(nthre),2)+'sigma__size'+strtrim(long(min_size),2)+'.sav'

  ;set below tatements to either of all_ or alt_ features

  features = alt_features
  mask = alt_masks

;+++++++++++++++++++++++++++++++++++
;            STATISTICS
;+++++++++++++++++++++++++++++++++++

;total flux
notfirst = features(where(features.step gt 0))
flux = total(notfirst.phi)
neg = notfirst(where(notfirst.sign eq -1))
negflux = -1*total(neg.phi)
pos = notfirst(where(notfirst.sign eq 1))
posflux = total(pos.phi)

;A feature that first appears at step t = k where no features previously existed (and
;is therefore not a result of fragmentation [see below]) has −k in its source field. (All
;features from the initial tracking step have .src set to 0.) Analogously, features that
;disappear at a step t = k have .trm set to −k. (All features from the final tracking
;step have .trm set to 0.) To find appearances, and calculate their contribution to the
;average flux (excluding all features from the first step), one could use

apps = where(features.src lt 0, n_apps)
if n_apps ne 0 then appflux = total(features(apps).phi)


;disapperance

disap = where(features.trm lt 0, n_disap)
if n_disap ne 0 then disapflux = total(features(disap).phi)

;A feature at step t = k that derived from a feature from step t = k − 1 (either as a
;one-to-one match, or as a fragment) has the label of the feature from t = k − 1 in its
;.src field. In a fragmentation, though, only the largest fragment (by flux) keeps its .src
;label, so fragments can be found by searching for .src fields that are positive and differ
;from .label fields. (Recall that all features from the initial step have a zero in their .src
;fields!) So, to find the average unsigned flux contained in a fragment, use:

where_frag = where(features.label ne features.src and features.src gt 0, n_frags)
if n_frags ne 0 then frags = features(where_frag)
if n_frags ne 0 then fragflux = total(frags.phi)/n_frags

;merging flux

where_merg = where(features.label ne features.trm and features.trm gt 0, n_mergs)
if n_mergs ne 0 then mergs = features(where_merg)
if n_mergs ne 0 then mergflux = total(mergs.phi)/n_mergs

;average unsigned flux of all elements

avgflux = total(features.phi)/n_elements(features.phi)

;print results
print,' '
if n_apps ne 0 then print, 'flux appeared = '+ strtrim(appflux,2)
if n_disap ne 0 then print, 'flux disappeared = '+strtrim(disapflux,2)
if n_apps ne 0 then print, 'flux fragmented = '+ strtrim(fragflux,2)
if n_disap ne 0 then print, 'flux merged = '+strtrim(mergflux,2)

print,' '
print,'Total Flux = ' + strtrim(flux,2)
print,'Avg. Unsigned Flux in All Features = '+ strtrim(avgflux,2)
print, 'Total negative flux = '+ strtrim(negflux,2)
print, 'Total positive flux = '+ strtrim(posflux,2)
print,' '

;The 1-D pixel addresses of each feature are also stored in that feature’s structure (as
;a string variable), and can be used to visualize a feature’s shape at
;a given step.

;number of all the features
nfeat = n_elements(features)

;number of individual features TRACKED
ntr_feat = max(features.label)

end
