FUNCTION wise_imaging_new_routine, obsid, fits_image_to_read_twelve, fits_image_to_read_twentytwo

READCOL, 'LOWMASS_catalog.dat', id, rac, decl, redshift, logM, SFR, r50, r90, D25, colour, format='L,F,F,F,F,F,F,F,F,F', comment='#'
AD = WHERE(id EQ obsid)
ra = rac[AD]
dec = decl[AD]
redshift = redshift[AD]
optical_uv_sfr = SFR[AD]
beamsize = 1.5*D25[AD]

;##################################The twelve micron image photometry#################################################################

img = mrdfits(fits_image_to_read_twelve, 0, hdr)
d=size(img,/dimension)
nx=d[0]
ny=d[1]

ra_sexagesimal = sixty(ra/15)   ;converting degree decimals into degree longtitude and latitude coordinates
dec_sexagesimal = sixty(dec) ;converting degree decimals into degree longtitude and latitude coordinates
degree_radius = ((beamsize/2.0)/3600.0) 
new_ra = ra + degree_radius 
new_dec = dec + degree_radius
ra_radius_sexagesimal = sixty(new_ra/15) 
dec_radius_sexagesimal = sixty(new_dec) 
coordinates = string([ra_radius_sexagesimal, dec_radius_sexagesimal], /PRINT)   ;co-ordinates of radius point
stringad, coordinates, x, y           ; cordinates of radius point in decimals
ADXY, hdr, ra, dec, x_center, y_center  ; coordinates of center in pixel coordinates
ADXY, hdr, x, y, x_radius, y_radius     ; coordinates of radius point in pixel coordinates
radius_in_pixels = abs(long(x_radius)-long(x_center))
print, 'radius in pixels' , radius_in_pixels
;radius_for_photometry = 3.0*radius_in_pixels		;watch with this, could introduce errors
print, x_center, y_center, 'centers of the image'

radin=1.0*radius_in_pixels
radout=5.0*radius_in_pixels
;areacircle = 3.14159*radin*radin
;areashell = (3.14159*radout*radout) - areacircle

totalflux=0.0
background=0.0
error = 0.0  
counter_shell = 0.0

for i=0, nx-1 do begin
	for j=0, ny-1 do begin
		if (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) LT radout) and (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) GT radin) then begin
			background = background + (img[i, j])
			counter_shell = counter_shell + 1.0
		endif else begin 
			totalflux = totalflux
			background = background
		endelse
	endfor
endfor

skylevel = ((background)/(counter_shell))
print, totalflux, background, skylevel, 'answers'


for i=0, nx-1 do begin
	for j=0, ny-1 do begin
		if (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) LT radout) and (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) GT radin) then begin
			error = error + (img[i, j] - skylevel)^2.0
		endif else begin 
			error = error 
		endelse
	endfor
endfor


counter_circle = 0.0
for i=0, nx-1 do begin
	for j=0, ny-1 do begin
		if (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) LT radin) then begin
			radius = sqrt(((abs(i-x_center))*(abs(i-x_center)))+((abs(j-y_center))*(abs(j-y_center))))
			weighting = 1.0  ;exp(-(radius/radius_in_pixels)*(radius/radius_in_pixels)*alog(2))
			totalflux = totalflux + ((img[i, j]-skylevel)*weighting)
			counter_circle = counter_circle + 1.0
		endif else begin 
			totalflux = totalflux
			background = background	
		endelse
	endfor
endfor

answer_twelve_error_flux = (sqrt(error)/sqrt(counter_shell))*sqrt(counter_circle) 

print, totalflux, background, answer_twelve_error_flux,'answers'


real_flux = totalflux

;To convert from WISE counts to magnitudes in the AB system (Oke 1990):
;http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html

mvega = - 2.5*alog10(real_flux) + 17.80 ;;;after testing with amelie this has changed;;18.0 - 0.665 ;aperture correction needed		;exposure time 8.8s
print, 'mvega equals' , mvega
;;;;;;;;;;;;;;;;;;;;;;;;;;;this is to convert the mvega magnitude into magnitudes in the AB system;;;;;;;;;;;;;
twelve_mag = mvega + 5.174

print, 'twelves mag equals', twelve_mag
twelve_flux = (10^((-48.60 - twelve_mag)/2.5));;;;*((3.0E+08)/(12.0E-06))         ; conversion see http://en.wikipedia.org/wiki/AB_magnitude units of F_UV ergs/s/cm2 * the frequency range taken from the WISE
print, twelve_flux					             ; this is the flux seen at our observational point.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;twelve_flux = (29.045/1.2) *(10^(-mvega/2.5))
;print, 'janskys', twelve_flux


distance_in_cm = lumdist(redshift)*(1E+06)*(3.08567758E+18) 		;distance of the galaxy in cm 
twelve_luminosity = (twelve_flux*4*3.14*distance_in_cm*distance_in_cm*((3.0E+08)/(12.0E-06)))/(3.846E+33)	; this now has units of solar luminosity
print, 'luminosity' , twelve_luminosity
answer_twelve = (4.91E-10)*twelve_luminosity        	;this is now the 12 micron star formation rate


answer_twelve_error_flux_mag = (0.434*answer_twelve_error_flux)/(real_flux)            ; at first done wrongly (-2.5*alog10(answer_twelve_error_flux/8.8)) + 18.0 - 0.665 + 5.174 
answer_twelve_error_flux_flux = twelve_flux*2.303*answer_twelve_error_flux_mag             ; at first done wrongly (10^((-48.60 - answer_twelve_error_flux_mag)/2.5))   ;;;;*((3E+08)/(12E-06))     ; *5.07E+14  or 6.429E+14          ; conversion see http://en.wikipedia.org/wiki/AB_magnitude units of F_UV ergs/s/cm2/hz range taken from the GALEX website


distance_in_cm = lumdist(redshift)*(1E+06)*(3.08567758e+18)
answer_twelve_error_flux_luminosity = (twelve_luminosity*answer_twelve_error_flux_flux)/(twelve_flux)        ; this was originally done wrongly (answer_twelve_error_flux_flux*4*3.1415*distance_in_cm*distance_in_cm*((3.0E+08)/(12.0E-06)))/(3.846E+33)		;this now has units of solar luminosity 
answer_twelve_error = (answer_twelve*answer_twelve_error_flux_luminosity)/(twelve_luminosity) +  (0.39E-10)*twelve_luminosity    ; this was originally done wrongly (4.91E-10)*answer_twelve_error_flux_luminosity

signal_to_noise_twelve = (totalflux/(answer_twelve_error_flux))

;###########################################################################################################################
;##################################The twentytwo  image photometry#################################################################
;###########################################################################################################################

READCOL, 'LOWMASS_catalog.dat', id, rac, decl, redshift, logM, SFR, r50, r90, D25, colour, format='L,F,F,F,F,F,F,F,F,F', comment='#'
AD = WHERE(id EQ obsid)
ra = rac[AD]
dec = decl[AD]
redshift = redshift[AD]
optical_uv_sfr = SFR[AD]
beamsize = 1.5*D25[AD]


img = mrdfits(fits_image_to_read_twentytwo, 0, hdr)
d=size(img,/dimension)
nx=d[0]
ny=d[1]

ra_sexagesimal = sixty(ra/15)   ;converting degree decimals into degree longtitude and latitude coordinates
dec_sexagesimal = sixty(dec) ;converting degree decimals into degree longtitude and latitude coordinates
degree_radius = ((beamsize/2.0)/3600.0) 
new_ra = ra + degree_radius 
new_dec = dec + degree_radius
ra_radius_sexagesimal = sixty(new_ra/15) 
dec_radius_sexagesimal = sixty(new_dec) 
coordinates = string([ra_radius_sexagesimal, dec_radius_sexagesimal], /PRINT)   ;co-ordinates of radius point
stringad, coordinates, x, y           ; cordinates of radius point in decimals
ADXY, hdr, ra, dec, x_center, y_center  ; coordinates of center in pixel coordinates
ADXY, hdr, x, y, x_radius, y_radius     ; coordinates of radius point in pixel coordinates
radius_in_pixels = abs(long(x_radius)-long(x_center))
print, 'radius in pixels' , radius_in_pixels
;radius_for_photometry = 3.0*radius_in_pixels		;watch with this, could introduce errors

radin=1.0*radius_in_pixels
radout=5.0*radius_in_pixels
;areacircle = 3.14159*radin*radin
;areashell = (3.14159*radout*radout) - areacircle

totalflux=0.0
background=0.0
error = 0.0  
counter_shell = 0.0

for i=0, nx-1 do begin
	for j=0, ny-1 do begin
		if (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) LT radout) and (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) GT radin) then begin
			background = background + (img[i, j])
			counter_shell = counter_shell + 1.0
		endif else begin 
			totalflux = totalflux
			background = background
		endelse
	endfor
endfor

skylevel = ((background)/(counter_shell))
print, totalflux, background, skylevel, 'answers'

for i=0, nx-1 do begin
	for j=0, ny-1 do begin
		if (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) LT radout) and (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) GT radin) then begin
			error = error + (img[i, j] - skylevel)^2.0
		endif else begin 
			error = error 
		endelse
	endfor
endfor


counter_circle = 0.0
for i=0, nx-1 do begin
	for j=0, ny-1 do begin
		if (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) LT radin) then begin
			radius = sqrt(((abs(i-x_center))*(abs(i-x_center)))+((abs(j-y_center))*(abs(j-y_center))))
			weighting = 1.0  ;exp(-(radius/radius_in_pixels)*(radius/radius_in_pixels)*alog(2))
			totalflux = totalflux + ((img[i, j]-skylevel)*weighting)
			counter_circle = counter_circle + 1.0
		endif else begin 
			totalflux = totalflux
			background = background	
		endelse
	endfor
endfor


answer_twentytwo_error_flux = (sqrt(error)/sqrt(counter_shell))*sqrt(counter_circle) 

print, totalflux, background, 'answers'

real_flux = totalflux



;To convert from WISE counts to magnitudes in the AB system (Oke 1990):
;http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html

mvega = - 2.5*alog10(real_flux) + 12.945    ;changed after testing with amelie 13.0  - 0.616 ;aperture correction not needed	
print, 'mvega equals' , mvega

twentytwo_mag = mvega + 6.62

print, '22 mag equals' , twentytwo_mag
twentytwo_flux = (10^((-48.60 - twentytwo_mag)/2.5));;;;*((3.0E+08)/(22.0E-06))   ; conversion see http://en.wikipedia.org/wiki/AB_magnitude units of F_UV ergs/s/cm2 * the frequency range taken from the WISE
print, twentytwo_flux					             ; this is the flux seen at our observational point.


distance_in_cm = lumdist(redshift)*(1E+06)*(3.08567758E+18) 		;distance of the galaxy in cm
twentytwo_luminosity = (twentytwo_flux*4*3.14*distance_in_cm*distance_in_cm*((3.0E+08)/(22.0E-06)))/(3.846E+33)	; this now has units of solar luminosity
print, 'luminosity', twentytwo_luminosity
answer_twentytwo = (7.50E-10)*twentytwo_luminosity        	;this is now the 22 micron star formation rate



answer_twentytwo_error_flux_mag = (0.434*answer_twentytwo_error_flux)/(real_flux)            ; at first done wrongly (-2.5*alog10(answer_twelve_error_flux/8.8)) + 18.0 - 0.665 + 5.174 
answer_twentytwo_error_flux_flux = twentytwo_flux*2.303*answer_twentytwo_error_flux_mag             ; at first done wrongly (10^((-48.60 - answer_twelve_error_flux_mag)/2.5))   ;;;;*((3E+08)/(12E-06))     ; *5.07E+14  or 6.429E+14          ; conversion see http://en.wikipedia.org/wiki/AB_magnitude units of F_UV ergs/s/cm2/hz range taken from the GALEX website


distance_in_cm = lumdist(redshift)*(1E+06)*(3.08567758e+18)
answer_twentytwo_error_flux_luminosity = (twentytwo_luminosity*answer_twentytwo_error_flux_flux)/(twentytwo_flux)        ; this was originally done wrongly (answer_twelve_error_flux_flux*4*3.1415*distance_in_cm*distance_in_cm*((3.0E+08)/(12.0E-06)))/(3.846E+33)		;this now has units of solar luminosity 
answer_twentytwo_error = (answer_twentytwo*answer_twentytwo_error_flux_luminosity)/(twentytwo_luminosity)  + (0.07E-10)*twentytwo_luminosity      ; this was originally done wrongly (4.91E-10)*answer_twelve_error_flux_luminosity


signal_to_noise_twentytwo = (totalflux/(answer_twentytwo_error_flux))


A = fltarr(12)
A = [obsid, answer_twelve, answer_twentytwo , optical_uv_SFR, answer_twelve_error, answer_twentytwo_error, twelve_mag, twentytwo_mag, answer_twelve_error_flux_mag, answer_twentytwo_error_flux_mag, signal_to_noise_twelve, signal_to_noise_twentytwo ]
return, A


end




