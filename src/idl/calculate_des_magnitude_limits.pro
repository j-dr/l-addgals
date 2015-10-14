function calculate_des_magnitude_limits, exposure_time

  es = [9.95, 10.01, 8.02, 6.18, 0.812]  ;flux for a 24th magnitude galaxy
  nk = [4.6, 12.9, 17.7, 45.1, 14.9] ;flux for sky at 24th magnitude


  zps = 24. + 2.5*alog10(es)  ;zerp-point for flux calculation of a galaxy
  apperture = 1.5   ;;angular size of a typical galaxy
  pixel = 0.27      ;;angular size of a detector pixel
  apperture_area = !PI*apperture^2/2
  pixels_area = pixel^2
  npixels = apperture_area/pixels_area
  
  sn = 11.
  flux_gal = 0.5*(sn+sqrt(sn*sn + sn*sn*nk*exposure_time*npixels))
  maglim = -2.5*alog10(flux_gal/exposure_time) + zps 

  return, maglim

end



