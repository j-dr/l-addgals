pro save_bcc_magnitudes, inmag, outfile, system=system, delete=delete

case (strupcase(system)) of
  'SDSS':begin
    sxaddpar, hdr, 'EXTNAME', 'SDSS'
    sxaddpar, hdr, 'IND_U',0
    sxaddpar, hdr, 'IND_G',1
    sxaddpar, hdr, 'IND_R',2
    sxaddpar, hdr, 'IND_I',3
    sxaddpar, hdr, 'IND_Z',4
  end
  'VISTA':begin
    sxaddpar, hdr, 'EXTNAME', 'Vista'
    sxaddpar, hdr, 'IND_J',0
    sxaddpar, hdr, 'IND_H',1
    sxaddpar, hdr, 'IND_K',2
  end
  'JOHNSON':begin
    sxaddpar, hdr, 'EXTNAME', 'Johnson'
    sxaddpar, hdr, 'IND_U',0
    sxaddpar, hdr, 'IND_B',1
    sxaddpar, hdr, 'IND_V',2
    sxaddpar, hdr, 'IND_R',3
    sxaddpar, hdr, 'IND_I',4
  end
  'DEEP':begin
    sxaddpar, hdr, 'EXTNAME', 'DEEP2'
    sxaddpar, hdr, 'IND_B',0
    sxaddpar, hdr, 'IND_R',1
    sxaddpar, hdr, 'IND_I',2
  end
  'FLAMEX':begin
    sxaddpar, hdr, 'EXTNAME', 'FLAMEX'
    sxaddpar, hdr, 'IND_Bw',0
    sxaddpar, hdr, 'IND_R',1
    sxaddpar, hdr, 'IND_I',2
    sxaddpar, hdr, 'IND_J',3
    sxaddpar, hdr, 'IND_H',4
    sxaddpar, hdr, 'IND_Ks',5
    sxaddpar, hdr, 'IND_3.6', 6
    sxaddpar, hdr, 'IND_4.5', 7
  end
  'CFHTLS':begin
    sxaddpar, hdr, 'EXTNAME', 'CFHTLS'
    sxaddpar, hdr, 'IND_U', 0
    sxaddpar, hdr, 'IND_G', 0
    sxaddpar, hdr, 'IND_R', 0
    sxaddpar, hdr, 'IND_I', 0
    sxaddpar, hdr, 'IND_Z', 0
  end
  'EUCLID':begin
    sxaddpar, hdr, 'EXTNAME', 'EUCLID'
    sxaddpar, hdr, 'IND_VIS', 0
    sxaddpar, hdr, 'IND_Y', 1
    sxaddpar, hdr, 'IND_J', 2
    sxaddpar, hdr, 'IND_H', 3
  end
  'HSC':begin
    sxaddpar, hdr, 'EXTNAME', 'HSC'
    sxaddpar, hdr, 'IND_G', 0
    sxaddpar, hdr, 'IND_R', 1
    sxaddpar, hdr, 'IND_I', 2
    sxaddpar, hdr, 'IND_Z', 3
    sxaddpar, hdr, 'IND_Y', 4
  end
  'LSST':begin
    sxaddpar, hdr, 'EXTNAME', 'LSST'
    sxaddpar, hdr, 'IND_U',0
    sxaddpar, hdr, 'IND_G',1
    sxaddpar, hdr, 'IND_R',2
    sxaddpar, hdr, 'IND_I',3
    sxaddpar, hdr, 'IND_Z',4
    sxaddpar, hdr, 'IND_Y2',5
    sxaddpar, hdr, 'IND_Y3',6
    sxaddpar, hdr, 'IND_Y4',7
  end
  'WFIRST':begin
    sxaddpar, hdr, 'EXTNAME', 'WFIRST'
    sxaddpar, hdr, 'IND_Y',0
    sxaddpar, hdr, 'IND_J',1
    sxaddpar, hdr, 'IND_H',2
    sxaddpar, hdr, 'IND_K',3
  end
  'WFIRST':begin
    sxaddpar, hdr, 'EXTNAME', 'WFIRST'
    sxaddpar, hdr, 'IND_J',0
    sxaddpar, hdr, 'IND_H',1
    sxaddpar, hdr, 'IND_Ks',2
  end
  'GALZ':begin
    sxaddpar, hdr, 'EXTNAME', 'GALZ'
    sxaddpar, hdr, 'IND_U',0
    sxaddpar, hdr, 'IND_G',1
    sxaddpar, hdr, 'IND_R',2
    sxaddpar, hdr, 'IND_I',3
    sxaddpar, hdr, 'IND_Z',4
  end 
  else: begin
    print,'System undefined: ',system
    print, 'Not creating any header information.'
  end
endcase


mag = rename_tags(inmag, 'omag', 'tmag')

mwrfits, mag, outfile, hdr, /create

if KEYWORD_SET(delete) then inmag = 0

end

