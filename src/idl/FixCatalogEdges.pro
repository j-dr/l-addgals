;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PROGRAM: FixCatalogEdges
; PURPOSE:
;         Reads in the PO octant data and locates halos with galaxies
;         in multiple files.  Updates the edge tag and ngals values.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

this one doesn't work!!!!!!!

look in addgalsPO

locfile = 'file_locations.dat' ;the file with the information about the 
                               ;regions contained in each section of 
                               ;the PO octant files
Cat_Name = 'hv_POv1.03_oct.'   ;base name of galaxy and halo catalogs

rdfloat, locfile, id, ramin, ramax, decmin, decmax, zmin, zmax
id = long(id)
nids = N_ELEMENTS(id)
UpdateNgals = lonarr(1000000) ;- how much does each ngals need to be updated by?
EdgeHalo = lonarr(1000000) ;- how much does each ngals need to be updated by?

for i = 0, nids - 1 do begin
  ThisGalFile = Cat_Name+num+'_galaxies.fit'
  ThisHaloFile = Cat_Name+num+'_halos.fit'
  g = mrdfits(ThisGalFile, 1)
  ngals = N_ELEMENTS(g)
  h = mrdfits(ThisHaloFile, 1)
  nhalos = N_ELEMENTS(h)
  ;;;---create a halo lookuparray
  HaloArray = lonarr(max(h.haloid)+1)
  for hi = 0L, nhalos - 1 do begin
    HaloArray(h(hi).haloid) = hi
  endfor
  for gi = 0L, ngals - 1 do begin
    ThisHalo = g(gi).haloid
    if ((h(ThisHalo).ra ge ramax(i)) or (h(ThisHalo).ra lt ramin(i)) or $
        (h(ThisHalo).dec ge decmax(i)) or (h(ThisHalo).dec lt decmin(i))) then begin
      ;;;---ok, we've found an edge galaxy.  Find where the halo
      ;;;---really is and update the nonedge array
      for hi = 0, nids - 1 do begin
          if ((h(ThisHalo).ra lt ramax(hi)) and (h(ThisHalo).ra ge ramin(hi)) and $
              (h(ThisHalo).dec lt decmax(hi)) and (h(ThisHalo).dec ge decmin(hi))) then break
      endfor
      if (hi eq nids) then begin
        print, "Error!  I couldn't find a file containing halo ", ThisHalo
        hi = 999
;        print, "stopping."
;        stop
      endif
      g(gi).nonedge = hi
      ;;;---do we need to update the ngals count for this halo?
      if (g(gi).rhalo lt g(gi).r200) then UpdateNgals(ThisHalo)++
      EdgeHalo(ThisHalo) = 1
    endif
  endfor
endfor

;;;---Now we update the halo information
for i = 0, nids - 1 do begin
  ThisHaloFile = Cat_Name+num+'_halos_Edge.fit'
  h = mrdfits(ThisHaloFile, 1)
  nhalos = N_ELEMENTS(h)
  ThisGalFile = Cat_Name+num+'_galaxies_Edge.fit'
  g = mrdfits(ThisGalFile, 1)
  ngals = N_ELEMENTS(g)
  for hi = 0L, nhalos - 1 do begin
    ThisHalo = h(hi).haloid
    h(hi).edge = EdgeHalo
    h(hi).ngals += UpdateNgals(ThisHalo)
  endfor
  for gi = 0L, ngals - 1 do begin
    ThisHalo = g(gi).haloid
    g(gi).ngals += UpdateNgals(ThisHalo)
  endfor
endfor
end

