pro mask_galaxies, g, outfile, mask=mask, sfx=sfx

if not KEYWORD_SET(mask) then mask = '/nfs/slac/g/ki/ki11/des/mswanson/DESdata/BCC/simplebccmask/simplebccmask_round82_t_i.pol'
;polyid = '/afs/slac.stanford.edu/u/ki/mbusha/mbusha.2/projects/mangle/dev10/bin/polyid -w'
polyid = '/afs/slac.stanford.edu/u/ki/mbusha/mbusha.2/projects/mangle/dev10/scripts/my_polyid_gals.sh'

if not KEYWORD_SET(sfx) then sfx = 'tmp'
tempfile = outfile+'.tmp'
openw,1,tempfile
niceprintf,1,g.ra, g.dec
close,1

cmd = polyid+' '+mask+' '+tempfile+' '+outfile+' 0 '+sfx
print, cmd
spawn, cmd
spawn, 'rm '+tempfile

end

