pro addgals_input_halos_to_fits, infile, outfile, offset=offset, omegam=omegam, omegal=omegal

if not KEYWORD_SET(omegam) then omegam = 0.25
if not KEYWORD_SET(omegal) then omegal = 1.0 - omegam
if not KEYWORD_SET(offset) then offset = 0L

table = generate_z_of_r_table(omegam, omegal)

readcol, infile, id, m200, sigma, ip, x, y, z, vx, vy, vz, siglos,r200, rdel, pid, format = 'l,f,f,l,f,f,f,f,f,f,f,f,f,l'

h1 = create_struct('haloID', 0L, 'M200', 0., 'ngals', 0L, 'r200', 0., 'lumtot', 0., 'lum20', 0., 'lbcg', 0., 'siglos', 0., 'halopx', 0., 'halopy', 0., 'halopz', 0., 'halovx', 0., 'halovy', 0., 'halovz', 0., 'ra', 0., 'dec', 0., 'z', 0., 'n18', 0L, 'n19', 0L, 'n20', 0L, 'n21', 0L, 'n22', 0L, 'edge', 0L)

nh = N_ELEMENTS(id)
h = replicate(h1, nh)
h.haloid = lindgen(nh)+offset
h.m200 = m200
h.ngals = 0
h.r200 = r200
h.lumtot = 0
h.lum20 = 0
h.lbcg = 0
h.siglos = siglos
h.halopx = x
h.halopy = y
h.halopz = z
h.halovx = vx
h.halovy = vy
h.halovz = vz
;h.ra = 180./!PI*atan(y/x)
;h.dec = 180./!PI*atan(z/sqrt(x*x + y*y))
h.ra = atan2(y,x, /deg)
h.dec = atan2(z, sqrt(x*x + y*y), /deg)
h.z = z_of_r(sqrt(x*x + y*y + z*z), table)
h.n18 = 0
h.n19 = 0
h.n20 = 0
h.n21 = 0
h.n22 = 0
h.edge = 0

mwrfits, h, outfile, /create

end
