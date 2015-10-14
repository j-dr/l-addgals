function get_aardvark_redshift_range, fname

boxlen = strmid(fname, 9,4)
num = strmid(fname, 35,3)

if (boxlen eq '1050') then begin
  if (num eq '000') then zrange = [0.0, 0.167160]
  if (num eq '001') then zrange = [0.167160, 0.215040]
  if (num eq '002') then zrange = [0.215040, 0.249760]
  if (num eq '003') then zrange = [0.249760, 0.277760]
  if (num eq '004') then zrange = [0.277760, 0.301840]
  if (num eq '005') then zrange = [0.301840, 0.36]
  if (num eq '006') then zrange = [0.35, 0.36]
  if (num eq '007') then zrange = [0.35, 0.36]
endif
if (boxlen eq '2600') then begin
  if (num eq '000') then zrange = [0.34, 0.374675]
  if (num eq '001') then zrange = [0.374675, 0.418820]
  if (num eq '002') then zrange = [0.418820, 0.456485]
  if (num eq '003') then zrange = [0.456485, 0.489695]
  if (num eq '004') then zrange = [0.489695, 0.520070]
  if (num eq '005') then zrange = [0.520070, 0.547610]
  if (num eq '006') then zrange = [0.547610, 0.573530]
  if (num eq '007') then zrange = [0.573530, 0.597425]
  if (num eq '008') then zrange = [0.597425, 0.620510]
  if (num eq '009') then zrange = [0.620510, 0.641975]
  if (num eq '010') then zrange = [0.641975, 0.662630]
  if (num eq '011') then zrange = [0.662630, 0.682070]
  if (num eq '012') then zrange = [0.682070, 0.701105]
  if (num eq '013') then zrange = [0.701105, 0.719330]
  if (num eq '014') then zrange = [0.719330, 0.737150]
  if (num eq '015') then zrange = [0.737150, 0.754160]
  if (num eq '016') then zrange = [0.754160, 0.770765]
  if (num eq '017') then zrange = [0.770765, 0.786965]
  if (num eq '018') then zrange = [0.786965, 0.802760]
  if (num eq '019') then zrange = [0.802760, 0.818150]
  if (num eq '020') then zrange = [0.818150, 0.833135]
  if (num eq '021') then zrange = [0.833135, 0.847715]
  if (num eq '022') then zrange = [0.847715, 0.861890]
  if (num eq '023') then zrange = [0.861890, 0.876065]
  if (num eq '024') then zrange = [0.876065, 0.889835]
  if (num eq '025') then zrange = [0.889835, 0.93]
  if (num eq '026') then zrange = [0.91, 0.93]
  if (num eq '027') then zrange = [0.91, 0.93]
endif
if (boxlen eq '4000') then begin
  if (num eq '000') then zrange = [0.87, 0.930280]
  if (num eq '001') then zrange = [0.930280, 0.985763]
  if (num eq '002') then zrange = [0.985763, 1.03851]
  if (num eq '003') then zrange = [1.03851, 1.08852]
  if (num eq '004') then zrange = [1.08852, 1.13578]
  if (num eq '005') then zrange = [1.13578, 1.18099]
  if (num eq '006') then zrange = [1.18099, 1.22483]
  if (num eq '007') then zrange = [1.22483, 1.26593]
  if (num eq '008') then zrange = [1.26593, 1.30634]
  if (num eq '009') then zrange = [1.30634, 1.34470]
  if (num eq '010') then zrange = [1.34470, 1.38170]
  if (num eq '011') then zrange = [1.38170, 1.41732]
  if (num eq '012') then zrange = [1.41732, 1.45225]
  if (num eq '013') then zrange = [1.45225, 1.48513]
  if (num eq '014') then zrange = [1.48513, 1.51732]
  if (num eq '015') then zrange = [1.51732, 1.54884]
  if (num eq '016') then zrange = [1.54884, 1.57966]
  if (num eq '017') then zrange = [1.57966, 1.62692]
  if (num eq '018') then zrange = [1.62692, 1.70913]
  if (num eq '019') then zrange = [1.70913, 1.78448]
  if (num eq '020') then zrange = [1.78448, 1.85435]
  if (num eq '021') then zrange = [1.85435, 1.92010]
  if (num eq '022') then zrange = [1.92010, 1.98176]
  if (num eq '023') then zrange = [1.98176, 2.04]
endif

return, zrange

end
