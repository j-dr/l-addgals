dir = '~/ki01/projects/addgals/Carmen/LC/v1.08a/v22/DES_catalog_2.12/'

file1 = dir+'DES_Mock_halos_lambda_mockgmr_07132010.fit'
file2 = dir+'DES_Mock_halos_lambda_mockrmi_07132010.fit'
file3 = dir+'DES_Mock_halos_lambda_mockimz_07132010.fit'
file4 = dir+'DES_Mock_halos_lambda_mockzmy_07142010.fit'
hfile = dir+'DES_Mock_halos.fit'
outfile = dir+'DES_Mock_halos_lambda.fit'

lambda1 = mrdfits(file1, 1)
lambda2 = mrdfits(file2, 1)
lambda3= mrdfits(file3, 1)
lambda4 = mrdfits(file4,1)
h = mrdfits(hfile,1)

add_tag, h, 'lambda', 0., h2
match, h2.haloid, lambda1.haloid, a, b
help, lambda1, b
h2[a].lambda = lambda1[b].lambda_mockgmr_00_090

match, h2.haloid, lambda2.haloid, a, b
help, lambda2, b
h2[a].lambda = lambda2[b].lambda_mockrmi_00_090

match, h2.haloid, lambda3.haloid, a, b
help, lambda3, b
h2[a].lambda = lambda3[b].lambda_mockimz_00_090

match, h2.haloid, lambda4.haloid, a, b
help, lambda4, b
h2[a].lambda = lambda4[b].lambda_mockzmy_00_090

mwrfits, h2, outfile, /create

plotsym, 0, /fill
plot, lambda1.m200, lambda1.lambda_mockgmr_00_090, psym = 1, /xlog,/ylog,$
      yrange = [0.1,1e4], xtitle='M!D200', ytitle = 'lambda'
oplot, lambda2.m200, lambda2.lambda_mockrmi_00_090, psym = 1, color = !blue
oplot, lambda3.m200, lambda3.lambda_mockimz_00_090, psym = 1, color = !green
oplot, lambda4.m200, lambda4.lambda_mockzmy_00_090, psym = 1, color = !red
legend, ['z < 0.35', '0.35 - 0.7', '0.7 - 0.9', '> 0.9'], psym = 8, $
        color=[!P.COLOR, !BLUE, !GREEN, !RED]

!P.MULTI = [0,2,2]
plot, lambda1.m200, lambda1.lambda_mockgmr_00_090, psym = 1, /xlog,/ylog,$
      yrange = [0.1,1e4], xtitle='M!D200', ytitle = 'lambda', title = 'z < 0.35'
oplot, lambda1.m200, lambda1.ngals, psym = 4, color = !red
legend, ['lambda', 'ngals'], psym = [1,4], color = [!P.COLOR, !RED]

plot, lambda2.m200, lambda2.lambda_mockrmi_00_090, psym = 1, /xlog,/ylog,$
      yrange = [0.1,1e4], xtitle='M!D200', ytitle = 'lambda', title = 'z = 0.35 - 0.7'
oplot, lambda2.m200, lambda2.ngals, psym = 4, color = !red
legend, ['lambda', 'ngals'], psym = [1,4], color = [!P.COLOR, !RED]

plot, lambda3.m200, lambda3.lambda_mockimz_00_090, psym = 1, /xlog,/ylog,$
      yrange = [0.1,1e4], xtitle='M!D200', ytitle = 'lambda', title = 'z = 0.7 - 0.9'
oplot, lambda3.m200, lambda3.ngals, psym = 4, color = !red
legend, ['lambda', 'ngals'], psym = [1,4], color = [!P.COLOR, !RED]

plot, lambda4.m200, lambda4.lambda_mockzmy_00_090, psym = 1, /xlog,/ylog,$
      yrange = [0.1,1e4], xtitle='M!D200', ytitle = 'lambda', title = 'z > 0.9'
oplot, lambda4.m200, lambda4.ngals, psym = 4, color = !red
legend, ['lambda', 'ngals'], psym = [1,4], color = [!P.COLOR, !RED]
!P.MULTI = 0



end
