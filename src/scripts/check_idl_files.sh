for PIX in 0 1 4 5 6 7 12 13 14 15 16 17 24 25 26 27 28 29 30 31 40 41 42 43 44 45 46 47 48 56 57 58 59 60 61 62 63 72 73 74 75 76 77 78 79 80 88 89 90 91 92 93 94 95
do
for Z in
do
if [ ! -f ${PIX}/${Z}/idl/PO__${PIX}.${Z}_twomass.fit ]
then
echo "Missing IDL files for pixel $PIX and redshift shell $Z"
fi
done
done
