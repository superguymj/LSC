:orz

python gen.py > a.in
std < a.in > a.out
force < a.in > as.out

fc a.out as.out
if errorlevel 1 pause

goto orz