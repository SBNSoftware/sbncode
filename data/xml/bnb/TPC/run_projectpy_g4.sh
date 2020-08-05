#!/bin/bash
for filename in *.xml 
do
#kx509; voms-proxy-init -noregen -rfc -voms 'fermilab:/fermilab/sbnd/Role=Analysis'
echo "running project.py on: $filename"
project.py --xml $filename --stage g4 --clean
#project.py --xml $filename --stage gen --check
#project.py --xml $filename --stage g4 --submit 
#project.py --xml $filename --stage g4 --check

done