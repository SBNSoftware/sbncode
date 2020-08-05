#!/bin/bash
for filename in *.xml 
do

echo "running project.py on: $filename"
kx509; voms-proxy-init -noregen -rfc -voms "fermilab:/fermilab/sbnd/Role=Analysis"
#project.py --xml $filename --stage gen --check
project.py --xml $filename --stage g4 --checkana
done