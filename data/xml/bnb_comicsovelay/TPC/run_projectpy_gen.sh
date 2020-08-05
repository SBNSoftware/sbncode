#!/bin/bash
for filename in *.xml 
do

echo "running project.py on: $filename"
project.py --xml $filename --stage gen --clean
project.py --xml $filename --stage gen --submit 
project.py --xml $filename --stage gen --check

done