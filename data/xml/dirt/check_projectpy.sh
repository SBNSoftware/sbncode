#!/bin/bash
for filename in *.xml 
do

echo "running project.py on: $filename"
project.py --xml $filename --stage gen --check
project.py --xml $filename --stage g4 --check
done