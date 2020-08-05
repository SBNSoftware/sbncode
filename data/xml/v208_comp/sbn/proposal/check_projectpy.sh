#!/bin/bash
for filename in *icarus.xml 
do

echo "running project.py on: $filename"
project.py --xml $filename --stage sbn --clean
#project.py --xml $filename --stage sbn --checkana
project.py --xml $filename --stage sbn --submit
sleep 600
done
