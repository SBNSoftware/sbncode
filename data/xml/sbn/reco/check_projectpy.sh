#!/bin/bash

for filename in $(find * -name "*.xml" -type f -print )
do
    echo "running project.py on: $filename"
#    project.py --xml $filename --stage sbn --clean
    project.py --xml $filename --stage sbn --submit

    sleep 20

done

