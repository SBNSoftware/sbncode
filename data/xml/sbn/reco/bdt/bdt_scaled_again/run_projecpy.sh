#!/bin/bash
for d in */ ; do
   cd $d
   echo "on directory $d" 
   bash check_projectpy.sh
   cd ../
   sleep 1000;

done
