#!/bin/bash

for file in "$@"
do
 sed -i -- 's/stroke-width="0.001"/stroke-width="0.01"/g' $file
 inkscape $file --export-png=$file".png" -d 120
 # -b ffffffff
 inkscape $file --export-pdf=$file".pdf"
done
