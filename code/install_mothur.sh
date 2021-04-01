#!/usr/bin/env bash

# author: Nick Lesniak
# input: none
# outputs: mothur installed in code/mothur
#
# The zip archive contains a director called "mothur" so we can extract it directly
# to code/

wget -P code/mothur/ -nc https://github.com/mothur/mothur/releases/download/v1.44.3/Mothur.linux.zip
unzip -n -d code/ code/mothur/Mothur.linux.zip

if [[ $? -eq 0 ]]
then
	touch code/mothur/mothur
else
	echo "FAIL: were not able to successfully install mothur"
fi