#!/bin/bash

# delete every other plotfile lol.
count=0
for p in plotfiles/plt*
do
	if [ $count = 0 ];
	then
		count=1
	elif [ $count = 1 ];
	then
		count=0
		echo "rm -r $p"
	fi
done
