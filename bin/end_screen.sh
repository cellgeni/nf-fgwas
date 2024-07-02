#!/bin/bash

while IFS="" read -r var || [ -n "$var" ]
do
	screen -S "$var" -X at "#" stuff $'\003'
	screen -S "$var" -X at "#" stuff $'\003'
	screen -S "$var" -X at "#" stuff $'\003'
	screen -S "$var" -X quit
done < <(screen -ls | grep Detached | cut -d. -f1 | awk '{print $1}')

