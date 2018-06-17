#!/bin/bash

for file in "$@"; do
	echo "python3 $file"
	python3 "$file"
done