#!/usr/bin/bash

file=$1

awk '$0 ~ "^>" { match($1, /^>([^:]+)/, id); filename=id[1]} {print >> filename".fa"}' $file
