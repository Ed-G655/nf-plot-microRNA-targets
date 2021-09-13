#!/usr/bin/env bash

#Find all *.mirmap.id files
find -L . \
  -type f \
  -name "*.changes" \
| sed "s#.changes#.changes.tsv#" \
| xargs mk
