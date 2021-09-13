#!/usr/bin/env bash

#Find all *.mirmap.id files
find -L . \
  -type f \
  -name "*.changes.tsv" \
| sed "s#.changes.tsv#.png#" \
| xargs mk
