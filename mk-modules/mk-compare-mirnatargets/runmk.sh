#!/usr/bin/env bash

#Find all *.mirmap.id files
find -L . \
  -type f \
  -name "*.ref" \
| sed "s#.ref#.changes#" \
| xargs mk
