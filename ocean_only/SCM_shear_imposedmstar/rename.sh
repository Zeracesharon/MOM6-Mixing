#!/bin/bash
shopt -s nullglob

for entry in -*; do
  # only directories
  if [[ -d $entry ]]; then
    new="S${entry#-}"
    echo "Renaming '$entry' ?~F~R '$new'"
    mv -- "$entry" "$new"
  fi
done
