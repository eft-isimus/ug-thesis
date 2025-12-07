#!/bin/bash

path="$1"
cd "$path" || { echo "Invalid path"; exit 1; }

for f in *_input.sim; do
    [ -f "$f" ] || continue
    idx=$(echo "$f" | sed -E 's/.*_([0-9]+)ri_input\.sim/\1/')
    mkdir -p "$idx"
    mv "$f" "$idx/"
done
