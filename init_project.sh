#!/bin/bash

# create directories
for dir in qc qc_final plots scripts initial gwas pca pgs; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        echo "Created directory: $dir"
    else
        echo "Directory $dir already exists, skipping."
    fi
done

# create symlinks
for file in ~/populationgenomics/project_data/GWAS/gwas_data* ~/populationgenomics/project_data/GWAS/*.txt; do
    target="initial/$(basename "$file")"
    if [ ! -e "$target" ]; then
        ln -s "$file" "$target"
        echo "Created symlink: $target"
    else
        echo "Link $target already exists, skipping."
    fi
done

echo "Initialization complete. Directory structure and data links checked/created."