#!/usr/bin/env bash
# set -Eeu

# Repo URL
bakcharak_repo="https://gitlab.com/bfr_bioinformatics/bakcharak.git"

# Commit hash to use
commit="89ab93f6271a20c323e089aed2e79b204ba15f2a" # ver 3.1.2

# Local directory to save the Repo
local_dir="$HOME/.nrw-geuebt/bakCharak"

# if already exists, wipe it clean and redo clone
[ -d "$local_dir" ] && rm -rf "$local_dir"

# clone and checkout repo
echo "Cloning BakCharak and checking out stable commit"
git clone -q "$bakcharak_repo" "$local_dir"
cd "$local_dir"
git checkout "$commit"

# Patch
# echo "Applying patches"
# patch -s --directory="$local_dir" --strip=1 << 'END'

# Update nev with bakcharak yaml
# then fix broken deps and update abricate since it doesn't find perl @INC otherwise
if [[ $(which mamba) == "" ]];then
  conda env update --prune --file "$local_dir/envs/bakcharak.yaml" -p $CONDA_PREFIX &> "$local_dir/deploy_log.txt"
#   conda update -y -p $CONDA_PREFIX openssl=3.1.4 libstdcxx-ng=13.2.0 libiconv=1.17 abricate &>> "$local_dir/deploy_log.txt"
  conda update y -p $CONDA_PREFIX ncbi-amrfinderplus'>=3.12.8' &>> "$local_dir/deploy_log.txt"
  amrfinder -u
else
  mamba env update --file "$local_dir/envs/bakcharak_extra.yaml" -p $CONDA_PREFIX &> "$local_dir/deploy_log.txt"
#   mamba update -y -p $CONDA_PREFIX openssl=3.1.4 libstdcxx-ng=13.2.0 libiconv=1.17 abricate &>> "$local_dir/deploy_log.txt"
  mamba update -y -p $CONDA_PREFIX ncbi-amrfinderplus'>=3.12.8' &>> "$local_dir/deploy_log.txt"
  amrfinder -u
fi

# Install databases
echo "Installing databases"
bash "$local_dir/scripts/bakcharak_setup.sh" --amrfinder_update --databases --status &>> "$local_dir/deploy_log.txt" 
