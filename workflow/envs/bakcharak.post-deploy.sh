#!/usr/bin/env bash
# set -Eeu

# Repo URL
bakcharak_repo="https://gitlab.com/bfr_bioinformatics/bakcharak.git"

# Commit hash to use
commit="58c29e5c22117a2a284a7689c9b7658510af3b1b"

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
echo "Applying patches"
patch -s --directory="$local_dir" --strip=1 << 'END'
diff --unified --recursive --no-dereference bakCharak-orig/bakcharak.py bakCharak/bakcharak.py
--- bakCharak-orig/bakcharak.py	2023-08-16 17:19:02.769515676 +0200
+++ bakCharak/bakcharak.py	2023-08-18 10:52:29.000000000 +0200
@@ -175,7 +175,7 @@
         'species'           : args.species,
 #
         'params'            : {    
-            'threads': args.threads_sample,
+            'threads': int(args.threads_sample),
             'plots': args.viz,
             'do_prokka': args.prokka,
             'do_bakta': args.bakta,
@@ -534,7 +534,7 @@
 
 
     # other   
-    parser.add_argument('-t', '--threads', help='Number of Threads to use. Note that samples can only be processed sequentially due to the required database access. However the allele calling can be executed in parallel for the different loci, default = 10', default=10, required=False)
+    parser.add_argument('-t', '--threads', help='Number of Threads to use. Note that samples can only be processed sequentially due to the required database access. However the allele calling can be executed in parallel for the different loci, default = 10', default=10, required=False, type=int)
     parser.add_argument('-n', '--dryrun', help='Snakemake dryrun. Only calculate graph without executing anything', default=False, action='store_true',  required=False)
     parser.add_argument('--threads_sample', help='Number of Threads to use per sample, default = 1', default=1, required=False)
     parser.add_argument('--forceall', help='Snakemake force. Force recalculation of all steps', default=False, action='store_true',  required=False)
diff --unified --recursive --no-dereference bakCharak-orig/database_setup.sh bakCharak/database_setup.sh
--- bakCharak-orig/database_setup.sh	2023-08-16 17:27:44.000000000 +0200
+++ bakCharak/database_setup.sh	2023-08-18 10:32:25.222250560 +0200
@@ -88,7 +88,7 @@
         echo "wget -O $REPO_PATH/databases/plasmidblast.tar.gz https://gitlab.bfr.berlin/bfr_bioinformatics/bakcharak_resources/-/raw/main/databases/plasmidblast.tar.gz"
         echo "tar -xzvf $REPO_PATH/databases/plasmidblast.tar.gz"
         wget -O $REPO_PATH/databases/plasmidblast.tar.gz https://gitlab.bfr.berlin/bfr_bioinformatics/bakcharak_resources/-/raw/main/databases/plasmidblast.tar.gz
-        tar -xzvf $REPO_PATH/databases/plasmidblast.tar.gz
+        tar -xzvf $REPO_PATH/databases/plasmidblast.tar.gz -C $REPO_PATH/databases/
     fi
 fi
 
@@ -109,7 +109,7 @@
         echo "wget -O $REPO_PATH/databases/mash.tar.gz https://gitlab.bfr.berlin/bfr_bioinformatics/bakcharak_resources/-/raw/main/databases/mash.tar.gz"
         echo "tar -xzvf $REPO_PATH/databases/mash.tar.gz"
         wget -O $REPO_PATH/databases/mash.tar.gz https://gitlab.bfr.berlin/bfr_bioinformatics/bakcharak_resources/-/raw/main/databases/mash.tar.gz
-        tar -xzvf $REPO_PATH/databases/mash.tar.gz
+        tar -xzvf $REPO_PATH/databases/mash.tar.gz -C $REPO_PATH/databases/
     fi
 fi
 
@@ -126,7 +126,7 @@
     echo "tar -xzf $REPO_PATH/databases/db.tar.gz -C $REPO_PATH/databases/ && mv $REPO_PATH/databases/db $REPO_PATH/databases/platon"
 
 else
-    mkdir -p $REPO_PATH/databases/platon
+    # mkdir -p $REPO_PATH/databases/platon
 
     if [[ -f $REPO_PATH/databases/platon/refseq-plasmids.nhr && $force == false ]]; then
         echo "Skip as data exists already. Use --force to force download"

END

# Update nev with bakcharak yaml 
if [[ $(which mamba) == "" ]];then
  conda env update --prune --file "$local_dir/envs/bakcharak.yaml" -p $CONDA_PREFIX &> "$local_dir/conda_log.txt"
  conda update -y -p $CONDA_PREFIX "python>=3.9" &>> "$local_dir/conda_log.txt"
else
  mamba env update --prune --file "$local_dir/envs/bakcharak.yaml" -p $CONDA_PREFIX &> "$local_dir/conda_log.txt"
  mamba update -y -p $CONDA_PREFIX "python>=3.9" &>> "$local_dir/conda_log.txt"
fi

# Install databases
echo "Installing databases"
bash database_setup.sh "$local_dir"
