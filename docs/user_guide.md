# Guide for users

## Installation

### Conda

Install conda from any distribution, i.e. miniconda.
You can follow the setup guide from the [Bioconda team](https://bioconda.github.io/).

We advise installing the mamaba solver in the base environement to speed up
environments creation.

```bash
conda install mamba -n base -c conda-forge
```

### Run environment

The running environement simply required a recent python version (>= 3.9) and snakemake.
If you followed the steps above just run:

```bash
mamba create -n snakemake snakemake
```

### Install module and databases

Download the [latest realease](https://github.com/NRW-GEUBT/geuebt-charak/releases/latest)
and unpack it.

If you're feeling brave, clone the repository form Github:

```bash
git clone https://github.com/NRW-GEUBT/geuebt-charak
```

Most software and databases dependencies will be installed during the first run.

## Configuration

The configuaration can be defined in two ways:

- either edit and locally save the `config/config.yaml` files and provide its path
  to the snakemake command with the `--configfile` argument

- or provide the parameters directly to the snakemake command with
  `--config <ARGNAME>=<VALUE>`

### User defined parameters

Following arguments must be provided for each run:

| Parameter | Type | Description |
| --- | --- | --- |
| `workdir` | path-like string | Path to the ouptut directory |
| `sample_sheet` | path-like string | Path to the sample sheet in TSV format |
| `species` | string | Species to analyze, can trigger additional species-specific results <br>("Ecoli", "Campylobacter", "Salmonella" or "Listeria"). |

### Optional parameters

Following parameters are optional and will revert to defaults if not set:

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `max_threads_per_job` | integer | 1 | Max number of threads assigned to a single job |
| `bakcharak_path` | path-like string | Default installation in `~/.nrw-geuebt/` | Path to the Bakcharak pipeline folder |

## Usage

The workflow can be started with:

```bash
snakemake --use-conda --conda-prefix <PATH TO CONDA ENVS> --configfile <PATH TO CONFIG> --cores <N>
```

See the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for more information.

## Input formats

See the input formatting specifications in the `config` folder.

### Sequence files

The assemblies must be provided as fasta file.
Wrapped and unwrapped fastas as well as multifastas are allowed.
There are no special requirements for sequence headers.

## Results

Results to be used for the next steps are located in the `staging` folder in the workdir.

### Isolate datasheets

Isolate information are summarized in single JSON files
Note that these are generated only for samples satisfying all filters.
They follow the same structure, here for a single entry:

```json
{
    "isolate_id": "16-LI00732-0",
    "characterization": {
        "amr": {
            "count_amr_genes": 2,
            "count_amr_genes_acquired": 2,
            "count_amr_genes_pointmutations": 0,
            "amr_genes": "fosX;vga(G)",
            "count_amr_bymobility": "chromosome (2)",
            "amrclassinfo": {
                "count_amr_classes": 2,
                "amr_classes": "FOSFOMYCIN;LINCOSAMIDE",
                "count_amr_subclasses": 2,
                "amr_subclasses": "FOSFOMYCIN;LINCOSAMIDE"
            },
            "amr_genes_byclass": "FOSFOMYCIN (fosX);LINCOSAMIDE (vga(G))",
            "count_amr_genes_byclass": "FOSFOMYCIN (1);LINCOSAMIDE (1)",
            "count_amr_genes_bysubclass": "FOSFOMYCIN (1);LINCOSAMIDE (1)",
            "amr_genes_bysubclass": "FOSFOMYCIN (fosX);LINCOSAMIDE (vga(G))",
            "count_amr_phenotypes": 0,
            "amr_phenotypes": "",
            "esbl_detected": {
                "esbl_detected": false
            },
            "stress_genes": {
                "METAL": {
                    "genes": {
                        "genes_by_type": "cadC"
                    },
                    "count": {
                        "count": 1
                    }
                }
            },
            "count_stress_genes_bytype": "METAL (1)",
            "stress_genes_bytype": "METAL (cadC)"
        },
        "virulence": {
            "count_virulencefactors": 35,
            "virulencefactors": "clpC;vip;inlF;pdgA;inlA;inlB;lntA;iap;inlC;fbpA;lspA;lpeA;bsh;hbp2;hbp1;prsA2;lapB;lap;oatA;inlK;inlJ;aut;clpE;lplA1;hpt;clpP;inlP;gtcA;ami;plcB;actA;mpl;hly;plcA;prfA",
            "count_virulencefactors_byclass": "Adherence (7);Exoenzyme (1);Exotoxin (3);Immune_modulation (5);Invasion (6);Motility (1);Nutritional/Metabolic_factor (4);Post-translational_modification (3);Regulation (1);Stress_survival (4)"
        },
        "plasmids": {
            "count_INC": 1,
            "INC_types": "rep25_2_M640p00130(J1776plasmid)",
            "cumulative_plasmid_length": 25677,
            "plasmid_contigs": 1,
            "circular_contigs": 1,
            "replication_elements": 1,
            "mobilization_elements": 0,
            "conjugation_elements": 0,
            "count_OriT": 0,
            "plasmid_fasta": "NA"
        },
        "reference_information": {
            "reference_species": "Listeria monocytogenes",
            "reference_genus": "Listeria",
            "reference_identity": 0.985,
            "reference_accession": "GCF_001565515.1",
            "reference_genome": " Listeria monocytogenes strain LM09-00521 genome assembly"
        },
        "Listeria": {
            "serotype": "1/2c, 3c",
            "PRS": "FULL",
            "LMO0737": "FULL",
            "LMO1118": "FULL",
            "ORF2110": "NONE",
            "ORF2819": "NONE",
            "comment": "NA"
        }
    }
}
```

Note that some fields are species-dependent
