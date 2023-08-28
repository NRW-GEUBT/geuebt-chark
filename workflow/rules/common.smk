import os
import time
from snakemake.utils import validate


# Pipeline setup --------------------------------------
version = open(os.path.join(workflow.basedir, "..", "VERSION"), "r").read()
pipe_log = os.path.join(os.getcwd(), "PIPELINE_STATUS")


# Validating config ----------------------------------
validate(config, schema="../schema/config.schema.yaml")


# Loading samples ---------------------
samples = pd.read_csv(
    config["sample_sheet"], index_col="sample", sep="\t", engine="python"
)
samples.index = samples.index.astype("str", copy=False)


# General puprose functions --------------------------
def get_local_time():
    return time.asctime(time.localtime(time.time()))


# Input functions ------------------------------------
def aggregate_summaries(wildcards):
    checkpoint_output = checkpoints.bakcharak.get(**wildcards).output["outdir"]
    return expand(
        "bakcharak/results/{isolate_id}/report/{isolate_id}.bakcharak.json",
        isolate_id=samples.index,
    )
