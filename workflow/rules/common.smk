import os
import time
from snakemake.utils import validate


# Pipeline setup --------------------------------------
version = open(os.path.join(workflow.basedir, "..", "VERSION"), "r").read()
pipe_log = os.path.join(os.getcwd(), "PIPELINE_STATUS")


# Validating config ----------------------------------
validate(config, schema="../schema/config.schema.yaml")


# Loading samples ---------------------
sample_path = config["sample_sheet"]
samples = pd.read_csv(sample_path, index_col="sample", sep="\t", engine="python")
samples.index = samples.index.astype("str", copy=False)


# General puprose functions --------------------------
def get_local_time():
    return time.asctime(time.localtime(time.time()))


# Input functions ------------------------------------
def aggregate_summaries(wildcards):
    checkpoint_output = checkpoints.bakcharak.get(**wildcards).output["outdir"]
    # ids_map = glob_wildcards(
        # os.path.join(checkpoint_output, "results/{isolate_id}")
    # ).isolate_id
    # return expand("staging/charak_sheets/{isolate_id}.json", isolate_id=ids_map)
    return expand("bakcharak/results/{isolate_id}/report/{isolate_id}.bakcharak.json", isolate_id=samples.index)