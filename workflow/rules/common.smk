import os
import time
from snakemake.utils import validate


# Pipeline setup --------------------------------------
version = open(os.path.join(workflow.basedir, "..", "VERSION"), "r").read()
pipe_log = os.path.join(os.getcwd(), "PIPELINE_STATUS")


# General puprose functions --------------------------
def get_local_time():
    return time.asctime(time.localtime(time.time()))


def validate_input_param(path, schema):
    try:
        df = pd.read_csv(path, index_col="sample", sep="\t", engine="python")
        validate(df, schema=schema)
    except FileNotFoundError:
        path = os.path.join(workflow.basedir, "..", ".tests", "integration", path)
        df = pd.read_csv(path, index_col="sample", sep="\t", engine="python")
        validate(df, schema=schema)


# Input functions ------------------------------------
def aggregate_summaries(wildcards):
    checkpoint_output = checkpoints.bakcharak.get(**wildcards).output["outdir"]
    return expand(
        "bakcharak/results/{isolate_id}/report/{isolate_id}.bakcharak.json",
        isolate_id=samples.index,
    )

# Validating config ----------------------------------
validate(config, schema="../schema/config.schema.yaml")
validate_input_param(config["sample_sheet"], schema="../schema/samples.schema.yaml")
