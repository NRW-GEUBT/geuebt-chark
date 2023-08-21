#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json


def main(summary, sheet_out, merged):
    os.makedirs(sheet_out, exist_ok=True)
    mergedlist = []
    for summ in summary:
        with open(summ, "r") as fi:
            report = json.load(fi)
        isolate_id = report["sample"]
        for key in ["sample", "meta_information", "version_information_list", "mlst"]:
            if key in report:
                del report[key]
        dict_out = {
            "isolate_id": isolate_id,
            "characterization": report
        }
        mergedlist.append(dict_out)
        with open(os.path.join(sheet_out, f"{isolate_id}.json"), "w") as fo:
            json.dump(dict_out, fo, indent=4)
    with open(merged, "w") as fo:
        json.dump(mergedlist, fo, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['summary'],
        snakemake.output['sheet_out'],
        snakemake.output['merged'],
    )
