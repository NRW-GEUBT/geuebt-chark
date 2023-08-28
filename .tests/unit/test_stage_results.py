import os
import sys

from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
from json import load


sys.path.insert(0, os.path.dirname(__file__))


def test_stage_results():
    with TemporaryDirectory() as tmpdir:
        # Modify paths to link to your test data and script
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/stage_results/data")
        expected_path = PurePosixPath(".tests/unit/stage_results/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/stage_results.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)
        
        # run function
        sys.path.insert(0, workdir)
        from stage_results import main # import main from your script
        main(
            summary=[os.path.join(workdir, '16-LI00732-0.bakcharak.json')],
            sheet_out=os.path.join(workdir),
            merged=os.path.join(workdir, 'merged.json'),
        )

        # Testing merged files
        with open(
            os.path.join(workdir, 'merged.json'), 'r'
        ) as res, open(
            os.path.join(expected_path, 'merged.json'), 'r'
        ) as expect:
            assert load(res) == load(expect)
        
        # Testing individual file
        with open(
            os.path.join(workdir, '16-LI00732-0.json'), 'r'
        ) as res, open(
            os.path.join(expected_path, '16-LI00732-0.json'), 'r'
        ) as expect:
            assert load(res) == load(expect)
