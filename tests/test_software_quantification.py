import filecmp
import json
import os
import re
import sys
import tempfile
from enum import Enum
from typing import Tuple

from main_test import Main_test


class InputType(Enum):
    hmnfusion = 1
    region = 2


class CompareQuantification(Enum):
    files = 1
    vaf = 2


def find_vaf(path: str):
    with open(path) as fid:
        text = fid.read()
    m = re.search(r'(VAF=\d{2},\d{2})', text)
    return m.group(1)


def quantification(sample, bed: str, input_type: InputType, compare: CompareQuantification):
    with tempfile.NamedTemporaryFile(delete=False) as fd:
        args = ["hmnfusion", "quantification"]
        if input_type == InputType.hmnfusion:
            args += ["--hmnfusion-json", sample.extractfusion]
        elif input_type == InputType.region:
            args += ["--region", sample.region]
        args += ["--input-bam", sample.bam]
        args += ["--input-bed", bed]
        args += ["--name", sample.name]
        args += ["--output-vcf", fd.name]

        ret = Main_test.launch(args)
        if ret.returncode > 0:
            print(ret.stderr)
            print(ret.stdout)
            sys.exit(1)

        if compare == CompareQuantification.files:
            sim = filecmp.cmp(fd.name, sample.quantification)
        elif compare == CompareQuantification.vaf:
            res = find_vaf(fd.name)
            theorical = find_vaf(sample.quantification)
            sim = res == theorical

    # Clean up
    os.remove(fd.name)
    return sim


class Test_software(Main_test):
    def test_quantification_hmnfusion(self):
        sim = quantification(
            self.sample_m, self.bed_bcr_path, InputType.hmnfusion, CompareQuantification.files
        )
        self.assertTrue(sim)
        sim = quantification(
            self.sample_p, self.bed_bcr_path, InputType.hmnfusion, CompareQuantification.files
        )
        self.assertTrue(sim)

    def test_quantification_region(self):
        sim = quantification(
            self.sample_p, self.bed_bcr_path, InputType.region, CompareQuantification.vaf
        )
        self.assertTrue(sim)
