import filecmp
import json
import re
import sys
import tempfile
from enum import Enum
from typing import Any, Dict

from main_test import Main_test


class InputType(Enum):
    hmnfusion = 1
    region = 2


class OutputFmt(Enum):
    vcf = 1
    json = 2


class CompareQuantification(Enum):
    files = 1
    vaf = 2


def find_vaf(path: str):
    with open(path) as fid:
        text = fid.read()
    m = re.search(r"(VAF=\d{2},\d{2})", text)
    if m is None:
        return None
    return m.group(1)


def simplify_json(path: str) -> Dict[Any, Any]:
    data = json.load(open(path))
    ix_software = None
    for ix, metadata in enumerate(data["data"]):
        if metadata[0] == "software":
            ix_software = ix
    if ix_software:
        del data["data"][ix_software]
    return data


def quantification(
    sample,
    bed: str,
    output_fmt: OutputFmt,
    input_type: InputType,
    compare: CompareQuantification,
) -> bool:
    with tempfile.NamedTemporaryFile() as fd:
        args = ["hmnfusion", "quantification"]
        if input_type == InputType.hmnfusion:
            args += ["--input-hmnfusion-json", sample.extractfusion]
        elif input_type == InputType.region:
            args += ["--region", sample.region]
        args += ["--input-sample-bam", sample.bam]
        args += ["--input-hmnfusion-bed", bed]
        args += ["--name", sample.name]
        if output_fmt == OutputFmt.vcf:
            args += ["--output-hmnfusion-vcf", fd.name]
        elif output_fmt == OutputFmt.json:
            args += ["--output-hmnfusion-json", fd.name]
        print(" ".join(args))
        ret = Main_test.launch(args)
        if ret.returncode > 0:
            print(ret.stderr)
            print(ret.stdout)
            sys.exit(1)

        if output_fmt == OutputFmt.vcf:
            if compare == CompareQuantification.files:
                sim = filecmp.cmp(fd.name, sample.quantification)
            elif compare == CompareQuantification.vaf:
                res = find_vaf(fd.name)
                theorical = find_vaf(sample.quantification)
                sim = res == theorical
        elif output_fmt == OutputFmt.json:
            if input_type == InputType.hmnfusion:
                res = simplify_json(fd.name) == simplify_json(
                    sample.quantification_1_json
                )
            elif input_type == InputType.region:
                res = simplify_json(fd.name) == simplify_json(
                    sample.quantification_2_json
                )
            sim = res
    return sim


class Test_software(Main_test):
    def test_quantification_hmnfusion_vcf(self):
        sim = quantification(
            self.sample_m,
            self.bed_bcr_path,
            OutputFmt.vcf,
            InputType.hmnfusion,
            CompareQuantification.files,
        )
        self.assertTrue(sim)
        sim = quantification(
            self.sample_p,
            self.bed_bcr_path,
            OutputFmt.vcf,
            InputType.hmnfusion,
            CompareQuantification.files,
        )
        self.assertTrue(sim)

    def test_quantification_hmnfusion_json(self):
        sim = quantification(
            self.sample_m,
            self.bed_bcr_path,
            OutputFmt.json,
            InputType.hmnfusion,
            CompareQuantification.files,
        )
        self.assertTrue(sim)
        sim = quantification(
            self.sample_p,
            self.bed_bcr_path,
            OutputFmt.json,
            InputType.hmnfusion,
            CompareQuantification.files,
        )
        self.assertTrue(sim)

    def test_quantification_region_vcf(self):
        sim = quantification(
            self.sample_p,
            self.bed_bcr_path,
            OutputFmt.vcf,
            InputType.region,
            CompareQuantification.vaf,
        )
        self.assertTrue(sim)

    def test_quantification_region_json(self):
        sim = quantification(
            self.sample_p,
            self.bed_bcr_path,
            OutputFmt.json,
            InputType.region,
            CompareQuantification.vaf,
        )
        self.assertTrue(sim)
