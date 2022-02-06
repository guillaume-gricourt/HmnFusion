import filecmp
import json
import os
import sys
import tempfile
from enum import Enum
from typing import Tuple

from main_test import Main_test


class GenefuseFmt(Enum):
    json = 1
    html = 2


def load_extractfusion(path: str) -> dict:
    data = {}
    with open(path) as fid:
        data = json.load(fid)
    data.pop("data", None)
    return data


def extract_fusion(sample, genefuse_fmt: GenefuseFmt) -> Tuple[dict, dict]:
    with tempfile.NamedTemporaryFile(delete=False) as fd:
        args = ["hmnfusion", "extractfusion"]
        if genefuse_fmt == GenefuseFmt.json:
            args += ["--genefuse-json", sample.genefuse_json]
        elif genefuse_fmt == GenefuseFmt.html:
            args += ["--genefuse-html", sample.genefuse_html]
        args += ["--lumpy-vcf", sample.lumpy]
        args += ["--output-json", fd.name]

        ret = Main_test.launch(args)
        if ret.returncode > 0:
            print(ret.stderr)
            print(ret.stdout)
            sys.exit(1)

    res = load_extractfusion(fd.name)
    # Theorical
    theorical = load_extractfusion(sample.extractfusion)

    # Clean up
    os.remove(fd.name)
    return res, theorical


class Test_software(Main_test):
    def test_extract_fusion_html(self):
        res, theorical = extract_fusion(self.sample_a, GenefuseFmt.html)
        self.assertEqual(res, theorical)
        res, theorical = extract_fusion(self.sample_m, GenefuseFmt.html)
        self.assertEqual(res, theorical)
        res, theorical = extract_fusion(self.sample_p, GenefuseFmt.html)
        self.assertEqual(res, theorical)

    def test_extract_fusion_json(self):
        res, theorical = extract_fusion(self.sample_m, GenefuseFmt.json)
        self.assertEqual(res, theorical)
        res, theorical = extract_fusion(self.sample_p, GenefuseFmt.json)
        self.assertEqual(res, theorical)
