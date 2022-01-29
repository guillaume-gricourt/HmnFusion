import json
import subprocess
import sys
import tempfile

from main_test import Main_test


def launch(args):
    if isinstance(args, str):
        args = args.split()
    ret = subprocess.run(args, capture_output=True, encoding="utf8")
    return ret


def load_extractfusion(path):
    data = {}
    with open(path) as fid:
        data = json.load(fid)
    data.pop("data", None)
    return data


class Test_functional(Main_test):

    def test_extract_fusion(self):
        # TEST-A
        with tempfile.NamedTemporaryFile(delete=False) as fd:
            args = ["hmnfusion", "extractfusion"]
            args += ["--genefuse-html", self.test_a_genefuse_html]
            args += ["--lumpy-vcf", self.test_a_lumpy]
            args += ["--output-json", fd.name]

            ret = launch(args)
            if ret.returncode > 0:
                print(ret.stderr)
                print(ret.stdout)
                sys.exit(1)

        res = load_extractfusion(fd.name)
        # Theorical
        theorical = load_extractfusion(fd.name)

        self.assertEqual(res, theorical)
