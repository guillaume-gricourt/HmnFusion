import sys
import tempfile
from typing import List

import pandas as pd
from tests.main_test import Main_test


def mmej_deletion(ref: str, samples: List[str], theorical: str) -> bool:
    res = False
    with tempfile.NamedTemporaryFile(suffix=".xlsx") as fd:
        args = ["hmnfusion", "mmej-deletion"]
        args += ["--input-sample-vcf"] + samples
        args += ["--input-reference-fasta", ref]
        args += ["--output-hmnfusion-xlsx", fd.name]

        ret = Main_test.launch(args)
        if ret.returncode > 0:
            print(ret.stderr)
            print(ret.stdout)
            sys.exit(1)
        # compare data
        df = pd.read_excel(fd.name)
        df_theorical = pd.read_excel(theorical)
        res = df.equals(df_theorical)
    return res


class Test_software(Main_test):
    def test_mmej_deletion_negative(self):
        self.assertTrue(
            mmej_deletion(
                ref=self.ref_mmej,
                samples=[self.n1_vcf],
                theorical=self.mmej_deletion_n1,
            )
        )
        self.assertTrue(
            mmej_deletion(
                ref=self.ref_mmej,
                samples=[self.n1_vcf, self.n2_vcf],
                theorical=self.mmej_deletion_n1n2,
            )
        )

    def test_mmej_positive(self):
        self.assertTrue(
            mmej_deletion(
                ref=self.ref_mmej,
                samples=[self.n1_vcf, self.p1_vcf, self.p2_vcf],
                theorical=self.mmej_deletion_n1p1p2,
            )
        )
