from typing import List

import pandas as pd
from hmnfusion import utils


class MmejFusion(object):
    """MmejFusion store useful methods to deal with specific attributes of fusion"""

    @classmethod
    def to_excel(
        cls, path: str, dfs: List[pd.DataFrame], sheet_name: str = "mmej_fusion"
    ) -> None:
        """Write dataframe to a file.
        Parameters
        ----------
        path: str
            Path of an output file (xlsx, excel format)
        dfs: List[pd.DataFrame]
            List of dataframe to write
        sheet_name: str (default: mmej_fusion)
            Sheet name in the excel file

        Return
        ------
        None
        """
        # Write empty.
        if len(dfs) == 0:
            df = pd.DataFrame()
            with pd.ExcelWriter(path) as writer:
                df.to_excel(writer, sheet_name=sheet_name, index=False)
            return
        # Write output.
        with pd.ExcelWriter(path) as writer:
            row = 0
            for df in dfs:
                df.to_excel(writer, sheet_name=sheet_name, startrow=row)
                row += df.shape[0] + 1

            # Access the workbook and worksheet
            worksheet = writer.sheets[sheet_name]

            # Format
            utils.adjust_dim_worksheet(worksheet)
