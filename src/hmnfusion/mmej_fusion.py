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
            writer = pd.ExcelWriter(path)
            df.to_excel(writer, sheet_name=sheet_name)
            writer.save()
            return
        # Write output.
        writer = pd.ExcelWriter(path)
        row = 0
        for df in dfs:
            df.to_excel(writer, sheet_name=sheet_name, startrow=row)
            row += df.shape[0] + 1
        ws = writer.sheets[sheet_name]  # pull worksheet object
        # Format
        utils.adjust_dim_worksheet(ws)
        # Write
        writer.save()
