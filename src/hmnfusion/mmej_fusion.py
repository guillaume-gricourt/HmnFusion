from typing import List

import pandas as pd
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.dimensions import ColumnDimension, DimensionHolder


def write(
    filename: str, dfs: List[pd.DataFrame], sheet_name: str = "mmej_fusion"
) -> None:
    """Write dataframe to a file.
    Parameters
    ----------
    filename: str
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
        writer = pd.ExcelWriter(filename)
        df.to_excel(writer, sheet_name=sheet_name)
        writer.save()
        return
    # Write output.
    writer = pd.ExcelWriter(filename)
    row = 0
    for df in dfs:
        df.to_excel(writer, sheet_name=sheet_name, startrow=row)
        row += df.shape[0] + 1
    ws = writer.sheets[sheet_name]  # pull worksheet object

    # Adjust width.
    dim_holder = DimensionHolder(worksheet=ws)
    for ix, col in enumerate(ws.columns):
        col_nb = ix + 1
        lengths = []
        for cell in col:
            if cell.value is not None:
                lengths.append(len(str(cell.value)))
        length = max(lengths) + 6
        dim_holder[get_column_letter(col_nb)] = ColumnDimension(
            ws, min=col_nb, max=col_nb, width=length
        )
    ws.column_dimensions = dim_holder

    writer.save()
