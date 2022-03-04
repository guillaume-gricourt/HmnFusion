from typing import Any, Dict, List

import snakemake


def run(
    snakefile: str, config: Dict[Any, Any], *args: List[Any], **kwargs: Dict[Any, Any]
) -> bool:
    """Run a Snakemake workflow

    Parameters
    ----------
    snakefile: str
        Path of a Snakefile
    config: Dict[Any, Any]
        A config dictionary
    *args, **kwargs

    Return
    ------
    bool
        True if the worklow complete without errors
    """
    res = snakemake.snakemake(snakefile, config=config, cores=kwargs.get("cores", 1))
    return res
