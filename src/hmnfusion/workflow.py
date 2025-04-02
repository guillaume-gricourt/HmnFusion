from pathlib import Path
from typing import Any, Dict, List

from snakemake.api import (
    ConfigSettings,
    DeploymentSettings,
    OutputSettings,
    ResourceSettings,
    SnakemakeApi,
    StorageSettings,
)
from snakemake.settings.types import DeploymentMethod


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
    with SnakemakeApi(OutputSettings(verbose=True, show_failed_logs=True, printshellcmds=True)) as snakemake_api:
        workflow_api = snakemake_api.workflow(
            snakefile=Path(snakefile),
            storage_settings=StorageSettings(),
            resource_settings=ResourceSettings(cores=1),
            deployment_settings=DeploymentSettings(deployment_method=[DeploymentMethod.CONDA]),
            config_settings=ConfigSettings(config=config),
        )
        dag_api = workflow_api.dag()
        dag_api.execute_workflow()
    return True
