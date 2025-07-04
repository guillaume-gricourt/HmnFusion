import logging
import os
import shutil
import stat
import tempfile
from typing import Optional

from hmnfusion import utils


class InstallSoftware(object):
    """Class to install required softwares: GeneFuse and Lumpy

    Class attributes
    ----------------
    GENEFUSE_URL: str
        Github url to download GeneFuse
    LUMPY_VERSION: str
        lumpy version to install
    LUMPY_ENV_NAME: str
        conda lumpy name environment
    SOFTWARE_REQUIRED: List[str]
        all software names required
    """

    GENEFUSE_URL = (
        "https://github.com/OpenGene/GeneFuse/archive/refs/tags/v0.6.1.tar.gz"
    )
    LUMPY_VERSION = "0.3.1"
    LUMPY_ENV_NAME = "lumpy-sv"
    SOFTWARE_REQUIRED = [
        "wget",
        "tar",
        "make",
        "conda",
    ]

    @classmethod
    def get_path_folder(cls) -> Optional[str]:
        """Get the first writable folder contained in the path

        Return
        ------
        str
            the path of the folder, None if not found
        """
        paths = os.environ.get("PATH", None)
        if paths:
            for path in paths.split(":"):
                if os.access(path=path, mode=os.W_OK):
                    return path
        return None

    @classmethod
    def required(cls) -> bool:
        """Check if all software defined as SOFTWARE_REQUIRED are in the PATH

        Raises
        ------
        ExecutableNotFound
            If exec is not found

        Return
        ------
        bool: True
        """
        for software in cls.SOFTWARE_REQUIRED:
            utils.find_executable(name=software)
        return True

    @classmethod
    def install_genefuse(cls, path: str) -> bool:
        """Install GeneFuse in the path provided

        Parameters
        ----------
        path: str
            the path folder

        Return
        ------
        bool
            True if install is done, False otherwise
        """
        if utils.find_executable(name="genefuse", toraise=False):
            logging.warning("GeneFuse is already in the PATH")
            return True
        tmpdir = tempfile.TemporaryDirectory()
        args = [
            "wget",
            cls.GENEFUSE_URL,
            "-O",
            os.path.join(tmpdir.name, "GeneFuse.tar.gz"),
        ]
        ret = utils.cmdline(args)
        if ret.returncode > 0:
            logging.error("Could not download GeneFuse")
            return False
        args = [
            "tar",
            "-xf",
            os.path.join(tmpdir.name, "GeneFuse.tar.gz"),
            "-C",
            tmpdir.name,
        ]
        ret = utils.cmdline(args)
        if ret.returncode > 0:
            logging.error("Could not deflate GeneFuse.tar.gz")
            return False
        args = ["make", "-C", os.path.join(tmpdir.name, "GeneFuse-0.6.1")]
        ret = utils.cmdline(args)
        if ret.returncode > 0:
            logging.error("Could not compile GeneFuse")
            return False
        shutil.copy(os.path.join(tmpdir.name, "GeneFuse-0.6.1", "genefuse"), path)
        os.chmod(os.path.join(tmpdir.name, "GeneFuse-0.6.1", "genefuse"), stat.S_IXUSR)

        # Verify install.
        utils.find_executable(name="genefuse")

        return True

    @classmethod
    def uninstall_genefuse(cls) -> bool:
        """UnInstall GeneFuse

        Return
        ------
        bool
            True if uninstall is done, False otherwise
        """
        path = shutil.which("genefuse")
        if path is not None:
            os.remove(path)
            return True
        return False

    @classmethod
    def check_env_conda(cls, name: str) -> bool:
        """Verify if an environment is installed

        Parameters
        ----------
        name: str
            Name of the conda environment to search

        Return
        ------
        bool
            True if the environment already exists, False otherwise
        """
        args = ["conda", "env", "list"]
        ret = utils.cmdline(args)
        env_found = False
        for line in ret.stdout.splitlines():
            if line.startswith(name):
                env_found = True
                break
        return env_found

    @classmethod
    def install_lumpy(cls) -> bool:
        """Install Lumpy in an environment conda

        Return
        ------
        bool
            True if install is done, False otherwise
        """
        if cls.check_env_conda(name=cls.LUMPY_ENV_NAME):
            logging.warning(
                "Environment 'lumpy-sv' already exists, install could not be done"
            )
            return True
        args = [
            "conda",
            "create",
            "-n",
            "lumpy-sv",
            "-c",
            "bioconda",
            "-c",
            "conda-forge",
            "--override-channels",
            "--yes",
            "lumpy-sv=" + cls.LUMPY_VERSION,
        ]
        ret = utils.cmdline(args)
        if ret.returncode > 0:
            logging.error("Could not create lumpy environment")
            return False
        return True

    @classmethod
    def uninstall_lumpy(cls) -> bool:
        """UnInstall conda environment

        Return
        ------
        bool
            True if uninstall is done, False otherwise
        """
        if not cls.check_env_conda(name=cls.LUMPY_ENV_NAME):
            logging.error("Environment: %s doesn't exists" % (cls.LUMPY_ENV_NAME,))
            return False
        args = [
            "conda",
            "remove",
            "-n",
            cls.LUMPY_ENV_NAME,
            "--all",
            "--yes",
        ]
        ret = utils.cmdline(args)
        if ret.returncode > 0:
            logging.error("Could not delete lumpy environment")
            return False
        return True
