import pytest
import unittest

from hmnfusion import install_software, utils


@pytest.mark.skipif(not install_software.InstallSoftware.required(), reason="software required are not available")
class Test_InstallSoftware(unittest.TestCase):
    def test_install_software(self):
        # Install
        args = ["hmnfusion", "install-software"]
        ret = utils.cmdline(args)
        self.assertEqual(ret.returncode, 0)
        # Check uninstall
        ret = utils.cmdline(["genefuse", "--help"])
        self.assertEqual(ret.returncode, 0)
        self.assertTrue(
            install_software.InstallSoftware.check_env_conda(
                name=install_software.InstallSoftware.LUMPY_ENV_NAME
            )
        )
        # Uninstall
        args += ["--uninstall"]
        ret = utils.cmdline(args)
        self.assertEqual(ret.returncode, 0)
        # Check uninstall
        with self.assertRaises(utils.ExecutableNotFound):
            utils.find_executable(name="genefuse")
        self.assertFalse(
            install_software.InstallSoftware.check_env_conda(
                name=install_software.InstallSoftware.LUMPY_ENV_NAME
            )
        )
