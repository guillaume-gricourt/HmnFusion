import os
import unittest

from hmnfusion import install_software, utils


class Test_install_software(unittest.TestCase):
    # def required() -> find_executable is tested
    def test_get_path_folder(self):
        """Test get_path_folder()"""
        path = install_software.InstallSoftware.get_path_folder()
        self.assertTrue(type(path) == str)
        self.assertTrue(os.path.isdir(path))

    def test_genefuse(self):
        """Test genefuse install/uninstall"""
        path = install_software.InstallSoftware.get_path_folder()
        self.assertTrue(install_software.InstallSoftware.install_genefuse(path=path))
        ret = utils.cmdline(["genefuse", "--help"])
        self.assertEqual(ret.returncode, 0)
        self.assertTrue(install_software.InstallSoftware.uninstall_genefuse())
        self.assertFalse(install_software.InstallSoftware.uninstall_genefuse())

    def test_lumpy(self):
        """Test lumpy install/uninstall"""
        self.assertTrue(install_software.InstallSoftware.install_lumpy())
        self.assertTrue(install_software.InstallSoftware.check_env_conda(name=install_software.InstallSoftware.LUMPY_ENV_NAME))
        self.assertTrue(install_software.InstallSoftware.uninstall_lumpy())
        self.assertFalse(install_software.InstallSoftware.check_env_conda(name=install_software.InstallSoftware.LUMPY_ENV_NAME))
