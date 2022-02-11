import glob

import setuptools
import yaml

# Version
version = ""
fversion = glob.glob("**/_version.py", recursive=True)[0]
with open(fversion) as fid:
    lines = fid.read().splitlines()
    version = lines[1].split("=")[-1].strip().replace('"', "")

# App name - dependencies
env = {}
with open("environment.yml") as fid:
    env = yaml.safe_load(fid)
name = env["name"]
install_requires = env["dependencies"]

setuptools.setup(
    name=name,
    version=version,
    description="Analysis fusion from DNA genomics",
    author=["guillaume-gricourt"],
    author_email=["guipagui@gmail.com"],
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    include_package_data=True,
    install_requires=install_requires,
    entry_points={"console_scripts": ["hmnfusion=hmnfusion.__main__:main"]},
)
