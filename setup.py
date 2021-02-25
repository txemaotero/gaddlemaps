#    Gaddlemaps python module.
#    Copyright (C) 2019-2021 José Manuel Otero Mato, Hadrián Montes Campos, Luis Miguel Varela Cabo
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
import os
import setuptools
from distutils.extension import Extension


def scandir(dir, files=[]):
    """
    Search for pyx files.
    """
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if os.path.isfile(path) and path.endswith(".pyx"):
            files.append(path.replace(os.path.sep, ".")[:-4])
        elif os.path.isdir(path):
            scandir(path, files)
    return files


def scancppdir(dir, files=[]):
    """
    Search for cpp files.
    """
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if os.path.isfile(path) and path.endswith(".cpp"):
            files.append(path)
    return files


# generate an Extension object from its dotted name
def makeExtension(extName, cpp, use_cython):
    """
    Create the extensions for external modules.
    """
    extPath = extName.replace(".", os.path.sep)
    if use_cython:
        extPath += ".pyx"
    else:
        extPath += ".cpp"
    return Extension(
        extName,
        [extPath]+cpp,
        language="c++",
        libraries = ['armadillo'],
        extra_link_args=["-L/usr/local/lib"],
        include_dirs = ["cython_backend", "cython_backend/c++_src",
                        "/usr/local/include", "."],
        extra_compile_args = ["-O3", "-Wall", "--std=c++11"],
        )


cmdclass = {}

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True
    cmdclass.update({'build_ext': build_ext})


extNames = scandir("cython_backend")
cpp = scancppdir("cython_backend/c++_src/")
ext_modules = [makeExtension(name, cpp, use_cython) for name in extNames]

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gaddlemaps",
    version="0.0.1",
    author='Jose Manuel Otero Mato, Hadrián Montes Campos, Luis Miguel Varela Cabo',
    author_email='hadrianmontes@gmail.com, josemanuel.otero.mato@gmail.com',
    description="Python package to apply the GADDLE-MAPS algorithm",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/txemaotero/gaddlemaps",
    packages=['gaddlemaps', 'gaddlemaps.parsers', 'gaddlemaps.components'],
    cmdclass=cmdclass,
    include_package_data=True,
    ext_modules=ext_modules,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    entry_points={
        "console_scripts": [
            "gaddlemaps=gaddlemaps._cli:main"
        ]
    }
)
