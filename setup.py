import os
import setuptools
from distutils.extension import Extension


def scandir(dir, files=[]):
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if os.path.isfile(path) and path.endswith(".pyx"):
            files.append(path.replace(os.path.sep, ".")[:-4])
        elif os.path.isdir(path):
            scandir(path, files)
    return files


def scancppdir(dir, files=[]):
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if os.path.isfile(path) and path.endswith(".cpp"):
            files.append(path)
    return files


# generate an Extension object from its dotted name
def makeExtension(extName, cpp, use_cython):
    extPath = extName.replace(".", os.path.sep)
    if use_cython:
        extPath += ".pyx"
    else:
        extPath += ".c"
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


# get the list of extensions
extNames = scandir("cython_backend")
cpp = scancppdir("cython_backend/c++_src/")
# and build up the set of Extension objects
ext_modules = [makeExtension(name, cpp, use_cython) for name in extNames]


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gaddlemaps", # Replace with your own username
    version="0.0.1",
    author='Hadrian Montes, Jose Manuel Otero Mato',
    author_email='hadrianmontes@gmail.com, josemanuel.otero.mato@gmail.com',
    description="Python package with to apply the GADDLE-MAPS algorithm",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=['gaddlemaps', 'gaddlemaps.parsers', 'gaddlemaps.components'],
    cmdclass=cmdclass,
    ext_modules=ext_modules,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
