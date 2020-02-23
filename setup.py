import setuptools

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
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
