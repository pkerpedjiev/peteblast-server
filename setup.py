import io

from setuptools import find_packages, setup

with io.open("README.rst", "rt", encoding="utf8") as f:
    readme = f.read()

setup(
    name="peteblast",
    version="0.1.0",
    maintainer="Peter Kerpedjiev",
    maintainer_email="pkerpedjiev@gmail.com",
    description="An extra fast BLAST server.",
    long_description=readme,
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=["flask"],
    extras_require={"test": ["pytest", "coverage"]},
)
