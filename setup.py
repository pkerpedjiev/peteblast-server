import io

from setuptools import find_packages, setup

setup(
    name="peteblast",
    version="0.2.0",
    maintainer="Peter Kerpedjiev",
    maintainer_email="pkerpedjiev@gmail.com",
    description="An extra fast BLAST server.",
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=["flask"],
    extras_require={"test": ["pytest", "coverage"]},
)
