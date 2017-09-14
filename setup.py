from setuptools import setup, find_packages

setup(
    name="crosslinked_monolayer",
    version="0.1.0",
    long_description=__doc__,
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'mbuild',
        'foyer',
    ],
    zip_safe=False,
)
