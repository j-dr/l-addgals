from distutils.core import setup

setup(
    name='rdel',
    version='0.1dev',
    packages=['rdel',],
    scripts=['bin/adg-train', 'bin/compress-hlists'],
    package_dir={'rdel': 'src/training/rdel'},
    long_description=open('README.md').read(),
    )
