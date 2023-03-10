import setuptools

with open("README.md", encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    name="cantera_adaptive_testing",
    version="0.0.1",
    author="Anthony Walker",
    author_email="walkanth@oregonstate.edu",
    license='MIT License',
    description="This package is used for testing preconditioner additions to Cantera and running studies on the additions.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/anthony-walker/cantera-adaptive-testing",
    entry_points={
        'console_scripts': ['adaptive-testing=cantera_adaptive_testing.commandline:cmd_line_main', 'adaptive-testing.mpi_run_same=cantera_adaptive_testing.commandline:mpi_run_loop', 'adaptive-utilities=cantera_adaptive_testing.commandline:cmd_line_utils', 'adaptive-yp=cantera_adaptive_testing.commandline:cli_yaml_plotter',
            'check_model=cantera_adaptive_testing.commandline:cli_check_model'
            ]
    },
    packages=setuptools.find_packages(),
    include_package_data=True,
    # scripts=["./scripts/study-1/adaptive-precon.sh",],
    package_data={'cantera_adaptive_testing': ['models/*',]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['cantera', 'mpi4py']
)
