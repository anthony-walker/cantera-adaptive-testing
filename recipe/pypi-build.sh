#!/bin/bash
rm -rf ./dist ./cantera_adaptive_testing.egg-info ./build
python3 setup.py sdist bdist_wheel
python3 -m twine upload dist/* #--repository testpypi dist/*
pip uninstall cantera_adaptive_testing
