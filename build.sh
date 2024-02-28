#activate conda environment
source activate knime-build
#build the plugin
rm -rf build
build_python_extension.py --knime-version 4.7 src build