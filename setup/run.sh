#tests the install for ssrviz

#upload to pip -> twine upload dist/*

python3 setup_pypi.py bdist_wheel
rm -R test_ssrviz
mkdir test_ssrviz
cd test_ssrviz
virtualenv -p python3 virtual_ssrviz
source virtual_ssrviz/bin/activate
pip install ../dist/ssrviz-0.1.2.12-py3-none-any.whl
ssrviz