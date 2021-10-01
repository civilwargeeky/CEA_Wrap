rmdir /s /q dist
python setup.py sdist
twine upload dist/*
pause