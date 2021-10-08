rmdir /s /q dist
FOR /d /r . %%d IN ("__pycache__") DO @IF EXIST "%%d" rd /s /q "%%d"
python setup.py sdist --formats=zip
rmdir /s /q CEA_Wrap.egg-info
pause