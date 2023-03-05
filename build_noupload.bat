rmdir /s /q dist
FOR /d /r . %%d IN ("__pycache__") DO @IF EXIST "%%d" rd /s /q "%%d"
python -m build . --wheel
rmdir /s /q CEA_Wrap.egg-info
rmdir /s /q build
pause