import sys, json
try:
  from PyQt5.QtWidgets import QApplication, QDialog, QMainWindow, QFileDialog, QMessageBox, QMenu, QAction
  from PyQt5 import QtCore
  from PyQt5.QtCore import Qt
  from PyQt5.QtGui import QPalette, QColor, QCursor
except ImportError:
  print("PyQt5 is not installed! Cannot run gui")
  print('Use "pip install CEA_Wrap[gui]" to install dependencies')
  print('Or just do "pip install PyQt5" yourself')
  sys.exit(1)
from .gui_base import Ui_MainWindow
from .thermo_lib import ThermoInterface
from .utils import _get_data_file


with open(_get_data_file("thermo_elements.json")) as file:
  thermo_elements = json.load(file)

def postLoad():
  elementNameWidgets = [ui.elem1name, ui.elem2name, ui.elem3name, ui.elem4name, ui.elem5name, ]
  elementAmntWidgets = [ui.elem1amnt, ui.elem2amnt, ui.elem3amnt, ui.elem4amnt, ui.elem5amnt, ]
  elementSymbols = [" "] + ["{: 2} {}".format(val["number"], val["symbol"]) for val in thermo_elements]
  for elem in elementNameWidgets:
    elem.addItems(elementSymbols)

if __name__ == "__main__":
  app = QApplication(sys.argv)
  window = QMainWindow()
  ui = Ui_MainWindow()
  ui.setupUi(window)
  window.ui = ui
  postLoad()
  window.show()
  sys.exit(app.exec_())