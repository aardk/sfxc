# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'progress.ui'
#
# Created: Tue Mar  1 13:19:53 2022
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Dialog1(object):
    def setupUi(self, Dialog1):
        Dialog1.setObjectName(_fromUtf8("Dialog1"))
        Dialog1.resize(996, 363)
        self.gridLayout = QtGui.QGridLayout(Dialog1)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.hboxlayout = QtGui.QHBoxLayout()
        self.hboxlayout.setSpacing(6)
        self.hboxlayout.setMargin(0)
        self.hboxlayout.setObjectName(_fromUtf8("hboxlayout"))
        self.scanLabel = QtGui.QLabel(Dialog1)
        self.scanLabel.setWordWrap(False)
        self.scanLabel.setObjectName(_fromUtf8("scanLabel"))
        self.hboxlayout.addWidget(self.scanLabel)
        self.scanEdit = QtGui.QLineEdit(Dialog1)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.scanEdit.sizePolicy().hasHeightForWidth())
        self.scanEdit.setSizePolicy(sizePolicy)
        self.scanEdit.setReadOnly(True)
        self.scanEdit.setObjectName(_fromUtf8("scanEdit"))
        self.hboxlayout.addWidget(self.scanEdit)
        spacerItem = QtGui.QSpacerItem(93, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.hboxlayout.addItem(spacerItem)
        self.gridLayout.addLayout(self.hboxlayout, 0, 0, 1, 1)
        self._2 = QtGui.QHBoxLayout()
        self._2.setSpacing(6)
        self._2.setMargin(0)
        self._2.setObjectName(_fromUtf8("_2"))
        self.jobLabel = QtGui.QLabel(Dialog1)
        self.jobLabel.setWordWrap(False)
        self.jobLabel.setObjectName(_fromUtf8("jobLabel"))
        self._2.addWidget(self.jobLabel)
        self.jobEdit = QtGui.QLineEdit(Dialog1)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.jobEdit.sizePolicy().hasHeightForWidth())
        self.jobEdit.setSizePolicy(sizePolicy)
        self.jobEdit.setReadOnly(True)
        self.jobEdit.setObjectName(_fromUtf8("jobEdit"))
        self._2.addWidget(self.jobEdit)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self._2.addItem(spacerItem1)
        self.gridLayout.addLayout(self._2, 0, 1, 1, 1)
        self.hboxlayout1 = QtGui.QHBoxLayout()
        self.hboxlayout1.setSpacing(6)
        self.hboxlayout1.setMargin(0)
        self.hboxlayout1.setObjectName(_fromUtf8("hboxlayout1"))
        self.timeLabel = QtGui.QLabel(Dialog1)
        self.timeLabel.setWordWrap(False)
        self.timeLabel.setObjectName(_fromUtf8("timeLabel"))
        self.hboxlayout1.addWidget(self.timeLabel)
        self.timeEdit = QtGui.QLineEdit(Dialog1)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.timeEdit.sizePolicy().hasHeightForWidth())
        self.timeEdit.setSizePolicy(sizePolicy)
        self.timeEdit.setReadOnly(True)
        self.timeEdit.setObjectName(_fromUtf8("timeEdit"))
        self.hboxlayout1.addWidget(self.timeEdit)
        spacerItem2 = QtGui.QSpacerItem(70, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.hboxlayout1.addItem(spacerItem2)
        self.gridLayout.addLayout(self.hboxlayout1, 1, 0, 1, 1)
        self._3 = QtGui.QHBoxLayout()
        self._3.setSpacing(6)
        self._3.setMargin(0)
        self._3.setObjectName(_fromUtf8("_3"))
        self.subjobLabel = QtGui.QLabel(Dialog1)
        self.subjobLabel.setWordWrap(False)
        self.subjobLabel.setObjectName(_fromUtf8("subjobLabel"))
        self._3.addWidget(self.subjobLabel)
        self.subjobEdit = QtGui.QLineEdit(Dialog1)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.subjobEdit.sizePolicy().hasHeightForWidth())
        self.subjobEdit.setSizePolicy(sizePolicy)
        self.subjobEdit.setReadOnly(True)
        self.subjobEdit.setObjectName(_fromUtf8("subjobEdit"))
        self._3.addWidget(self.subjobEdit)
        spacerItem3 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self._3.addItem(spacerItem3)
        self.gridLayout.addLayout(self._3, 1, 1, 1, 1)
        self.progressBar = QtGui.QProgressBar(Dialog1)
        self.progressBar.setProperty("value", 24)
        self.progressBar.setObjectName(_fromUtf8("progressBar"))
        self.gridLayout.addWidget(self.progressBar, 2, 0, 1, 2)
        self.logEdit = QtGui.QPlainTextEdit(Dialog1)
        self.logEdit.setReadOnly(True)
        self.logEdit.setObjectName(_fromUtf8("logEdit"))
        self.gridLayout.addWidget(self.logEdit, 3, 0, 1, 2)
        self.hboxlayout2 = QtGui.QHBoxLayout()
        self.hboxlayout2.setSpacing(6)
        self.hboxlayout2.setMargin(0)
        self.hboxlayout2.setObjectName(_fromUtf8("hboxlayout2"))
        self.buttonDebug = QtGui.QPushButton(Dialog1)
        self.buttonDebug.setEnabled(True)
        self.buttonDebug.setAutoDefault(False)
        self.buttonDebug.setObjectName(_fromUtf8("buttonDebug"))
        self.hboxlayout2.addWidget(self.buttonDebug)
        spacerItem4 = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.hboxlayout2.addItem(spacerItem4)
        self.buttonAbort = QtGui.QPushButton(Dialog1)
        self.buttonAbort.setAutoDefault(False)
        self.buttonAbort.setDefault(False)
        self.buttonAbort.setObjectName(_fromUtf8("buttonAbort"))
        self.hboxlayout2.addWidget(self.buttonAbort)
        self.gridLayout.addLayout(self.hboxlayout2, 4, 0, 1, 2)

        self.retranslateUi(Dialog1)
        QtCore.QObject.connect(self.buttonAbort, QtCore.SIGNAL(_fromUtf8("clicked()")), Dialog1.abort)
        QtCore.QObject.connect(self.buttonDebug, QtCore.SIGNAL(_fromUtf8("clicked()")), Dialog1.create_debug)
        QtCore.QMetaObject.connectSlotsByName(Dialog1)

    def retranslateUi(self, Dialog1):
        Dialog1.setWindowTitle(_translate("Dialog1", "Progress", None))
        self.scanLabel.setText(_translate("Dialog1", "Scan:", None))
        self.jobLabel.setText(_translate("Dialog1", "      Job ID:", None))
        self.timeLabel.setText(_translate("Dialog1", "Time:", None))
        self.subjobLabel.setText(_translate("Dialog1", "Subjob ID:", None))
        self.buttonDebug.setToolTip(_translate("Dialog1", "Let all SFXC processes dump their state", None))
        self.buttonDebug.setText(_translate("Dialog1", "Create &Debug", None))
        self.buttonDebug.setShortcut(_translate("Dialog1", "Alt+D", None))
        self.buttonAbort.setText(_translate("Dialog1", "&Abort", None))
        self.buttonAbort.setShortcut(_translate("Dialog1", "Alt+A", None))

