# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'S:\_Lib\python\slab\widgets\SweepDialog.ui'
#
# Created: Mon Jun 18 20:55:33 2012
#      by: PyQt4 UI code generator 4.8.5
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_SweepDialog(object):
    def setupUi(self, SweepDialog):
        SweepDialog.setObjectName(_fromUtf8("SweepDialog"))
        SweepDialog.resize(421, 609)
        SweepDialog.setWindowTitle(QtGui.QApplication.translate("SweepDialog", "Dialog", None, QtGui.QApplication.UnicodeUTF8))
        self.verticalLayoutWidget = QtGui.QWidget(SweepDialog)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(10, 0, 401, 571))
        self.verticalLayoutWidget.setObjectName(_fromUtf8("verticalLayoutWidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.sweep1_groupBox = QtGui.QGroupBox(self.verticalLayoutWidget)
        self.sweep1_groupBox.setTitle(QtGui.QApplication.translate("SweepDialog", "Sweep 1", None, QtGui.QApplication.UnicodeUTF8))
        self.sweep1_groupBox.setObjectName(_fromUtf8("sweep1_groupBox"))
        self.gridLayoutWidget = QtGui.QWidget(self.sweep1_groupBox)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(9, 19, 361, 151))
        self.gridLayoutWidget.setObjectName(_fromUtf8("gridLayoutWidget"))
        self.gridLayout = QtGui.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.param_sweep1 = QtGui.QComboBox(self.gridLayoutWidget)
        self.param_sweep1.setObjectName(_fromUtf8("param_sweep1"))
        self.gridLayout.addWidget(self.param_sweep1, 0, 1, 1, 1)
        self.label_4 = QtGui.QLabel(self.gridLayoutWidget)
        self.label_4.setText(QtGui.QApplication.translate("SweepDialog", "Parameter", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout.addWidget(self.label_4, 0, 0, 1, 1)
        self.sweep1_stackedWidget = QtGui.QStackedWidget(self.gridLayoutWidget)
        self.sweep1_stackedWidget.setObjectName(_fromUtf8("sweep1_stackedWidget"))
        self.page_3 = QtGui.QWidget()
        self.page_3.setObjectName(_fromUtf8("page_3"))
        self.formLayoutWidget = QtGui.QWidget(self.page_3)
        self.formLayoutWidget.setGeometry(QtCore.QRect(10, 20, 341, 91))
        self.formLayoutWidget.setObjectName(_fromUtf8("formLayoutWidget"))
        self.formLayout = QtGui.QFormLayout(self.formLayoutWidget)
        self.formLayout.setMargin(0)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.label_9 = QtGui.QLabel(self.formLayoutWidget)
        self.label_9.setText(QtGui.QApplication.translate("SweepDialog", "Start", None, QtGui.QApplication.UnicodeUTF8))
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.label_9)
        self.label_10 = QtGui.QLabel(self.formLayoutWidget)
        self.label_10.setText(QtGui.QApplication.translate("SweepDialog", "Stop", None, QtGui.QApplication.UnicodeUTF8))
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_10)
        self.param_sweep1_start_int = QtGui.QSpinBox(self.formLayoutWidget)
        self.param_sweep1_start_int.setObjectName(_fromUtf8("param_sweep1_start_int"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.param_sweep1_start_int)
        self.param_sweep1_stop_int = QtGui.QSpinBox(self.formLayoutWidget)
        self.param_sweep1_stop_int.setObjectName(_fromUtf8("param_sweep1_stop_int"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.param_sweep1_stop_int)
        self.label_3 = QtGui.QLabel(self.formLayoutWidget)
        self.label_3.setText(QtGui.QApplication.translate("SweepDialog", "Step", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.label_3)
        self.param_sweep1_step_int = QtGui.QSpinBox(self.formLayoutWidget)
        self.param_sweep1_step_int.setObjectName(_fromUtf8("param_sweep1_step_int"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.param_sweep1_step_int)
        self.sweep1_stackedWidget.addWidget(self.page_3)
        self.page_4 = QtGui.QWidget()
        self.page_4.setObjectName(_fromUtf8("page_4"))
        self.formLayoutWidget_2 = QtGui.QWidget(self.page_4)
        self.formLayoutWidget_2.setGeometry(QtCore.QRect(10, 20, 341, 91))
        self.formLayoutWidget_2.setObjectName(_fromUtf8("formLayoutWidget_2"))
        self.formLayout_2 = QtGui.QFormLayout(self.formLayoutWidget_2)
        self.formLayout_2.setMargin(0)
        self.formLayout_2.setObjectName(_fromUtf8("formLayout_2"))
        self.label_11 = QtGui.QLabel(self.formLayoutWidget_2)
        self.label_11.setText(QtGui.QApplication.translate("SweepDialog", "Start", None, QtGui.QApplication.UnicodeUTF8))
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.formLayout_2.setWidget(0, QtGui.QFormLayout.LabelRole, self.label_11)
        self.label_12 = QtGui.QLabel(self.formLayoutWidget_2)
        self.label_12.setText(QtGui.QApplication.translate("SweepDialog", "Stop", None, QtGui.QApplication.UnicodeUTF8))
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.formLayout_2.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_12)
        self.label_13 = QtGui.QLabel(self.formLayoutWidget_2)
        self.label_13.setText(QtGui.QApplication.translate("SweepDialog", "Step", None, QtGui.QApplication.UnicodeUTF8))
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.formLayout_2.setWidget(2, QtGui.QFormLayout.LabelRole, self.label_13)
        self.param_sweep1_start_float = QtGui.QDoubleSpinBox(self.formLayoutWidget_2)
        self.param_sweep1_start_float.setObjectName(_fromUtf8("param_sweep1_start_float"))
        self.formLayout_2.setWidget(0, QtGui.QFormLayout.FieldRole, self.param_sweep1_start_float)
        self.param_sweep1_stop_float = QtGui.QDoubleSpinBox(self.formLayoutWidget_2)
        self.param_sweep1_stop_float.setObjectName(_fromUtf8("param_sweep1_stop_float"))
        self.formLayout_2.setWidget(1, QtGui.QFormLayout.FieldRole, self.param_sweep1_stop_float)
        self.param_sweep1_step_float = QtGui.QDoubleSpinBox(self.formLayoutWidget_2)
        self.param_sweep1_step_float.setObjectName(_fromUtf8("param_sweep1_step_float"))
        self.formLayout_2.setWidget(2, QtGui.QFormLayout.FieldRole, self.param_sweep1_step_float)
        self.sweep1_stackedWidget.addWidget(self.page_4)
        self.gridLayout.addWidget(self.sweep1_stackedWidget, 2, 0, 1, 2)
        self.verticalLayout.addWidget(self.sweep1_groupBox)
        self.param_sweep2_enabled = QtGui.QCheckBox(self.verticalLayoutWidget)
        self.param_sweep2_enabled.setText(QtGui.QApplication.translate("SweepDialog", "Enable Sweep 2", None, QtGui.QApplication.UnicodeUTF8))
        self.param_sweep2_enabled.setChecked(True)
        self.param_sweep2_enabled.setObjectName(_fromUtf8("param_sweep2_enabled"))
        self.verticalLayout.addWidget(self.param_sweep2_enabled)
        self.sweep2_groupBox = QtGui.QGroupBox(self.verticalLayoutWidget)
        self.sweep2_groupBox.setTitle(QtGui.QApplication.translate("SweepDialog", "Sweep 2", None, QtGui.QApplication.UnicodeUTF8))
        self.sweep2_groupBox.setObjectName(_fromUtf8("sweep2_groupBox"))
        self.gridLayoutWidget_4 = QtGui.QWidget(self.sweep2_groupBox)
        self.gridLayoutWidget_4.setGeometry(QtCore.QRect(9, 19, 361, 151))
        self.gridLayoutWidget_4.setObjectName(_fromUtf8("gridLayoutWidget_4"))
        self.gridLayout_4 = QtGui.QGridLayout(self.gridLayoutWidget_4)
        self.gridLayout_4.setMargin(0)
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.param_sweep2 = QtGui.QComboBox(self.gridLayoutWidget_4)
        self.param_sweep2.setObjectName(_fromUtf8("param_sweep2"))
        self.gridLayout_4.addWidget(self.param_sweep2, 0, 1, 1, 1)
        self.label_14 = QtGui.QLabel(self.gridLayoutWidget_4)
        self.label_14.setText(QtGui.QApplication.translate("SweepDialog", "Parameter", None, QtGui.QApplication.UnicodeUTF8))
        self.label_14.setObjectName(_fromUtf8("label_14"))
        self.gridLayout_4.addWidget(self.label_14, 0, 0, 1, 1)
        self.sweep2_stackedWidget = QtGui.QStackedWidget(self.gridLayoutWidget_4)
        self.sweep2_stackedWidget.setObjectName(_fromUtf8("sweep2_stackedWidget"))
        self.page_5 = QtGui.QWidget()
        self.page_5.setObjectName(_fromUtf8("page_5"))
        self.formLayoutWidget_3 = QtGui.QWidget(self.page_5)
        self.formLayoutWidget_3.setGeometry(QtCore.QRect(10, 20, 341, 91))
        self.formLayoutWidget_3.setObjectName(_fromUtf8("formLayoutWidget_3"))
        self.formLayout_3 = QtGui.QFormLayout(self.formLayoutWidget_3)
        self.formLayout_3.setMargin(0)
        self.formLayout_3.setObjectName(_fromUtf8("formLayout_3"))
        self.label_15 = QtGui.QLabel(self.formLayoutWidget_3)
        self.label_15.setText(QtGui.QApplication.translate("SweepDialog", "Start", None, QtGui.QApplication.UnicodeUTF8))
        self.label_15.setObjectName(_fromUtf8("label_15"))
        self.formLayout_3.setWidget(0, QtGui.QFormLayout.LabelRole, self.label_15)
        self.label_16 = QtGui.QLabel(self.formLayoutWidget_3)
        self.label_16.setText(QtGui.QApplication.translate("SweepDialog", "Stop", None, QtGui.QApplication.UnicodeUTF8))
        self.label_16.setObjectName(_fromUtf8("label_16"))
        self.formLayout_3.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_16)
        self.param_sweep2_start_int = QtGui.QSpinBox(self.formLayoutWidget_3)
        self.param_sweep2_start_int.setObjectName(_fromUtf8("param_sweep2_start_int"))
        self.formLayout_3.setWidget(0, QtGui.QFormLayout.FieldRole, self.param_sweep2_start_int)
        self.param_sweep2_stop_int = QtGui.QSpinBox(self.formLayoutWidget_3)
        self.param_sweep2_stop_int.setObjectName(_fromUtf8("param_sweep2_stop_int"))
        self.formLayout_3.setWidget(1, QtGui.QFormLayout.FieldRole, self.param_sweep2_stop_int)
        self.label_17 = QtGui.QLabel(self.formLayoutWidget_3)
        self.label_17.setText(QtGui.QApplication.translate("SweepDialog", "Step", None, QtGui.QApplication.UnicodeUTF8))
        self.label_17.setObjectName(_fromUtf8("label_17"))
        self.formLayout_3.setWidget(2, QtGui.QFormLayout.LabelRole, self.label_17)
        self.param_sweep2_step_int = QtGui.QSpinBox(self.formLayoutWidget_3)
        self.param_sweep2_step_int.setObjectName(_fromUtf8("param_sweep2_step_int"))
        self.formLayout_3.setWidget(2, QtGui.QFormLayout.FieldRole, self.param_sweep2_step_int)
        self.sweep2_stackedWidget.addWidget(self.page_5)
        self.page_6 = QtGui.QWidget()
        self.page_6.setObjectName(_fromUtf8("page_6"))
        self.formLayoutWidget_4 = QtGui.QWidget(self.page_6)
        self.formLayoutWidget_4.setGeometry(QtCore.QRect(10, 20, 341, 91))
        self.formLayoutWidget_4.setObjectName(_fromUtf8("formLayoutWidget_4"))
        self.formLayout_4 = QtGui.QFormLayout(self.formLayoutWidget_4)
        self.formLayout_4.setMargin(0)
        self.formLayout_4.setObjectName(_fromUtf8("formLayout_4"))
        self.label_18 = QtGui.QLabel(self.formLayoutWidget_4)
        self.label_18.setText(QtGui.QApplication.translate("SweepDialog", "Start", None, QtGui.QApplication.UnicodeUTF8))
        self.label_18.setObjectName(_fromUtf8("label_18"))
        self.formLayout_4.setWidget(0, QtGui.QFormLayout.LabelRole, self.label_18)
        self.label_19 = QtGui.QLabel(self.formLayoutWidget_4)
        self.label_19.setText(QtGui.QApplication.translate("SweepDialog", "Stop", None, QtGui.QApplication.UnicodeUTF8))
        self.label_19.setObjectName(_fromUtf8("label_19"))
        self.formLayout_4.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_19)
        self.label_20 = QtGui.QLabel(self.formLayoutWidget_4)
        self.label_20.setText(QtGui.QApplication.translate("SweepDialog", "Step", None, QtGui.QApplication.UnicodeUTF8))
        self.label_20.setObjectName(_fromUtf8("label_20"))
        self.formLayout_4.setWidget(2, QtGui.QFormLayout.LabelRole, self.label_20)
        self.param_sweep2_start_float = QtGui.QDoubleSpinBox(self.formLayoutWidget_4)
        self.param_sweep2_start_float.setObjectName(_fromUtf8("param_sweep2_start_float"))
        self.formLayout_4.setWidget(0, QtGui.QFormLayout.FieldRole, self.param_sweep2_start_float)
        self.param_sweep2_stop_float = QtGui.QDoubleSpinBox(self.formLayoutWidget_4)
        self.param_sweep2_stop_float.setObjectName(_fromUtf8("param_sweep2_stop_float"))
        self.formLayout_4.setWidget(1, QtGui.QFormLayout.FieldRole, self.param_sweep2_stop_float)
        self.param_sweep1_step_float_2 = QtGui.QDoubleSpinBox(self.formLayoutWidget_4)
        self.param_sweep1_step_float_2.setObjectName(_fromUtf8("param_sweep1_step_float_2"))
        self.formLayout_4.setWidget(2, QtGui.QFormLayout.FieldRole, self.param_sweep1_step_float_2)
        self.sweep2_stackedWidget.addWidget(self.page_6)
        self.gridLayout_4.addWidget(self.sweep2_stackedWidget, 2, 0, 1, 2)
        self.verticalLayout.addWidget(self.sweep2_groupBox)
        self.groupBox = QtGui.QGroupBox(self.verticalLayoutWidget)
        self.groupBox.setTitle(QtGui.QApplication.translate("SweepDialog", "Actions", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.gridLayoutWidget_2 = QtGui.QWidget(self.groupBox)
        self.gridLayoutWidget_2.setGeometry(QtCore.QRect(10, 20, 379, 124))
        self.gridLayoutWidget_2.setObjectName(_fromUtf8("gridLayoutWidget_2"))
        self.gridLayout_2 = QtGui.QGridLayout(self.gridLayoutWidget_2)
        self.gridLayout_2.setMargin(0)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.addActionButton = QtGui.QPushButton(self.gridLayoutWidget_2)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.addActionButton.sizePolicy().hasHeightForWidth())
        self.addActionButton.setSizePolicy(sizePolicy)
        self.addActionButton.setText(QtGui.QApplication.translate("SweepDialog", "add", None, QtGui.QApplication.UnicodeUTF8))
        self.addActionButton.setObjectName(_fromUtf8("addActionButton"))
        self.gridLayout_2.addWidget(self.addActionButton, 0, 1, 1, 1)
        self.actions_comboBox = QtGui.QComboBox(self.gridLayoutWidget_2)
        self.actions_comboBox.setObjectName(_fromUtf8("actions_comboBox"))
        self.gridLayout_2.addWidget(self.actions_comboBox, 0, 0, 1, 1)
        self.clearActionsButton = QtGui.QPushButton(self.gridLayoutWidget_2)
        self.clearActionsButton.setText(QtGui.QApplication.translate("SweepDialog", "Clear", None, QtGui.QApplication.UnicodeUTF8))
        self.clearActionsButton.setObjectName(_fromUtf8("clearActionsButton"))
        self.gridLayout_2.addWidget(self.clearActionsButton, 1, 1, 1, 1)
        self.param_actions = QtGui.QTextEdit(self.gridLayoutWidget_2)
        self.param_actions.setObjectName(_fromUtf8("param_actions"))
        self.gridLayout_2.addWidget(self.param_actions, 1, 0, 1, 1)
        self.verticalLayout.addWidget(self.groupBox)
        self.buttonBox = QtGui.QDialogButtonBox(SweepDialog)
        self.buttonBox.setGeometry(QtCore.QRect(50, 570, 341, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))

        self.retranslateUi(SweepDialog)
        self.sweep1_stackedWidget.setCurrentIndex(1)
        self.sweep2_stackedWidget.setCurrentIndex(0)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), SweepDialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), SweepDialog.reject)
        QtCore.QObject.connect(self.param_sweep2_enabled, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.sweep2_groupBox.setVisible)
        QtCore.QMetaObject.connectSlotsByName(SweepDialog)

    def retranslateUi(self, SweepDialog):
        pass

