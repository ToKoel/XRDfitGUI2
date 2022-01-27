import sys
import matplotlib
from PyQt5 import QtCore
from PyQt5 import QtWidgets as QW
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from fit import Data, Fit
import numpy as np
matplotlib.use('Qt5Agg')


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.figure = fig
        self.axes = fig.add_subplot(111)
        self.axes.set_xlabel(r"Q [$\AA^{-1}$]", fontsize=15)
        self.axes.set_ylabel("I [a.u.]", fontsize=15)
        self.axes.axhline(0, color=(0, 0, 0), lw=1)
        super(MplCanvas, self).__init__(fig)


class Window(QW.QDialog):

    def __init__(self, data_obj, parent=None):
        super(Window, self).__init__(parent)
        self.setWindowTitle("XRDfitGUI")
        self.canvas = MplCanvas(self, width=8, height=6, dpi=100)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.Data = data_obj
        self.y_data_original = self.Data.intensity
        self.x_data_original = self.Data.Q
        self.y_data = self.y_data_original
        self.x_data = self.x_data_original
        self.Fit = Fit(self.x_data, self.y_data)
        self.background_points = []

        self.textbox_D = QW.QLineEdit(self)
        self.textbox_D.setText("153.6")
        self.label_D = QW.QLabel(self)
        self.label_D.setText("D [A]: ")
        self.vary_D = QW.QCheckBox(self)

        self.textbox_delta = QW.QLineEdit(self)
        self.textbox_delta.setText("0.05")
        self.label_delta = QW.QLabel(self)
        self.label_delta.setText("delta: ")
        self.vary_delta = QW.QCheckBox(self)

        self.textbox_I1 = QW.QLineEdit(self)
        self.textbox_I1.setText("0.034")
        self.label_I1 = QW.QLabel(self)
        self.label_I1.setText("I1: ")
        self.textbox_P1 = QW.QLineEdit(self)
        self.textbox_P1.setText("2.126")
        self.label_P1 = QW.QLabel(self)
        self.label_P1.setText("P1: ")
        self.vary_P1 = QW.QCheckBox(self)
        self.vary_I1 = QW.QCheckBox(self)

        self.textbox_I2 = QW.QLineEdit(self)
        self.textbox_I2.setText("0.099")
        self.label_I2 = QW.QLabel(self)
        self.label_I2.setText("I2: ")
        self.textbox_P2 = QW.QLineEdit(self)
        self.textbox_P2.setText("2.488")
        self.label_P2 = QW.QLabel(self)
        self.label_P2.setText("P2: ")
        self.vary_P2 = QW.QCheckBox(self)
        self.vary_I2 = QW.QCheckBox(self)

        self.textbox_I3 = QW.QLineEdit(self)
        self.textbox_I3.setText("0.002")
        self.label_I3 = QW.QLabel(self)
        self.label_I3.setText("I3: ")
        self.textbox_P3 = QW.QLineEdit(self)
        self.textbox_P3.setText("2.598")
        self.label_P3 = QW.QLabel(self)
        self.label_P3.setText("P3: ")
        self.vary_P3 = QW.QCheckBox(self)
        self.vary_I3 = QW.QCheckBox(self)

        self.textbox_I4 = QW.QLineEdit(self)
        self.textbox_I4.setText("0.017")
        self.label_I4 = QW.QLabel(self)
        self.label_I4.setText("I4: ")
        self.textbox_P4 = QW.QLineEdit(self)
        self.textbox_P4.setText("3.006")
        self.label_P4 = QW.QLabel(self)
        self.label_P4.setText("P4: ")
        self.vary_P4 = QW.QCheckBox(self)
        self.vary_I4 = QW.QCheckBox(self)

        self.textbox_I5 = QW.QLineEdit(self)
        self.textbox_I5.setText("0.0")
        self.label_I5 = QW.QLabel(self)
        self.label_I5.setText("I5: ")
        self.textbox_P5 = QW.QLineEdit(self)
        self.textbox_P5.setText("3.216")
        self.label_P5 = QW.QLabel(self)
        self.label_P5.setText("P5: ")
        self.vary_P5 = QW.QCheckBox(self)
        self.vary_I5 = QW.QCheckBox(self)

        self.textbox_I6 = QW.QLineEdit(self)
        self.textbox_I6.setText("0.012")
        self.label_I6 = QW.QLabel(self)
        self.label_I6.setText("I6: ")
        self.textbox_P6 = QW.QLineEdit(self)
        self.textbox_P6.setText("3.682")
        self.label_P6 = QW.QLabel(self)
        self.label_P6.setText("P6: ")
        self.vary_P6 = QW.QCheckBox(self)
        self.vary_I6 = QW.QCheckBox(self)

        self.textbox_I7 = QW.QLineEdit(self)
        self.textbox_I7.setText("0.031")
        self.label_I7 = QW.QLabel(self)
        self.label_I7.setText("I7: ")
        self.textbox_P7 = QW.QLineEdit(self)
        self.textbox_P7.setText("3.908")
        self.label_P7 = QW.QLabel(self)
        self.label_P7.setText("P7: ")
        self.vary_P7 = QW.QCheckBox(self)
        self.vary_I7 = QW.QCheckBox(self)

        self.textbox_I8 = QW.QLineEdit(self)
        self.textbox_I8.setText("0.0")
        self.label_I8 = QW.QLabel(self)
        self.label_I8.setText("I8: ")
        self.textbox_P8 = QW.QLineEdit(self)
        self.textbox_P8.setText("3.908")
        self.label_P8 = QW.QLabel(self)
        self.label_P8.setText("P8: ")
        self.vary_P8 = QW.QCheckBox(self)
        self.vary_I8 = QW.QCheckBox(self)

        self.textbox_I9 = QW.QLineEdit(self)
        self.textbox_I9.setText("0.046")
        self.label_I9 = QW.QLabel(self)
        self.label_I9.setText("I9: ")
        self.textbox_P9 = QW.QLineEdit(self)
        self.textbox_P9.setText("4.252")
        self.label_P9 = QW.QLabel(self)
        self.label_P9.setText("P9: ")
        self.vary_P9 = QW.QCheckBox(self)
        self.vary_I9 = QW.QCheckBox(self)

        # self.textbox_bkgm = QW.QLineEdit(self)
        # self.textbox_bkgm.setText("0.02339")
        # self.label_bkgm = QW.QLabel(self)
        # self.label_bkgm.setText("bkg m: ")
        # self.textbox_bkgt = QW.QLineEdit(self)
        # self.textbox_bkgt.setText("-0.03875")
        # self.label_bkgt = QW.QLabel(self)
        # self.label_bkgt.setText("bkg t: ")
        # self.vary_bkgt = QW.QCheckBox(self)
        # self.vary_bkgm = QW.QCheckBox(self)

        self.textbox_eta = QW.QLineEdit(self)
        self.textbox_eta.setText("0.0")
        self.label_eta = QW.QLabel(self)
        self.label_eta.setText("eta: ")
        self.textbox_sig = QW.QLineEdit(self)
        self.textbox_sig.setText("0.0")
        self.label_sig = QW.QLabel(self)
        self.label_sig.setText("sig: ")
        self.vary_sig = QW.QCheckBox(self)
        self.vary_eta = QW.QCheckBox(self)

        parameters = [[self.label_I1, self.textbox_I1,
                       self.vary_I1,
                       self.label_P1, self.textbox_P1,
                       self.vary_P1],
                      [self.label_I2, self.textbox_I2,
                       self.vary_I2,
                       self.label_P2, self.textbox_P2,
                       self.vary_P2],
                      [self.label_I3, self.textbox_I3,
                       self.vary_I3,
                       self.label_P3, self.textbox_P3,
                       self.vary_P3],
                      [self.label_I4, self.textbox_I4,
                       self.vary_I4,
                       self.label_P4, self.textbox_P4,
                       self.vary_P4],
                      [self.label_I5, self.textbox_I5,
                       self.vary_I5,
                       self.label_P5, self.textbox_P5,
                       self.vary_P5],
                      [self.label_I6, self.textbox_I6,
                       self.vary_I6,
                       self.label_P6, self.textbox_P6,
                       self.vary_P6],
                      [self.label_I7, self.textbox_I7,
                       self.vary_I7,
                       self.label_P7, self.textbox_P7,
                       self.vary_P7],
                      [self.label_I8, self.textbox_I8,
                       self.vary_I8,
                       self.label_P8, self.textbox_P8,
                       self.vary_P8],
                      [self.label_I9, self.textbox_I9,
                       self.vary_I9,
                       self.label_P9, self.textbox_P9,
                       self.vary_P9],
                     # [self.label_bkgm, self.textbox_bkgm,
                     #  self.vary_bkgm,
                     #  self.label_bkgt, self.textbox_bkgt,
                     #  self.vary_bkgt],
                      [self.label_eta, self.textbox_eta,
                       self.vary_eta,
                       self.label_sig, self.textbox_sig,
                       self.vary_sig]
                      ]

        self.button = QW.QPushButton(self)
        self.button.setText("OK")
        self.button.clicked.connect(self.on_clicked)

        self.reset_button = QW.QPushButton(self)
        self.reset_button.setText("reset")
        self.reset_button.clicked.connect(self.reset)

        self.bkg_button = QW.QPushButton(self)
        self.bkg_button.setText("select bkg points")
        self.bkg_button.clicked.connect(self.bkg_clicked)
        self.bkg_button2 = QW.QPushButton(self)
        self.bkg_button2.setText("fit bkg")
        self.textbox_deg = QW.QLineEdit(self)
        self.textbox_deg.setText("3")
        self.label_deg = QW.QLabel(self)
        self.label_deg.setText("Cheb. deg.: ")
        self.bkg_button2.clicked.connect(self.bkg_clicked_finish)

        self.textbox_qmin = QW.QLineEdit(self)
        self.textbox_qmin.setText("1.8")
        self.label_qmin = QW.QLabel(self)
        self.label_qmin.setText("Q_min [A^(-1)]: ")
        self.textbox_qmax = QW.QLineEdit(self)
        self.textbox_qmax.setText("4.5")
        self.label_qmax = QW.QLabel(self)
        self.label_qmax.setText("Q_max [A^(-1)]: ")

        self.fit_button = QW.QPushButton(self)
        self.fit_button.setText("Fit")
        self.fit_button.clicked.connect(self.on_clicked_fit)

        self.simulate_button = QW.QPushButton(self)
        self.simulate_button.setText("Simulate")
        self.simulate_button.clicked.connect(lambda: self.on_clicked_fit(True))

        # layout
        layout = QW.QGridLayout(self)

        layout.addWidget(self.toolbar, 0, 0)
        layout.addWidget(self.canvas, 1, 0)

        q_range_layout = QW.QHBoxLayout()
        q_range_layout.addWidget(self.button)
        q_range_layout.addWidget(self.reset_button)
        q_range_layout.addWidget(self.label_qmin)
        q_range_layout.addWidget(self.textbox_qmin)
        q_range_layout.addWidget(self.label_qmax)
        q_range_layout.addWidget(self.textbox_qmax)

        bkg_layout = QW.QVBoxLayout()
        bkg_layout.addWidget(self.bkg_button)
        bkg_fit_layout = QW.QHBoxLayout()
        bkg_fit_layout.addWidget(self.bkg_button2)
        bkg_fit_layout.addWidget(self.label_deg)
        bkg_fit_layout.addWidget(self.textbox_deg)
        bkg_layout.addLayout(bkg_fit_layout)
        q_range_layout.addLayout(bkg_layout)

        layout.addLayout(q_range_layout, 2, 0)

        parameters_layout = QW.QVBoxLayout()
        par1_layout = QW.QHBoxLayout()
        par1_layout.addWidget(self.label_D)
        par1_layout.addWidget(self.textbox_D)
        par1_layout.addWidget(self.vary_D)
        parameters_layout.addLayout(par1_layout)

        par2_layout = QW.QHBoxLayout()
        par2_layout.addWidget(self.label_delta)
        par2_layout.addWidget(self.textbox_delta)
        par2_layout.addWidget(self.vary_delta)
        parameters_layout.addLayout(par2_layout)

        for par in parameters:
            par_layout = QW.QHBoxLayout()
            par_layout.addWidget(par[0])
            par_layout.addWidget(par[1])
            par_layout.addWidget(par[2])
            par_layout.addWidget(par[3])
            par_layout.addWidget(par[4])
            par_layout.addWidget(par[5])
            parameters_layout.addLayout(par_layout)

        layout.addLayout(parameters_layout, 1, 1)
        fit_button_layout = QW.QHBoxLayout()
        fit_button_layout.addWidget(self.fit_button)
        fit_button_layout.addWidget(self.simulate_button)
        layout.addLayout(fit_button_layout, 2, 1)
        layout.setColumnStretch(0, 2)
        layout.setColumnStretch(1, 1)

        # We need to store a reference to the plotted line
        # somewhere, so we can apply the new data to it.
        self._plot_ref = None
        self._plot_ref2 = None
        self._plot_ref3 = None
        self._plot_ref4 = None
        self._plot_ref5 = None
        self.update_plot()

        self.show()

    def update_plot(self):
        if self._plot_ref is None and self._plot_ref2 is None:
            plot_refs = self.canvas.axes.plot(self.x_data, self.y_data, marker='o', lw=0,
                                              color=(0.1804, 0.4235, 0.72157),
                                              markersize=3, markerfacecolor=(1, 1, 1),
                                              label="Data")
            plot_refs2 = self.canvas.axes.plot(self.x_data, self.Fit.fit_curve, lw=1.2,
                                               color=(0.988, 0.2667, 0.2706),
                                               label="Fit")
            plot_refs3 = self.canvas.axes.plot(self.x_data, np.zeros_like(self.x_data), lw=1.2,
                                               color=(0.0, 0.8, 0.0),
                                               label="Difference")

            self.canvas.axes.set_xlim(self.x_data.min(), self.x_data.max())
            self.canvas.axes.set_ylim(-self.y_data.max() * 0.1,
                                      self.y_data.max() + 0.1 * self.y_data.max())

            self.canvas.axes.legend(prop={'size': 13})
            self._plot_ref = plot_refs[0]
            self._plot_ref2 = plot_refs2[0]
            self._plot_ref3 = plot_refs3[0]
        else:
            self.canvas.axes.set_xlim(self.x_data.min(), self.x_data.max())
            self.canvas.axes.set_ylim(-self.y_data.max() * 0.1,
                                      self.y_data.max() + 0.1 * self.y_data.max())
            self._plot_ref.set_xdata(self.x_data)
            self._plot_ref.set_ydata(self.y_data)
            self._plot_ref2.set_xdata(self.x_data)
            self._plot_ref2.set_ydata(self.Fit.fit_curve)
            if self._plot_ref3 is not None:
                self._plot_ref3.set_xdata(self.x_data)
                self._plot_ref3.set_ydata(self.y_data - self.Fit.fit_curve)

        # Trigger the canvas to update and redraw.
        self.canvas.draw()

    def on_clicked(self):
        input_number_qmin = self.textbox_qmin.text()
        input_number_qmax = self.textbox_qmax.text()
        self.q_range(float(input_number_qmin),
                     float(input_number_qmax))
        self.Fit = Fit(self.x_data, self.y_data)
        self.update_plot()

    def on_clicked_fit(self, simulate=False):
        params = [(self.textbox_D.text(), self.vary_D.isChecked()),
                  (self.textbox_delta.text(), self.vary_delta.isChecked()),
                  (self.textbox_I1.text(), self.vary_I1.isChecked()),
                  (self.textbox_P1.text(), self.vary_P1.isChecked()),
                  (self.textbox_I2.text(), self.vary_I2.isChecked()),
                  (self.textbox_P2.text(), self.vary_P2.isChecked()),
                  (self.textbox_I3.text(), self.vary_I3.isChecked()),
                  (self.textbox_P3.text(), self.vary_P3.isChecked()),
                  (self.textbox_I4.text(), self.vary_I4.isChecked()),
                  (self.textbox_P4.text(), self.vary_P4.isChecked()),
                  (self.textbox_I5.text(), self.vary_I5.isChecked()),
                  (self.textbox_P5.text(), self.vary_P5.isChecked()),
                  (self.textbox_I6.text(), self.vary_I6.isChecked()),
                  (self.textbox_P6.text(), self.vary_P6.isChecked()),
                  (self.textbox_I7.text(), self.vary_I7.isChecked()),
                  (self.textbox_P7.text(), self.vary_P7.isChecked()),
                  (self.textbox_I8.text(), self.vary_I8.isChecked()),
                  (self.textbox_P8.text(), self.vary_P8.isChecked()),
                  (self.textbox_I9.text(), self.vary_I9.isChecked()),
                  (self.textbox_P9.text(), self.vary_P9.isChecked()),
                  (self.textbox_eta.text(), self.vary_eta.isChecked()),
                  (self.textbox_sig.text(), self.vary_sig.isChecked())]
            # ,
             #     (self.textbox_bkgm.text(), self.vary_bkgm.isChecked()),
             #     (self.textbox_bkgt.text(), self.vary_bkgt.isChecked())]
        param_dict = {'D': (0, False), 'delta': (0, False),
                      'I1': (0, False), 'P1': (0, False),
                      'I2': (0, False), 'P2': (0, False),
                      'I3': (0, False), 'P3': (0, False),
                      'I4': (0, False), 'P4': (0, False),
                      'I5': (0, False), 'P5': (0, False),
                      'I6': (0, False), 'P6': (0, False),
                      'I7': (0, False), 'P7': (0, False),
                      'I8': (0, False), 'P8': (0, False),
                      'I9': (0, False), 'P9': (0, False),
                      'eta': (0, False), 'sig': (0, False)} #,
                    #  'bkgm': (0, False), 'bkgt': (-0.5, False)}
        for key, par in zip(param_dict.keys(), params):
            if par[0] == '' or par[1] == '':
                msg = QW.QMessageBox()
                msg.setIcon(QW.QMessageBox.Critical)
                msg.setText("Error")
                msg.setInformativeText("Parameter is empty")
                msg.exec_()
                break
            else:
                param_dict[key] = (float(par[0]), par[1])

        self.Fit.update_parameters(param_dict)

        if simulate:
            self.Fit.simulate()
            self.update_plot()
        else:
            count = 0
            for par in params:
                if par[1] == False:
                    count += 1
            if count == len(params):
                msg = QW.QMessageBox()
                msg.setIcon(QW.QMessageBox.Critical)
                msg.setText("Error")
                msg.setInformativeText("Nothing to refine")
                msg.exec_()
                return

            self.Fit.fit()
            self.update_plot()

            parameters = [self.textbox_D, self.textbox_delta,
                          self.textbox_I1, self.textbox_P1,
                          self.textbox_I2, self.textbox_P2,
                          self.textbox_I3, self.textbox_P3,
                          self.textbox_I4, self.textbox_P4,
                          self.textbox_I5, self.textbox_P5,
                          self.textbox_I6, self.textbox_P6,
                          self.textbox_I7, self.textbox_P7,
                          self.textbox_I8, self.textbox_P8,
                          self.textbox_I9, self.textbox_P9,
                          self.textbox_eta, self.textbox_sig] #,
                          # self.textbox_bkgm, self.textbox_bkgt]

            for par, res in zip(parameters, self.Fit.result_params.valuesdict().values()):
                par.setText("%.5f" % res)

    def reset(self):
        self.x_data = self.x_data_original
        self.y_data = self.y_data_original
        self.Fit = Fit(self.x_data, self.y_data)
        self.update_plot()

    def bkg_clicked(self):
        if self._plot_ref4 is None:
            crosses, = self.canvas.axes.plot([], [], marker='x', lw=0, color=(0.9, 0, 0))
            self._plot_ref4 = crosses
        else:
            self._plot_ref4.set_xdata([])
            self._plot_ref4.set_ydata([])
            self.background_points.clear()
        bkg_crosses = BackgroundCrosses(self._plot_ref4)
        self.background_points.append((bkg_crosses.xs, bkg_crosses.ys))

    def bkg_clicked_finish(self):
        print(self.background_points[0][0])
        x = np.array(self.background_points[0][0])
        y = np.array(self.background_points[0][1])
        cheb_coeffs = np.polynomial.chebyshev.Chebyshev.fit(x, y, int(self.textbox_deg.text()))
        if self._plot_ref5 is None:
            bkg_fit, = self.canvas.axes.plot(self.x_data, cheb_coeffs(self.x_data), lw =1, color=(0.9,0,0))
            self._plot_ref5 = bkg_fit
        else:
            self._plot_ref5.set_xdata(self.x_data)
            self._plot_ref5.set_ydata(cheb_coeffs(self.x_data))
        self.update_plot()
        self.Fit.bkg_poly = cheb_coeffs(self.x_data)



    def q_range(self, qmin, qmax, step=1):
        idx_l = (np.abs(self.x_data - qmin)).argmin()
        idx_u = (np.abs(self.x_data - qmax)).argmin()
        self.y_data = self.y_data[idx_l:idx_u:step]
        self.x_data = self.x_data[idx_l:idx_u:step]
        print("Q-range adjusted")


class BackgroundCrosses:

    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        if event.inaxes != self.line.axes:
            return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()


class XRDfitGUI:

    def __init__(self):
        app = QW.QApplication(sys.argv)

        file = "OC15b2_G50_0p9mm_E28p7keV_T14350eV_M_WAXS.xye"
        bkg = "MT_G50_0p9mm_E28p7keV_T14350eV_M_WAXS.xye"
        data = Data(file, 0.4329, bkg)

        w = Window(data)
        app.exec_()
