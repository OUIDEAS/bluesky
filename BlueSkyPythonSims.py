import bluesky as bs
import numpy as np
# from bluesky.simulation import ScreenIO
import matplotlib.pyplot as plt
from bluesky.network.client import Client
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QTextEdit

# class ScreenDummy(ScreenIO):
#     """
#     Dummy class for the screen. Inherits from ScreenIO to make sure all the
#     necessary methods are there. This class is there to reimplement the echo
#     method so that console messages are printed.
#     """
#     def echo(self, text='', flags=0):
#         """Just print echo messages"""
#         print("BlueSky console:", text)




# The echo textbox, command line, and bluesky network client as globals
echobox = None
cmdline = None
bsclient = None


class TextClient(Client):
    '''
        Subclassed Client with a timer to periodically check for incoming data,
        an overridden event function to handle data, and a stack function to
        send stack commands to BlueSky.
    '''
    def __init__(self):
        super().__init__()
        self.timer = QTimer()
        self.timer.timeout.connect(self.update)
        self.timer.start(20)

    def event(self, name, data, sender_id):
        ''' Overridden event function to handle incoming ECHO commands. '''
        if name == b'ECHO' and echobox is not None:
            echobox.echo(**data)

    def stack(self, text):
        ''' Stack function to send stack commands to BlueSky. '''
        self.send_event(b'STACK', text)

    def echo(self, text, flags=None):
        ''' Overload Client's echo function. '''
        if echobox is not None:
            echobox.echo(text, flags)
    def get_trajectories(self):
        return self.send_event(b'GET_TRAJECTORIES')

class Echobox(QTextEdit):
    ''' Text box to show echoed text coming from BlueSky. '''
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(150)
        self.setReadOnly(True)
        # self.setFocusPolicy(Qt.NoFocus)

    def echo(self, text, flags=None):
        ''' Add text to this echo box. '''
        self.append(text)
        self.verticalScrollBar().setValue(self.verticalScrollBar().maximum())


# class Cmdline(QTextEdit):
#     ''' Wrapper class for the command line. '''
#     def __init__(self, parent=None):
#         super().__init__(parent)
#         self.setMaximumHeight(21)
#         # self.setFocusPolicy(Qt.StrongFocus)
#         # self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

#     def keyPressEvent(self, event):
#         ''' Handle Enter keypress to send a command to BlueSky. '''
#         if event.key() == Qt.Key.Key_Enter or event.key() == Qt.Key.Key_Return:
#             if bsclient is not None:
#                 bsclient.stack(self.toPlainText())
#                 echobox.echo(self.toPlainText())
#             self.setText('')
#         else:
#             super().keyPressEvent(event)





#Creating Aircraft
# bs.stack.stack('CRE KL001 B737 KCMH 0 FL150 300')
# bs.sim.step()
# print('Lat Lon', bs.traf.lat, bs.traf.lon)

if __name__ == '__main__':
    # Construct the Qt main object
    app = QApplication([])

    # Create a window with a stack text box and a command line
    win = QWidget()
    win.setWindowTitle('External Client with a GUI')
    layout = QVBoxLayout()
    win.setLayout(layout)

    echobox = Echobox(win)
    # cmdline = Cmdline(win)
    layout.addWidget(echobox)
    layout.addWidget(cmdline)
    # win.show()

    # Create and start BlueSky client
    bsclient = TextClient()
    bsclient.connect(event_port=11000, stream_port=11001)
    # app.exec()
    bs.init(mode ='sim')
    # bs.scr = ScreenDummy()

    n = 3

    bs.traf.mcre(n, actype ="F150")
    bsclient.get_trajectories()

    for acid in bs.traf.id:
        # set the origin (not needed if initialized in flight),
        # and add some waypoints, here only the altitude (in m) is passed to the
        # function, but you can additionally pass a speed as well
        # finally turn on VNAV for each flight
        bsclient.stack(f'ORIG {acid} KCMH;'
                    f'ADDWPT {acid} APE FL150;'
                    f'ADDWPT {acid} BUD FL150;'
                    f'ADDWPT {acid} TVT FL150;'
                    f'ADDWPT {acid} JPU FL150;'
                    f'VNAV {acid} ON')
    bs.stack.stack('DT 1;FF')

    # we'll run the simulation for up to 4000 seconds
    t_max = 5000

    ntraf = bs.traf.ntraf
    n_steps = int(t_max + 1)
    t = np.linspace(0, t_max, n_steps)

    # allocate some empty arrays for the results
    res = np.zeros((n_steps, 4, ntraf))

    # iteratively simulate the traffic
    for i in range(n_steps):
        # Perform one step of the simulation
        # bs.sim.step()

        # save the results from the simulator in the results array,
        # here we keep the latitude, longitude, altitude and TAS
        res[i] = [bs.traf.lat,
                    bs.traf.lon,
                    bs.traf.alt,
                    bs.traf.tas]

    # for idx, acid in enumerate(bs.traf.id):
    #     fig = plt.figure(figsize=(10, 15))
    #     ax1 = plt.subplot2grid((4, 1), (0, 0), rowspan=2)
    #     ax2 = plt.subplot2grid((4, 1), (2, 0))
    #     ax3 = plt.subplot2grid((4, 1), (3, 0))

    #     ax1.plot(res[:, 1, idx], res[:, 0, idx])
    #     ax1.set_xlabel('lon')
    #     ax1.set_ylabel('lat')

    #     ax2.plot(t, res[:, 2, idx])
    #     ax2.set_xlabel('t [s]')
    #     ax2.set_ylabel('alt [m]')

    #     ax3.plot(t, res[:, 3, idx])
    #     ax3.set_xlabel('t [s]')
    #     ax3.set_ylabel('TAS [m/s]')
        
    #     fig.suptitle(f'Trajectory {acid}')

    # plt.show()