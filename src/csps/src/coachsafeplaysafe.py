from kivy.app import App
from kivy.uix.screenmanager import ScreenManager, Screen, SwapTransition
from kivy.properties import ListProperty
import threespace_api as ts_api
from random import gauss
import os

no_sensor_api = False
no_sensor = False
notEstablished = True

if os.name == 'posix':
    no_sensor_api = True


class ScreenOne(Screen):
    pass


class ScreenTwo(Screen):
    gyroX = ListProperty([])
    gyroY = ListProperty([])
    gyroZ = ListProperty([])

    accelX = ListProperty([])
    accelY = ListProperty([])
    accelZ = ListProperty([])

    values = ListProperty([])

    def __init__(self, **kwargs):
        super(ScreenTwo, self).__init__(**kwargs)
        global no_sensor, no_sensor_api
        self.tssensor = None
        if not no_sensor_api:
            self.com_port = None
            self.device_list = ts_api.getComPorts()
            if len(self.device_list) > 0:
                for self.device in self.device_list:
                    self.cp = self.device.com_port
                    self.port_info = ts_api.getDeviceInfoFromComPort(self.cp, poll_device=True)
                    if self.port_info.dev_type == 'BT-H3' or self.port_info.dev_type == 'BT':
                        self.com_port = self.cp
                        self.tssensor = ts_api.TSBTSensor(self.com_port)
                        break
                    if self.port_info.dev_type == 'DNG':
                        self.com_port = self.cp
                        self.tsdongle = ts_api.TSDongle(self.com_port)
                        if not self.tsdongle:
                            no_sensor = True
                            break
                        self.tssensor = self.tsdongle.getSensorFromDongle(0)
                        if not self.tssensor:
                            no_sensor = True
                            break
                        break
                    if self.port_info.dev_type == 'WL-H3' or self.port_info.dev_type == 'WL':
                        self.com_port = self.cp
                        self.tssensor = ts_api.TSWLSensor(self.com_port)
                        break
                if not self.com_port:
                    no_sensor = True
            else:
                no_sensor = True
        else:
            no_sensor = True
            no_sensor_api = True
            print ("No Sensor API")

    def add_running_values(self, dt):
        global notEstablished, no_sensor

        if not no_sensor:
            print ("zxc")
            data = self.tssensor.getCorrectedGyroRate()
            self.gyroX.append(data[0])
            self.gyroX = self.gyroX[-100:]
            self.gyroY.append(data[1])
            self.gyroY = self.gyroY[-100:]
            self.gyroZ.append(data[2])
            self.gyroZ = self.gyroZ[-100:]

            data = self.tssensor.getCorrectedAccelerometerVector()
            self.accelX.append(data[0])
            self.accelX = self.accelX[-100:]
            self.accelY.append(data[1])
            self.accelY = self.accelY[-100:]
            self.accelZ.append(data[2])
            self.accelZ = self.accelZ[-100:]
        else:
            self.values.append(gauss(.5, .1))
            self.values = self.values[-100:]
            if notEstablished:
                print ("Sensor not Connected")
                notEstablished = False


class CoachSafePlaySafeApp(App):
    def build(self):
        screen_manager = ScreenManager(transition=SwapTransition())
        screen_manager.add_widget(ScreenOne(name="MPV"))
        screen_manager.add_widget(ScreenTwo(name="VidSensor"))
        return screen_manager


if __name__ == "__main__":
    CoachSafePlaySafeApp().run()
