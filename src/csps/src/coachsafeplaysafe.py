from kivy.lang import Builder
from kivy.properties import ListProperty
from kivy.app import App
import threespace_api as ts_api
import os

no_sensor_api = False
no_sensor = False
notEstablished = True

if os.name == 'posix':
	no_sensor_api = True

KV = '''
#:import Camera kivy.uix.camera
#:import chain itertools.chain
#:import Clock kivy.clock.Clock

<Graph@Widget>:
	max: 100
	canvas:
		Color:
			rgba: .4, .4, 1, 1
		Line:
			points:
				list(chain(*
				[[
				self.x + x * self.width / len(app.gyroX),
				self.y + (self.height * 0.3) + y * self.height / (self.max * 100)
				] for x, y in enumerate(app.gyroX)])) if app.gyroX else []

		Color:
			rgba: 1, .4, .4, 1
		Line:
			points:
				list(chain(*
				[[
				self.x + x * self.width / len(app.gyroY),
				self.y + (self.height * 0.3) + y * self.height / (self.max * 100)
				] for x, y in enumerate(app.gyroY)])) if app.gyroY else []

		Color:
			rgba: .4, 1, .4, 1
		Line:
			points:
				list(chain(*
				[[
				self.x + x * self.width / len(app.gyroZ),
				self.y + (self.height * 0.3) + y * self.height / (self.max * 100)
				] for x, y in enumerate(app.gyroZ)])) if app.gyroZ else []

		Color:
			rgba: .4, .4, 1, 1
		Line:
			points:
				list(chain(*
				[[
				self.x + x * self.width / len(app.accelX),
				self.y + (self.height * 0.7) + y * self.height / (self.max * 100)
				] for x, y in enumerate(app.accelX)])) if app.accelX else []

		Color:
			rgba: 1, .4, .4, 1
		Line:
			points:
				list(chain(*
				[[
				self.x + x * self.width / len(app.accelY),
				self.y + (self.height * 0.7) + y * self.height / (self.max * 100)
				] for x, y in enumerate(app.accelY)])) if app.accelY else []

		Color:
			rgba: .4, 1, .4, 1
		Line:
			points:
				list(chain(*
				[[
				self.x + x * self.width / len(app.accelZ),
				self.y + (self.height * 0.7) + y * self.height / (self.max * 100)
				] for x, y in enumerate(app.accelZ)])) if app.accelZ else []

GridLayout:
	rows: 2
	GridLayout:
		cols: 2
		Camera:
			id: 'Cam'
			index: 0
			resolution: (1920, 1080)
			size: (self.width, self.height)
		Graph:
			max: 1

	ToggleButton:
		text: 'Start/Stop'
		size_hint_y: None
		on_state:
			if self.state == 'down': Clock.schedule_interval(app.add_running_values, 0)
			else: Clock.unschedule(app.add_running_values)
'''


class Graph(App):
	gyroX = ListProperty([])
	gyroY = ListProperty([])
	gyroZ = ListProperty([])

	accelX = ListProperty([])
	accelY = ListProperty([])
	accelZ = ListProperty([])

	def __init__(self, **kwargs):
		super(Graph, self).__init__(**kwargs)
		global no_sensor
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
			print ("No sensor API")

	def build(self):
		return Builder.load_string(KV)

	def add_running_values(self, dt):
		global notEstablished

		if not no_sensor:
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
			if notEstablished:
				print ("No sensor")
				notEstablished = False


if __name__ == '__main__':
	Graph().run()
