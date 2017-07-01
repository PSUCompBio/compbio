from __future__ import print_function, division, absolute_import
import threespace_api as ts_api
from PIL import Image, ImageTk, ImageFont, ImageDraw, ImageGrab
import matplotlib.pyplot as plt
from serial import Serial
import Tkinter as tk
import tkMessageBox
import numpy as np
import time
import sys
import cv2
import math
import random
import timeit
import tinys3
import smtplib
import os

no_sensor_api = False

if os.name == 'posix':
    no_sensor_api = True

# constants for framerates replay duration etc - Ishan
CONST_videoRateMs = 20
CONST_replayDuration = 1
CONST_processingSlack = 0.7
CONST_slowDown = 1.0
CONST_cacheLimit = 1000 / CONST_videoRateMs * CONST_replayDuration * CONST_processingSlack
CONST_panelWidth = 300
CONST_panelHeight = 300
CONST_videoWidth = 382
CONST_videoHeight = 866
CONST_ImpactGroupTimeThreshold = 1

# flags to check for replay and recordings - Ishan
replay_on = False
record_on = False

# videocache and replay cache and current frame in replay - Ishan
videoCache = []
replayCache = []
timeCache = []
replayTimeCache = []
replay_frame = 0
graph_video = []
graphCacheReplay = []
impactsAt = []
start_time = timeit.default_timer()

gyroCache = []
accelCache = []
gyroCacheReplay = []
accelCacheReplay = []

num = 1

class CSPS(tk.Frame):
    def __init__(self, parent):
        global flag
        tk.Frame.__init__(self, parent)

        self.tssensor = None

        self.x = [0 for i in range(100)]
        self.y = [0 for i in range(100)]
        self.z = [0 for i in range(100)]

        self.x2 = [0 for i in range(100)]
        self.y2 = [0 for i in range(100)]
        self.z2 = [0 for i in range(100)]

        self.canvas = tk.Canvas(self, background="black")
        self.canvas.bind("<Configure>", self.on_resize)
        self.canvas.grid(sticky="news")
        self.cwidth = self.canvas.winfo_width()
        self.cheight = self.canvas.winfo_height()
        self.heading = self.canvas.create_text(200, 25, tag='heading', fill="white", text="Sensor Data")
        self.canvas.create_line((0, 0, 0, 0), tag='X', fill='red', width=1)
        self.canvas.create_line((0, 0, 0, 0), tag='Y', fill='blue', width=1)
        self.canvas.create_line((0, 0, 0, 0), tag='Z', fill='green', width=1)
        self.canvas.create_line((0, 0, 0, 0), tag='T1_1', fill='white', dash=(1,), width=1)
        self.canvas.create_line((0, 0, 0, 0), tag='T1_2', fill='white', dash=(1,), width=1)
        self.text1 = self.canvas.create_text(200, 275, tag='text1', fill="white", text="Rotational Acceleration \n (Radians/sec)")
        self.canvas.create_line((0, 0, 0, 0), tag='X2', fill='red', width=1)
        self.canvas.create_line((0, 0, 0, 0), tag='Y2', fill='blue', width=1)
        self.canvas.create_line((0, 0, 0, 0), tag='Z2', fill='green', width=1)
        self.canvas.create_line((0, 0, 0, 0), tag='T2_1', fill='white', dash=(1,), width=1)
        self.canvas.create_line((0, 0, 0, 0), tag='T2_2', fill='white', dash=(1,), width=1)
        self.text2 = self.canvas.create_text(200, 575, tag='text2', fill="white", text="Linear Acceleration \n (Meters/sec^2)")

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.grid(sticky="news")
        parent.grid_rowconfigure(0, weight=1)
        parent.grid_columnconfigure(0, weight=1)

        if no_sensor_api == False:
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
                            flag = 1
                            break
                        self.tssensor = self.tsdongle.getSensorFromDongle(0)
                        if not self.tssensor:
                            flag = 1
                            break
                        break
                    if self.port_info.dev_type == 'WL-H3' or self.port_info.dev_type == 'WL':
                        self.com_port = self.cp
                        self.tssensor = ts_api.TSWLSensor(self.com_port)
                        break
                if not self.com_port:
                    flag = 1
            else:
                flag = 1

    def on_resize(self, event):
        self.replot()

    def bg1(self):
        global button
        button.configure(bg="green", state="normal")

    def video(self):
        global record_on, videoCache, replayCache, timeCache, replayTimeCache

        _, frame = capture.read()
        # frame = cv2.flip(frame, 1)

        # keeps caching for preimpact recording - Ishan
        if len(videoCache) > CONST_cacheLimit / 2:
            videoCache.pop(0)
            videoCache.append(frame)
            timeCache.pop(0)
            timeCache.append(timeit.default_timer() - start_time)
            # print("limit reached")
        else:
            videoCache.append(frame)
            timeCache.append(timeit.default_timer() - start_time)

        # in the event of an impact triggers postimpact recording into replayCache - Ishan
        if record_on:
            if len(replayCache) < CONST_cacheLimit:
                replayCache.append(frame)
                replayTimeCache.append(timeit.default_timer() - start_time)
            else:
                # turn off recording as entire replay has been recorded
                record_on = False

        curWidth = video.winfo_width()
        curHeight = video.winfo_height()
        maxsize = (curWidth, curHeight)
        frame = cv2.resize(frame, maxsize)
        cv2image = cv2.cvtColor(frame, cv2.COLOR_BGR2RGBA)
        img = Image.fromarray(cv2image)
        imgtk = ImageTk.PhotoImage(image=img)
        video.imgtk = imgtk
        video.configure(image=imgtk)

    def read_serial(self):
        global button, replayCache, record_on, gyroCacheReplay, accelCacheReplay, replay_frame, menu, fileMenu, replayTimeCache, impactsAt, graphCacheReplay, graph_video
        data = self.tssensor.getCorrectedGyroRate()
        data2 = self.tssensor.getCorrectedAccelerometerVector()
        x, y, z = data[0], data[1], data[2]
        x2, y2, z2 = data2[0], data2[1], data2[2]

        if len(gyroCache) > CONST_cacheLimit / 2:
            gyroCache.pop(0)
            accelCache.pop(0)
            gyroCache.append(data)
            accelCache.append(data2)
        else:
            gyroCache.append(data)
            accelCache.append(data2)

        if record_on:
            if len(gyroCacheReplay) < CONST_cacheLimit:
                gyroCacheReplay.append(data)
                accelCacheReplay.append(data2)
            else:
                record_on = False

        if x > 3 or x < -3 or y > 3 or y < -3 or z > 3 or z < -3 or x2 > 3 or x2 < -3 or y2 > 3 or y2 < -3 or z2 > 3 or z2 < -3:
            button.configure(bg="red", state="disabled")
            button.after(3000, self.bg1)
            # shifts pre impact to replay cache and enables after impact recording - Ishan
            record_on = True
            impactsAt.append(timeit.default_timer() - start_time)

            # discart impacts too far back - Ishan
            while float(impactsAt[len(impactsAt) - 1] - impactsAt[0]) > CONST_ImpactGroupTimeThreshold:
                impactsAt.pop(0)

            # copy precoreded data cache at impact                        
            replayCache = videoCache[:]
            replayTimeCache = timeCache[:]
            graphCacheReplay = graph_video[:]
            gyroCacheReplay = gyroCache[:]
            accelCacheReplay = accelCache[:]
            replay_frame = 0

        self.add(data, data2)
        self.after_idle(self.replot)

    def add(self, data, data2):
        self.x.append(float(data[0]))
        self.x = self.x[-100:]

        self.y.append(float(data[1]))
        self.y = self.y[-100:]

        self.z.append(float(data[2]))
        self.z = self.z[-100:]

        self.x2.append(float(data2[0]))
        self.x2 = self.x2[-100:]

        self.y2.append(float(data2[1]))
        self.y2 = self.y2[-100:]

        self.z2.append(float(data2[2]))
        self.z2 = self.z2[-100:]

    def replot(self):
        global record_on, graph_video, graphCacheReplay, num
        w = self.winfo_width()
        h = self.winfo_height()

        coordsX, coordsY, coordsZ = [], [], []
        coordsX2, coordsY2, coordsZ2 = [], [], []
        coordsT1_1, coordsT1_2, coordsT2_1, coordsT2_2 = [], [], [], []

        for n in range(0, 100):
            x = (w * n) / 100

            coordsX.append(x)
            coordsX.append(h - ((h * (self.x[n] + 150)) / 200.0))

            coordsY.append(x)
            coordsY.append(h - ((h * (self.y[n] + 150)) / 200.0))

            coordsZ.append(x)
            coordsZ.append(h - ((h * (self.z[n] + 150)) / 200.0))

            coordsT1_1.append(x)
            coordsT1_1.append(h - ((h * (147)) / 200.0))

            coordsT1_2.append(x)
            coordsT1_2.append(h - ((h * (153)) / 200.0))

            coordsX2.append(x)
            coordsX2.append(h - ((h * (self.x2[n] + 50)) / 200.0))
            coordsY2.append(x)
            coordsY2.append(h - ((h * (self.y2[n] + 50)) / 200.0))
            coordsZ2.append(x)
            coordsZ2.append(h - ((h * (self.z2[n] + 50)) / 200.0))

            coordsT2_1.append(x)
            coordsT2_1.append(h - ((h * (47)) / 200.0))

            coordsT2_2.append(x)
            coordsT2_2.append(h - ((h * (53)) / 200.0))

        self.canvas.coords('X', *coordsX)
        self.canvas.coords('Y', *coordsY)
        self.canvas.coords('Z', *coordsZ)

        self.canvas.coords('X2', *coordsX2)
        self.canvas.coords('Y2', *coordsY2)
        self.canvas.coords('Z2', *coordsZ2)

        self.canvas.coords('T1_1', *coordsT1_1)
        self.canvas.coords('T1_2', *coordsT1_2)

        self.canvas.coords('T2_1', *coordsT2_1)
        self.canvas.coords('T2_2', *coordsT2_2)

        w = CONST_videoWidth
        h = CONST_videoHeight

        coordsX, coordsY, coordsZ = [], [], []
        coordsX2, coordsY2, coordsZ2 = [], [], []
        coordsT1_1, coordsT1_2, coordsT2_1, coordsT2_2 = [], [], [], []

        for n in range(0, 100):
            x = (w * n) / 100

            coordsX.append(x)
            coordsX.append(h - ((h * (self.x[n] + 150)) / 200.0))

            coordsY.append(x)
            coordsY.append(h - ((h * (self.y[n] + 150)) / 200.0))

            coordsZ.append(x)
            coordsZ.append(h - ((h * (self.z[n] + 150)) / 200.0))

            coordsT1_1.append(x)
            coordsT1_1.append(h - ((h * (147)) / 200.0))

            coordsT1_2.append(x)
            coordsT1_2.append(h - ((h * (153)) / 200.0))

            coordsX2.append(x)
            coordsX2.append(h - ((h * (self.x2[n] + 50)) / 200.0))

            coordsY2.append(x)
            coordsY2.append(h - ((h * (self.y2[n] + 50)) / 200.0))

            coordsZ2.append(x)
            coordsZ2.append(h - ((h * (self.z2[n] + 50)) / 200.0))

            coordsT2_1.append(x)
            coordsT2_1.append(h - ((h * (47)) / 200.0))

            coordsT2_2.append(x)
            coordsT2_2.append(h - ((h * (53)) / 200.0))

        graphImage = Image.new("RGB", (CONST_videoWidth, CONST_videoHeight), "black")
        draw = ImageDraw.Draw(graphImage)
        font = ImageFont.truetype("arial.ttf", 16)
        draw.text((141, 50), "Sensor Data", fill="white", font=font)
        draw.line(coordsX, fill="blue", width=1)
        draw.line(coordsY, fill="red", width=1)
        draw.line(coordsZ, fill="green", width=1)
        draw.line(coordsT1_1, fill="white", width=1)
        draw.line(coordsT1_2, fill="white", width=1)
        draw.text((141, 350), "Rotational Acceleration \n (Radians/sec)", fill="white", font=font)
        draw.line(coordsX2, fill="blue", width=1)
        draw.line(coordsY2, fill="red", width=1)
        draw.line(coordsZ2, fill="green", width=1)
        draw.line(coordsT2_1, fill="white", width=1)
        draw.line(coordsT2_2, fill="white", width=1)
        draw.text((141, 800), "Linear Acceleration \n (Meters/sec^2)", fill="white", font=font)

        if len(graph_video) > CONST_cacheLimit / 2:
            graph_video.pop(0)
            graph_video.append(graphImage)
        else:
            graph_video.append(graphImage)

        if record_on:
            if len(graphCacheReplay) < CONST_cacheLimit:
                graphCacheReplay.append(graphImage)
            else:
                record_on = False

    def replay(self):
        global replay_video, replay_on, replay_frame, replay_sensor_graph, tx1, ty1, tz1, tx2, ty2, tz2, timer, simul

        rframe = replayCache[int(replay_frame)]

        curWidth = replay_video.winfo_width()
        curHeight = replay_video.winfo_height()
        maxsize = (curWidth, curHeight)
        rframe = cv2.resize(rframe, maxsize)
        cv2image = cv2.cvtColor(rframe, cv2.COLOR_BGR2RGBA)

        img = Image.fromarray(cv2image)
        imgtk = ImageTk.PhotoImage(image=img)
        replay_video.imgtk = imgtk
        replay_video.configure(image=imgtk)

        w = replay_sensor_graph.winfo_width()
        h = replay_sensor_graph.winfo_height()

        coordsX, coordsY, coordsZ = [], [], []
        coordsX2, coordsY2, coordsZ2 = [], [], []
        coordsT1_1, coordsT1_2, coordsT2_1, coordsT2_2 = [], [], [], []

        tx1.pop(0)
        ty1.pop(0)
        tz1.pop(0)

        tx2.pop(0)
        ty2.pop(0)
        tz2.pop(0)

        tx1.append(gyroCacheReplay[int(replay_frame)][0])
        ty1.append(gyroCacheReplay[int(replay_frame)][1])
        tz1.append(gyroCacheReplay[int(replay_frame)][2])

        tx2.append(accelCacheReplay[int(replay_frame)][0])
        ty2.append(accelCacheReplay[int(replay_frame)][1])
        tz2.append(accelCacheReplay[int(replay_frame)][2])

        for n in range(0, 100):
            x = (w * n) / 100

            coordsX.append(x)
            coordsX.append(h - ((h * (tx1[n] + 150)) / 200.0))

            coordsY.append(x)
            coordsY.append(h - ((h * (ty1[n] + 150)) / 200.0))

            coordsZ.append(x)
            coordsZ.append(h - ((h * (tz1[n] + 150)) / 200.0))

            coordsT1_1.append(x)
            coordsT1_1.append(h - ((h * (147)) / 200.0))

            coordsT1_2.append(x)
            coordsT1_2.append(h - ((h * (153)) / 200.0))

            coordsX2.append(x)
            coordsX2.append(h - ((h * (tx2[n] + 50)) / 200.0))

            coordsY2.append(x)
            coordsY2.append(h - ((h * (ty2[n] + 50)) / 200.0))

            coordsZ2.append(x)
            coordsZ2.append(h - ((h * (tz2[n] + 50)) / 200.0))

            coordsT2_1.append(x)
            coordsT2_1.append(h - ((h * (47)) / 200.0))

            coordsT2_2.append(x)
            coordsT2_2.append(h - ((h * (53)) / 200.0))

        replay_sensor_graph.coords('X', *coordsX)
        replay_sensor_graph.coords('Y', *coordsY)
        replay_sensor_graph.coords('Z', *coordsZ)

        replay_sensor_graph.coords('X2', *coordsX2)
        replay_sensor_graph.coords('Y2', *coordsY2)
        replay_sensor_graph.coords('Z2', *coordsZ2)

        replay_sensor_graph.coords('T1_1', *coordsT1_1)
        replay_sensor_graph.coords('T1_2', *coordsT1_2)

        replay_sensor_graph.coords('T2_1', *coordsT2_1)
        replay_sensor_graph.coords('T2_2', *coordsT2_2)

        msTimer = str(int(float(replayTimeCache[int(replay_frame)] - replayTimeCache[0]) * 1000)) + "ms"
        impactStartTime = str(int(float(impactsAt[0] - replayTimeCache[0]) * 1000)) + "ms"
        impactEndTime = str(int(float(impactsAt[len(impactsAt) - 1] - replayTimeCache[0]) * 1000)) + "ms"
        timer.configure(text="Impact from " + impactStartTime + " to " + impactEndTime + " (Time elapsed: " + msTimer + ")")

        replay_frame += 1 / CONST_slowDown

        bframeStart = 15
        bframeEnd = 104

        number = bframeStart + replay_frame / len(replayCache) * (bframeEnd - bframeStart)
        filename = "./Simulation/frame" + str(int(number)) + ".jpg"
        frame = cv2.imread(filename)
        cv2image = cv2.cvtColor(frame, cv2.COLOR_BGR2RGBA)
        img = Image.fromarray(cv2image)
        img = img.resize((247, 650), Image.ANTIALIAS)
        imgtk = ImageTk.PhotoImage(image=img)
        simul.imgtk = imgtk
        simul.configure(image=imgtk)

        # once youve replayed the recording stop and reset - Ishan
        if replay_frame == len(replayCache):
            replay_on = False
            replay_frame = 0
            tx1 = [0 for i in range(100)]
            ty1 = [0 for i in range(100)]
            tz1 = [0 for i in range(100)]

            tx2 = [0 for i in range(100)]
            ty2 = [0 for i in range(100)]
            tz2 = [0 for i in range(100)]


    def show(self):
        global flag, flag2

        self.video()

        if no_sensor_api == False:
            if flag:
                if flag2:
                    flag2 = 0
                    sensor.destroy()
                    tkMessageBox.showerror("Error", "No sensor data")
            else:
                self.read_serial()
        else:
            sensor.destroy()
            tkMessageBox.showerror("Error", "No sensor API")

        if flag3 and replay_on:
            self.replay()

        self.after(1, self.show)


def hello():
    hello = tk.Tk()
    hello.title("Welcome")

    text = tk.Label(hello)
    text.config(text="Welcome to CSPS")
    text.pack()
    hello.after(5000, lambda: hello.destroy())


def about():
    about = tk.Tk()
    about.title("About")

    text = tk.Label(about)
    text.config(text="You are using CSPS v1.0")
    text.pack()
    about.after(5000, lambda: about.destroy())


def exit():
    root.destroy()


def exit_sub():
    subRoot.destroy()


def setFlag():
    global flag3, temp, embox
    flag3 = 0
    if tkMessageBox.askyesno("Video", "Do you want to generate a video? If yes, kindly wait for some time..."):
        genVideo()
    subRoot.destroy()

def gen():
    global temp, embox, videoFileName
    
    xz = embox.get()

    if xz:
        sendEmail(xz)
    else:
        sendEmail()
    #temp.destroy()

def pickRandomImage():
    j = random.randint(1, 11)
    img = cv2.imread("./BS/Capture_" + str(j) + ".png")
    cv2.imshow('Img', img)


def call():
    global flag3, replay_video, subRoot, replay_sensor_graph, replay_sensor, replay_on, simul, timer, button1

    if flag3 == 1:
        setFlag()

    # flag indicating replay is on - Ishan
    replay_on = True
    flag3 = 1
    subMenu = tk.Menu()
    subRoot = tk.Toplevel()
    subRoot.title("Details")
    subRoot.config(menu=subMenu, background='black')
    subRoot.columnconfigure(0, weight=1)
    subRoot.rowconfigure(0, weight=1)
    subRoot.protocol('WM_DELETE_WINDOW', setFlag)
    subRoot.wm_state("zoomed")
    # subRoot.attributes("-fullscreen", True)
    file = tk.Menu(subMenu)
    subMenu.add_cascade(label="File", menu=file)
    file.add_command(label="Exit", command=exit_sub)

    replay_video = tk.Label(subRoot)
    replay_video.grid(row=1, column=0, sticky="news")
    replay_video.configure(width=CONST_panelWidth, height=CONST_panelHeight)

    replay_sensor = tk.Label(subRoot)
    replay_sensor.grid(row=1, column=1, sticky="news")
    replay_sensor.configure(width=CONST_panelWidth, height=CONST_panelHeight)
    replay_sensor.grid_rowconfigure(0, weight=1)
    replay_sensor.grid_columnconfigure(0, weight=1)

    replay_sensor_graph = tk.Canvas(replay_sensor, background="black")
    replay_sensor_graph.grid(sticky="news")
    replay_sensor_graph.grid_rowconfigure(0, weight=1)
    replay_sensor_graph.grid_columnconfigure(0, weight=1)
    rwidth = replay_sensor_graph.winfo_width()
    rheight = replay_sensor_graph.winfo_height()
    heading = replay_sensor_graph.create_text(200, 25, tag='heading', fill="white", text="Sensor Data")
    replay_sensor_graph.create_line((0, 0, 0, 0), tag='X', fill='red', width=1)
    replay_sensor_graph.create_line((0, 0, 0, 0), tag='Y', fill='blue', width=1)
    replay_sensor_graph.create_line((0, 0, 0, 0), tag='Z', fill='green', width=1)
    replay_sensor_graph.create_line((0, 0, 0, 0), tag='T1_1', fill='white', dash=(1,), width=1)
    replay_sensor_graph.create_line((0, 0, 0, 0), tag='T1_2', fill='white', dash=(1,), width=1)
    text1 = replay_sensor_graph.create_text(200, 275, tag='text1', fill="white", text="Rotational Acceleration \n (Radians/sec)")
    replay_sensor_graph.create_line((0, 0, 0, 0), tag='X2', fill='red', width=1)
    replay_sensor_graph.create_line((0, 0, 0, 0), tag='Y2', fill='blue', width=1)
    replay_sensor_graph.create_line((0, 0, 0, 0), tag='Z2', fill='green', width=1)
    replay_sensor_graph.create_line((0, 0, 0, 0), tag='T2_1', fill='white', dash=(1,), width=1)
    replay_sensor_graph.create_line((0, 0, 0, 0), tag='T2_2', fill='white', dash=(1,), width=1)
    text2 = replay_sensor_graph.create_text(200, 575, tag='text2', fill="white", text="Linear Acceleration \n (Meters/sec^2)")

    simul = tk.Label(subRoot)
    simul.grid(row=1, column=2, sticky="news")
    simul.configure(width=247, height=650)
    
    timer = tk.Label(subRoot, height=4, fg="red", bg="black", justify="center", font=("Courier", 16))
    timer.grid(row=0, column=0, columnspan=3, sticky="news")

    button1 = tk.Button(subRoot, text="Detailed Brain Reuslts", height=2, bg='black', fg='white', command=pickRandomImage)
    button1.grid(row=2, column=0, columnspan=3, sticky="news")

#send email
def sendEmail(email):
    global videoFileName
    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.starttls()
    server.login("coachsafeplaysafe@gmail.com", "QAZwsx123!@#")
    msg = """The video of the latest impact is at
    https://s3.amazonaws.com/cspsvideoupload/"""+ videoFileName 
    msg2 = "\r\n".join([
    "Subject: Link to video of impact",
    "",
    msg
    ])
    server.sendmail("coachsafeplaysafe@gmail.com", "coachsafeplaysafe@gmail.com", msg2)
    server.sendmail("coachsafeplaysafe@gmail.com", email, msg2)
    server.quit()


#upload to s3 server
def uploadFileS3(name):
    aws_access_key_id='AKIAI7XOEFVYPOUQKGYQ'
    aws_secret_access_key='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'
    conn = tinys3.Connection(aws_access_key_id, aws_secret_access_key)    
    f = open(name, 'rb')
    conn.upload(name, f, 'cspsvideoupload')


def genVideo():
    global temp, embox, videoFileName
    #maxsize = (CONST_videoWidth, CONST_videoHeight)
    videosize = (1206, 866)
    graphsize = (382, 866)
    simulationsize = (324, 866)
    fourcc = cv2.VideoWriter_fourcc(*'X264')
    out = cv2.VideoWriter(videoFileName, fourcc, 10.0, (1912, 866))

    curFrame = 0

    while (curFrame < len(replayCache)):

        videoFrame = replayCache[curFrame]
        videoFrame = cv2.resize(videoFrame, videosize)

        sensorFrame = np.array(graphCacheReplay[curFrame], dtype=np.uint8)
        sensorFrame = cv2.resize(sensorFrame, graphsize)

        offset = 25
        number = offset + curFrame
        if number > 104:
            number = 104
        filename = "./Simulation/frame" + str(int(number)) + ".jpg"
        simulFrame = cv2.imread(filename)
        simulFrame = cv2.resize(simulFrame, simulationsize)

        combinedFrame = np.concatenate((videoFrame, sensorFrame, simulFrame), axis=1)
        out.write(combinedFrame)

        if cv2.waitKey(1) & 0xFF == ord('q'):
            break

        curFrame += 1

    out.release()
    uploadFileS3(videoFileName)
    temp = tk.Tk()
    temp.title("Email Address")
    email = tk.Label(temp, text="Email")
    email.grid(row=0, column=0, sticky="news")
    embox = tk.Entry(temp)
    embox.grid(row=0, column=1, sticky="news")
    b = tk.Button(temp, text="Submit", command=gen)
    b.grid(row=1, column=0, columnspan=2, sticky="news")


root = tk.Tk()
menu = tk.Menu(root)
flag = 0
flag2 = 1
flag3 = 0
tx1 = [0 for i in range(100)]
ty1 = [0 for i in range(100)]
tz1 = [0 for i in range(100)]
tx2 = [0 for i in range(100)]
ty2 = [0 for i in range(100)]
tz2 = [0 for i in range(100)]
replay_video = None
replay_sensor = None
replay_sensor_graph = None
simul = None
timer = None
subRoot = None
button1 = None
embox = None
temp = None
videoFileName = "videoofimpact.mp4"
gray15 = (40, 40, 40)
root.title("CoachSafe PlaySafe")
root.iconbitmap(default='CSPS_HR.ico')
root.config(menu=menu, background='black')
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)
root.wm_state("zoomed")
# root.attributes("-fullscreen", True)
fileMenu = tk.Menu(menu)
menu.add_cascade(label="File", menu=fileMenu)
fileMenu.add_command(label="Hello", command=hello)
fileMenu.add_separator()
fileMenu.add_command(label="Exit", command=exit)

helpMenu = tk.Menu(menu)
menu.add_cascade(label="Help", menu=helpMenu)
helpMenu.add_command(label="About", command=about)

button = tk.Button(text="Show", height=8, bg='green', fg='black', state="disabled", command=call)
button.grid(row=1, column=0, columnspan=3, sticky="news")

video = tk.Label(root)
video.grid(row=0, column=0, sticky="news")
video.configure(width=CONST_panelWidth, height=CONST_panelHeight)
capture = cv2.VideoCapture(0)

sensor = tk.Label(root)
sensor.grid(row=0, column=1, sticky="news")
sensor.configure(width=CONST_panelWidth, height=CONST_panelHeight)
sensor.grid_rowconfigure(0, weight=1)
sensor.grid_columnconfigure(0, weight=1)
obj = CSPS(sensor)
obj.show()
root.mainloop()
