import numpy as np
import time
import sys
import cv2
import math
import random
import timeit
import tinys3
import os
import smtplib
import datetime

def makeS3Conn(aws_access_key_id,aws_secret_access_key,bucketName):
    conn = tinys3.Connection(aws_access_key_id, aws_secret_access_key,default_bucket = bucketName)
    return conn

def uploadFileS3(s3Conn, fName,cloudLimit):
    log = s3Conn.list()
    
    filenames = []
    filesizes = []
    cloudUsage = 0
    for i in log:
        filenames.append(str(i['key']))
        filesizes.append(int(i['size'])/1024)
        cloudUsage += filesizes[-1]

    cSize = float(os.stat(fName).st_size/1024)

    print("currrent log is")
    for i in range(len(filenames)):
        print(filenames[i],filesizes[i])
    print("total usage is",cloudUsage)

    print("going to add",fName,"with size",cSize)
    while (cSize + cloudUsage) > cloudLimit:
        print("usage = ",cloudUsage,"limit =",cloudLimit)
        print("so removing",filenames[0],"with size",filesizes[0]) 
        s3Conn.delete(filenames[0])
        cloudUsage -= filesizes[0]
        filenames.pop(0)
        filesizes.pop(0)

    f = open(fName, 'rb')
    s3Conn.upload(fName, f)
    
def genVideo():
        CONST_videoWidth = 300
        CONST_videoHeight = 300

        maxsize = (CONST_videoWidth, CONST_videoHeight)
        fourcc = cv2.VideoWriter_fourcc(*'XVID')
        cdatetime = str(datetime.datetime.now())
        cdate = cdatetime.split(" ")[0]

        ctime = cdatetime.split(" ")[1]
        ctime = ctime.split(".")[0]
        ctime = ctime.replace(":","-")
        
        videoFileName = "Date"+cdate+"Time"+ctime+"video.avi"
        print(videoFileName)
        out = cv2.VideoWriter(videoFileName, fourcc, 10.0, maxsize)
        curFrame = 0

        while (curFrame < 104):
                offset = 15
                number = offset + curFrame
                if number > 104:
                        number = 104
                filename = "./Simulation/frame" + str(int(number)) + ".jpg"
                simulFrame = cv2.imread(filename)
                simulFrame = cv2.resize(simulFrame, maxsize)

                out.write(simulFrame)
                curFrame += 1

        out.release()
        return videoFileName
    
