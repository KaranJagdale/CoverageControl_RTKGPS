#!/usr/bin/env python
import rospy
import serial
import time
import struct
from fb5_torque_ctrl.msg import PwmInput
from fb5_torque_ctrl.msg import encoderData
ser=serial.Serial('/dev/ttyUSB0',9600)

ser.flushInput()
ser.flushOutput()
ser.write('8')
def callback(pwmInput):
	ser.write('A') # ~ is \x7E
	rightPWM=pwmInput.rightInput
    	leftPWM=pwmInput.leftInput

    	#Ensuring the wheels are set to rotate in the correct direction for each scenario.
    	if (rightPWM>=0) and (leftPWM>=0):
        	ser.write('8')
    	elif (rightPWM<0) and (leftPWM>=0):
        	ser.write('6')
    	elif (rightPWM>=0) and (leftPWM<0):
        	ser.write('4')
    	elif (rightPWM<0) and (leftPWM<0):
        	ser.write('2')
    	else:
        	ser.write('5') #Should never occur
	lowPWM_R=90
    	#To ensure the signal sent is never larger than a byte. Probably happens automatically but safegaurd.
    	if (rightPWM>255) or (rightPWM<-255):
        	rightPWM=255
	elif (abs(rightPWM)<lowPWM_R):  #Below 90 PWM we observe that bot doesn't move at all hence we force it along
		rightPWM=lowPWM_R
	lowPWM_L=90
    	if (leftPWM>255) or (leftPWM<-255):
        	leftPWM=255
	elif (abs(leftPWM)<lowPWM_L):  #Below 90 PWM we observe that bot doesn't move at all hence we force it along
		leftPWM=lowPWM_L

    	#The C code on the bot accepts left motor velocity first.
    	ser.write(struct.pack('>B',abs(leftPWM)))
    	ser.write(struct.pack('>B',abs(rightPWM)))

def encoderOut():
	rospy.Subscriber('pwmCmd', PwmInput, callback)
	pub_encoder=rospy.Publisher('encoderData1',encoderData,queue_size=10)
	rospy.init_node('encoderOut',anonymous=True)
	rate = rospy.Rate(30) #Since bot sends data at 25hz the publisher will be forced to slow down
	while not rospy.is_shutdown():
		i=0
		while i<10:
			encoder=encoderData()
			encoder.interval=0
			encoder.encoderR=0
			encoder.encoderL=0
			encoder.chksum=0
			startTime=rospy.get_time()
			byte1=ser.read()
                	if byte1 is 'A':
                        	encoderRightH=ser.read()
                        	encoderRightL=ser.read()
                        	encoderLeftH=ser.read()
                        	encoderLeftL=ser.read()
                        	chksum=ser.read()
                        	encoder.encoderR=ord(encoderRightH)*256+ord(encoderRightL)
                        	encoder.encoderL=ord(encoderLeftH)*256+ord(encoderLeftL)
                        	sum=ord(encoderRightH)+ord(encoderRightL)+ord(encoderLeftH)+ord(encoderLeftL)
                        	#Byte2 and byte3 are the higher significant and lower significant bits respectively
                        	encoder.chksum=ord(chksum)==sum%256
        		encoder.interval = rospy.get_time()-startTime
        		#rospy.loginfo(encoder)
        		pub_encoder.publish(encoder)
			rate.sleep()

if __name__ == '__main__':
    try:
        encoderOut()
    except rospy.ROSInterruptException:
        pass
