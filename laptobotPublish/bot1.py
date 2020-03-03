import rospy
import csv
from fb5_torque_ctrl.msg import laptobot,PwmInput,encoderData
import numpy as np
import time
import math
from scipy.stats import multivariate_normal
from numpy.linalg import inv
from geometry_msgs.msg import TransformStamped


N_bots = 4

allData = [False,False,False,False]

wR=0
wL=0
wdotR=0
wdotL=0
encRPrev=0
encLPrev=0
wRprev=0
wLprev=0
wdotRprev=0
wdotLprev=0
codeStartTime=0
BotNumber = 0
PWM = PwmInput()
q=np.array([[0.0],[0.0],[0.0]])
q_prev=np.array([[0.0],[0.0],[0.0]])
previousCallbackTime=0
F_bar = np.matrix([[0.002],[0.002]])
tau = np.array([[0.0],[0.0]])

m=1.72 #mass in kg
d=0.065 #distance from CG to axis
R=0.085	#Semi-Distance between wheels
r=0.025 #Radius of wheels

K_wdot=0.088508			#(b+RJ/Kt). This is the coefficient of w_dot
K_w=0.0002705
K_tau_R=18000
K_tau_L=18000

k1=0.25
k2=0.05

wantData = True

origin = np.array([0.0,0.0])

vel = np.array([[0.0],[0.0]])

grid_Size = 5.0 #in meters
grid_Res = 0.02 #in meters. Each square is 2 cm width.

N_Grid_Points = int(grid_Size/grid_Res) #We gonna assume square grids. Else even the Voronoi partition function has to change.

x_grid=np.arange(-grid_Size/2+origin[0], grid_Size/2+origin[0], grid_Res)+grid_Res/2 #+grid_Res/2 gives the centroid
y_grid=np.arange(-grid_Size/2+origin[1], grid_Size/2+origin[1], grid_Res)+grid_Res/2

X_grid, Y_grid = np.meshgrid(x_grid,y_grid)

#pos_grid is a three 250x250 matrix with each point holding the coordinate of the centroid of a given grid square.
#eg. pos_grid[0,0,:] will give the 1st squares centroid as [-2.49,-2.49]
pos_grid = np.empty(X_grid.shape + (2,))
pos_grid[:, :, 0] = X_grid; pos_grid[:, :, 1] = Y_grid

bot_loc = np.empty((2,N_bots))
#print (np.shape(bot_loc),bot_loc[0][:])
#print(bot_loc)

mu0 = np.array([5*0.25, 5*0.25]) 	#The basis is formed by taking two gaussians on each diagnol and a constant function. 
mu1 = np.array([-5*0.25, -5*0.25])
mu2 = np.array([5*0.25, -5*0.25])
mu3 = np.array([-5*0.25, 5*0.25])
Sigma = np.array([[150.0/N_Grid_Points, 0],[0, 150.0/N_Grid_Points]])

rv0 = multivariate_normal(mu0, Sigma)
rv1 = multivariate_normal(mu1, Sigma)
rv2 = multivariate_normal(mu2, Sigma)
rv3 = multivariate_normal(mu3, Sigma)

#Adding all the basis ventors in 1 variable.
K=np.dstack((rv0.pdf(pos_grid), rv1.pdf(pos_grid),rv2.pdf(pos_grid),rv3.pdf(pos_grid),np.ones(pos_grid[:,:,0].shape)))

a_hat=np.array([3.0/1,1.0/100,1.0/1,1.0/100,1.0/100])

def quat2eul(qu):
	global phi
	global theta
	global psi
	#Computing Phi
	phi = math.atan2(2*(qu[0]*qu[1]+qu[2]*qu[3]),1-2*(qu[1]**2+qu[2]**2))

	#Computing Theta
	#Introducing a check to avoid numerical errors
	sinTheta=2*(qu[0]*qu[2]-qu[1]*qu[3])
	if sinTheta>=1:
		theta=math.pi/2
	else:
		theta=math.asin(sinTheta)

	#Computing Psi
	psi=math.atan2(2*(qu[0]*qu[3]+qu[2]*qu[1]),1-2*(qu[2]**2+qu[3]**2))
	return np.array([[phi],[theta],[psi]])


def cartesianDist(a,b):
	return math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

def voronoi(grid, Nbots, BotNo, locations):
	VPartition=[]
	N_Grid_Points = len(grid[:,0,0])
	for i in range(N_Grid_Points):
		for j in range(N_Grid_Points):	#This iterates over all points in the domain.
			inPartition=True 			#This stays one as long as the point is closer to the botIn question than any other point.
			for N in range(Nbots):
				if N!=BotNo and inPartition:
					inPartition = inPartition and cartesianDist(grid[i,j,:],locations[:,BotNo])<cartesianDist(grid[i,j,:],locations[:,N])
			if(inPartition):
				VPartition.append(np.array([i,j]))
	return VPartition

def Lv(partition,grid,K,a_hat,grid_res):
	phi_hat=a_hat[0]*K[:,:,0]+a_hat[1]*K[:,:,1]+a_hat[2]*K[:,:,2]+a_hat[3]*K[:,:,3]+a_hat[4]*K[:,:,4]
	LV=np.array([0.0,0.0])
	dq=grid_res*grid_res
	for point in partition:
		LV=LV+phi_hat[point[0],point[1]]*grid[point[0],point[1],:]*dq #integral(q*phi(q)*dq). dq is a constant equal to area of grid square.
	return LV

def Mv(partition,K,a_hat,grid_res):
	phi_hat=a_hat[0]*K[:,:,0]+a_hat[1]*K[:,:,1]+a_hat[2]*K[:,:,2]+a_hat[3]*K[:,:,3]+a_hat[4]*K[:,:,4]
	MV=0.0
	dq=grid_res*grid_res
	for point in partition:
		MV=MV+phi_hat[point[0],point[1]]*dq
	return MV

def S(q):
	global d
	s=np.matrix([[math.cos(q[2]), -d*math.sin(q[2])],[math.sin(q[2]), d*math.cos(q[2])],[0, 1]])
	return s

def S1(q):
	global d
	s1=np.matrix([[math.cos(q[2]), -d*math.sin(q[2])],[math.sin(q[2]), d*math.cos(q[2])]])
	return s1

def B(q):
	global R
	global r
	cont_transform=1/r*np.matrix([[math.cos(q[2]), math.cos(q[2])],[math.sin(q[2]), math.sin(q[2])],[R, -R]])
	return cont_transform

#By calculation we obtain B_bar as [[1/r,1/r],[R/r -R/r]]
def B_bar(q):
	b_bar=S(q).transpose()*B(q)
	return b_bar


def torque(q,v,Mv,Cv):
	global F_bar
	global k2
	global k1
	return inv(B_bar(q))*(-k1*Mv*S1(q).transpose()*(np.matrix(q[0:2])-np.matrix(Cv).transpose())-k2*np.matrix(v)+F_bar)


def pub():
   # global Lv
    global q,vel,wR,wL,wdotR,wdotL,pos_grid,N_bots,bot_loc,PWM
    pub = rospy.Publisher('pwmCmd1',PwmInput,queue_size = 10)   
   # print(bot_loc)
    partition = voronoi(pos_grid, N_bots, 1, bot_loc)
    #print(partition)
    mv = Mv(partition,K,a_hat,grid_Res) 
    lv = Lv(partition,pos_grid,K,a_hat,grid_Res)
   # print(mv,lv)
    Cv = lv/mv
   # print('vel',vel,'wr',wR,'wL',wL)
    tau=torque(q,vel,mv,Cv)
   # print('tau',tau)
    PWM.rightInput=K_tau_R*tau[0][0]+K_wdot*wdotR+K_w*wR
    PWM.leftInput=K_tau_L*tau[1][0]+K_wdot*wdotL+K_w*wL

    print('pwm_original',PWM)
    if(PWM.rightInput> 255):
	PWM.rightInput = 255
    if (PWM.rightInput< -255):
	PWM.rightInput = -255
    if(PWM.leftInput> 255):
	PWM.leftInput = 255
    if (PWM.leftInput< -255):
	PWM.leftInput = -255

    if(PWM.rightInput< 90 and PWM.rightInput>0):
	PWM.rightInput = 90
    if (PWM.leftInput< 90 and PWM.leftInput>0):
	PWM.leftInput = 90

    if(PWM.rightInput> -90 and PWM.rightInput<0):
	PWM.rightInput = -90
    if (PWM.leftInput> -90 and PWM.leftInput<0):
	PWM.leftInput = -90

    pub.publish(PWM)
    print('pwm',PWM)
    

def callback(data):
	global encRPrev
	global pwmInput
	global encLPrev
	global wRprev
	global wLprev
	global wdotRprev
	global wdotLprev
	global pwmInput
	global codeStartTime
	#rospy.loginfo(rospy.get_caller_id''() + "I heard %s", data)
	#['interval','encR','encL','wR','wL','wdotR','wdotL','wddotR','wddotL']
	#30 counts is one rotation of the wheel i.e. 2*pi radian
	duration=0.04 #25 Hz
	wR=(data.encoderR-encRPrev)*2*math.pi/30/duration
        wL=(data.encoderL-encLPrev)*2*math.pi/30/duration
        wdotR=(wR-wRprev)/duration
        wdotL=(wL-wLprev)/duration
        wddotR=(wdotR-wdotRprev)/duration
        wddotL=(wdotL-wdotLprev)/duration
	encRPrev=data.encoderR
	encLPrev=data.encoderL
	wRprev=wR
	wLprev=wL
	wdotRprev=wdotR
	wdotLprev=wdotL
        logTime=rospy.get_time()-codeStartTime
        row=[logTime,data.interval,data.encoderR,data.encoderL,wR,wL,wdotR,wdotL,wddotR,wddotL,pwmInput.rightInput,pwmInput.leftInput]
        with open('data.csv','ab') as myfile:
                writer=csv.writer(myfile)
                writer.writerow(row)
  
def callbackVICON(data,args):
	global q
        global allData 
	global BotNumber
	global bot_loc
	global wantData
	global q_prev
	global q_dot
	global vel
	global previousCallbackTime
	global wR
	global wL
	global wdotR
	global wdotL
	global w_prev
        global pos_grid
	
        if wantData:
                
                bot_loc[:,args]=np.array([data.transform.translation.x,data.transform.translation.y])
                allData[args] = True
                
                if args==BotNumber:
			if(allData[0] and allData[2] and allData[3]):
                     		wantData = False
                        q_prev[0][0]=q[0][0]
                        q_prev[1][0]=q[1][0]
                        q_prev[2][0]=q[2][0]
			delT=rospy.get_time()-previousCallbackTime
                        previousCallbackTime=rospy.get_time()
			#print(delT)
                        q[0][0]=data.transform.translation.x
                        q[1][0]=data.transform.translation.y
			eulerAng=quat2eul([data.transform.rotation.x, data.transform.rotation.y, data.transform.rotation.z, data.transform.rotation.w])
                        #The order of rotation so happens that phi is actually psi. Hence eulerAng[0] is used.
                        q[2][0]=eulerAng[0]
			q_dot = (q - q_prev)/delT
			vel = np.array(inv((S(q).transpose()*S(q)))*(S(q).transpose())*np.matrix(q_dot))
                        
                        if (wantData == False):
                            pub()
			    allData[0],allData[3],allData[2] = False,False,False
                        wantData = True  
        
                
    

def main():
    



    rospy.init_node('pub',anonymous = True)

    pub = rospy.Publisher('pwmCmd1',PwmInput,queue_size = 10)

    previousCallbackTime=rospy.get_time()

    rospy.Subscriber('encoderData1', encoderData, callback)

    rospy.Subscriber("/vicon/Karan0/Karan0", TransformStamped, callbackVICON,0)
    rospy.Subscriber("/vicon/Karan1/Karan1", TransformStamped, callbackVICON,1)
    rospy.Subscriber("/vicon/Karan2/Karan2", TransformStamped, callbackVICON,2)
    rospy.Subscriber("/vicon/Karan3/Karan3", TransformStamped, callbackVICON,3)
    rospy.spin()

   # Lv = np.array([0.0,0.0])
   # rate = rospy.Rate(10)
   # while not rospy.is_shutdown():
   #     t1 = time.time()
   #     pub()
   #     t2 = time.time()
   #     print(t2-t1)  
   #     rate.sleep()


if __name__ == '__main__':
    try:
        main()
    except rospy.ROSInterruptException:
        pass
    

