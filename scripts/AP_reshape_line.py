# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 14:22:31 2019

@author: karan
"""



import rospy
from fb5_torque_ctrl.msg import Array
import numpy as np
import time
import math

N_bots = 4
#BotN = 2

origin = np.array([0.0,0.0])

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
#print(X_grid)

bot_loc = np.empty((2,N_bots))


# Bot locations
bot_loc[:,0] = np.array([0.2 ,0.3])
bot_loc[:,1] = np.array([-1.2 ,2])
bot_loc[:,2] = np.array([-0.7 ,-0.3])
bot_loc[:,3] = np.array([0.8 ,2.1])


def checkside(x,y,n,m):  #(x,y) is the point to check and n,m is the pair of bots
    x1 = bot_loc[0,n]
    y1 = bot_loc[1,n]
    x2 = bot_loc[0,m]
    y2 = bot_loc[1,m]   
    if y1!=y2:
        return (x1-x2)/(y2-y1)*(x-(x1+x2)/2)-y+(y1+y2)/2      
    else:
        return x - (x1+x2)/2

X = []
Y = []
X1 = []
Y1 = []
X2 = []
Y2 = []
partition = []

def pointappend(grid,BotNo):
    del1 = []
    del2 = []
    global partition,X,Y,X1,Y1,X2,Y2,bot_loc
    if BotNo == 1:
        Rb1,Rb2,Rb3 = 0,2,3
    if BotNo == 2:
        Rb1,Rb2,Rb3 = 0,1,3
    if BotNo == 0:
        Rb1,Rb2,Rb3 = 1,2,3
    if BotNo == 3:
        Rb1,Rb2,Rb3 = 0,2,1

    botSign = checkside(bot_loc[0,BotNo],bot_loc[1,BotNo],BotNo,Rb1)
    #print(botSign)
    for i in range(N_Grid_Points):
        for j in range(N_Grid_Points):
            
            if checkside(grid[i,j,0],grid[i,j,1],BotNo,Rb1)*botSign > 0:
                partition.append(np.array([i,j]))
    
    #Plotting partition in stage1

#    for m in partition:
#        
#        X1.append(grid[m[0],m[1],0])
#        Y1.append(grid[m[0],m[1],1])
#    
#    X1 = np.array(X1)
#    Y1 = np.array(Y1)
#    plt.figure(1)    
#    plt.plot(X1,Y1,'b.')
   
#    plt.show() 
    
    
    
    botSign = checkside(bot_loc[0,BotNo],bot_loc[1,BotNo],BotNo,Rb2)
    for k in range(len(partition)):
        if checkside(grid[partition[k][0],partition[k][1],0],grid[partition[k][0],partition[k][1],1],BotNo,Rb2)*botSign < 0:
            del1.append(int(k))
    count = 0
    #print(del1)
    for i in del1:
        #print(i-count,'PL',len(partition))
        
        del partition[i-count]
        #partition.remove(i-count)
        count = count+1
    
    
    
    
    #Partition in stage 2

#    for m in partition1:
        
#        X2.append(grid[m[0],m[1],0])
#        Y2.append(grid[m[0],m[1],1])
   
#    X2 = np.array(X2)
#    Y2 = np.array(Y2)
#    plt.figure(4)    
#    plt.plot(X2,Y2,'g.')
#   
#    plt.show()
    
#    print(len(partition))
    botSign = checkside(bot_loc[0,BotNo],bot_loc[1,BotNo],BotNo,Rb3)
   
    for l in range(len(partition)):
        if checkside(grid[partition[l][0],partition[l][1],0],grid[partition[l][0],partition[l][1],1],BotNo,Rb3)*botSign < 0:
            del2.append(int(l))
    #print(len(partition))
    #print(del2)
    count = 0
    for i in del2:
        del partition[i-count]
        #partition.remove(i-count)
        count = count+1
    
        
   # for i in partition:
    #    i = np.array(i)
#    print(len(partition))
    
    
    return partition

def pub():
    rospy.init_node('pub',anonymous = True)
    t1 = time.time()
    partition0 = pointappend(pos_grid, 0)
    partition1 = pointappend(pos_grid, 1)
    partition2 = pointappend(pos_grid, 2)
    partition3 = pointappend(pos_grid, 3)

    pub0 = rospy.Publisher('AP0',Array,queue_size = 10)
    pub1 = rospy.Publisher('AP1',Array,queue_size = 10) 
    pub2 = rospy.Publisher('AP2',Array,queue_size = 10)
    pub3 = rospy.Publisher('AP3',Array,queue_size = 10)
    
    a0 = Array()
    a1 = Array()
    a2 = Array()
    a3 = Array()
    #data = np.array([[1.2,2.3],[3.4,4.5],[3.5,5.2]])

    data0 = np.array(partition0)
    data1 = np.array(partition1)
    data2 = np.array(partition2)
    data3 = np.array(partition3)

    a0.len = len(data0)
    a1.len = len(data1)
    a2.len = len(data2)
    a3.len = len(data3)

    a0.data = np.reshape(data0,2*(a0.len))
    a1.data = np.reshape(data1,2*(a1.len))
    a2.data = np.reshape(data2,2*(a2.len))
    a3.data = np.reshape(data3,2*(a3.len))

    rate = rospy.Rate(1)
    #while not rospy.is_shutdown():
    #print(data)
    t2 = time.time()
    print("time to publish - ",t2 - t1)
    while not rospy.is_shutdown():
        #pub0.publish(a0)
        #pub1.publish(a1)
        pub2.publish(a2) 
        #pub3.publish(a3)
        rate.sleep()
        
if __name__ == '__main__':
    try:
        pub()
    except rospy.ROSInterruptException:
        pass
