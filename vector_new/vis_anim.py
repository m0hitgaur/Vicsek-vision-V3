import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import glob
from itertools import product
import os
import matplotlib.animation as animation
from  matplotlib.animation import FFMpegWriter
current_directory = os.getcwd()

import pandas as pd     



# Parameters

noise="0.05"#input("Noise : ")
angle="45"#input("Angle : ")
trial="0"#input("Trial : ")


detail=pd.read_csv("data/parameters.csv")
numberoftrial=int(detail["trial"][0])
maxiter=int(detail["maxiter"][0])
Number_of_agents=int(detail["N"][0])
Lx=float(detail["Lx"][0])
Ly=float(detail["Ly"][0])
v0=float(detail["v0"][0])
dt=float(detail["dt"][0])



times=[]
for i in range( maxiter):
    if(i<10):tf=1
    if(i>10):tf=10
    if(i>100):tf=50
    if(i>1000):tf=100
    if(i%tf==0):times.append(i)


# Set up the figure and axis
fig, ax = plt.subplots()
ax.set_aspect('equal')
quiver = ax.quiver(0, 0, 0, 0, angles='xy', scale_units='xy', scale=1)




colour=['r','b','g','m']




def update(t):
    global quiver
    quiver.remove()  # Remove the previous quiver 


    # Load position and angle data
    data_config=pd.read_csv(f"data/config_data/trial_{trial}/config_{times[t]}.csv")
    # Extract data for the specific time step
    px = data_config["x"]
    py = data_config["y"]
    
    # Compute velocity components
    vx = data_config["vx"]
    vy = data_config["vy"]
    
    # Plot the quiver plot
    #ax.set_xlim(Lx)
    #ax.set_ylim(Ly)
    ax.set_title("N="+str(Number_of_agents)+r"|$\eta=$"+str(noise)+r"|$\alpha=$"+str(angle)+r"|$t=$"+str(times[t]))
    quiver=ax.quiver(px, py, vx, vy,angles='xy', scale_units='xy', scale=0.01)
    
    return quiver

plt.show()



    

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=len(times), blit=False)
writer = FFMpegWriter(fps=10, bitrate=-1)   # -1 = auto bitrate
ani.save(f"anim.mp4", writer=writer)


