import numpy as np
import matplotlib.pyplot as plt
R0=6.378e8
d=1e3
theta=70
theta_rad=theta/180*np.pi
slope=-np.tan(theta_rad)
det=[0,R0-d] #x,y coords
ys=np.linspace(-R0+d,R0-d,10**5)
def f(y):
    return np.abs(y-slope*(R0**2-y**2)**(1/2)-(R0-d))

fs=np.zeros(np.size(ys))
for i in range(np.size(fs)):
    fs[i]=f(ys[i])
min_y=ys[np.argmin(fs)]
print(np.min(fs))
print(np.argmin(fs))
print(ys[np.argmin(fs)])
x=(R0**2-min_y**2)**(1/2)
print(x,min_y)
print((x**2+min_y**2)**(1/2),180/np.pi*np.arctan(min_y/x))

'''
plt.figure(1)
plt.plot(ys,fs)
plt.show()

'''

# it Works to find the initial starting point!!!!
#how about stepping along the path

chord_length=((x-0)**2+(min_y-R0-d)**2)**(1/2)
print(chord_length)

step_length=chord_length/20
x_stepped=x-step_length*np.cos(theta_rad)
y_stepped=min_y+step_length*np.sin(theta_rad)
print(0,R0-d)
print(x_stepped,y_stepped)
print(x,min_y)
def step(x0,y0,step_length):
    x1=x0-step_length*np.cos(theta_rad)
    y1=y0+step_length*np.sin(theta_rad)
    return x1,y1
x_s=[x]
y_s=[min_y]
i=0
stepped_length=0
while(stepped_length<chord_length):
    x_temp,y_temp=step(x_s[i],y_s[i],step_length)
    x_s.append(x_temp)
    y_s.append(y_temp)
    stepped_length=stepped_length+step_length
    i=i+1



print()
print(R0-d,90)
print((x_stepped**2+y_stepped**2)**(1/2),180/np.pi*np.arctan(y_stepped/x_stepped))
print((x**2+min_y**2)**(1/2),180/np.pi*np.arctan(min_y/x))
y_circle=np.linspace(-R0,R0,10**5)
x_circle=np.zeros(np.size(y_circle))
for i in range(np.size(y_circle)):
    x_circle[i]=(R0**2-y_circle[i]**2)**(1/2)
plt.figure(1)
plt.axis([-R0-10**8,R0+10**8,-R0-10**8,R0+10**8])
plt.scatter([0,x_stepped,x],[R0-d,y_stepped,min_y])
plt.plot(x_circle,y_circle)
plt.scatter(x_s,y_s)
plt.show()

# worked again!!!

