import numpy as np
from scipy.signal import convolve2d
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy.fft import fft, fftfreq
import math

# Width, height of the image. For good enough res and fast processing 360x360 is enough
nx, ny = 360,360
# Reaction rates. I studied (1.21,1,1) and (1,1,1)
alpha, beta, gamma = 1,1,1
N = 2000
def update(p,arr):
    # 1 - 2 = 0
    # 2 - 0 = 1
    # 0 - 1 = 2


    # Count the average amount of each species in the 9 cells around each cell
    # by convolution with the 3x3 array m.
    # basically averaging out over the immediate neighbours, convolution by a matrix of 1/9s is a faster method
    q = (p+1) % 2
    s = np.zeros((3, ny,nx))
    m = np.ones((3,3)) / 9
    for k in range(3):
        s[k] = convolve2d(arr[p,k], m, mode='same', boundary='wrap')
        
    # Apply the reaction difference equations
    arr[q,0] = s[0] + s[0]*(alpha*s[1] - gamma*s[2])
    arr[q,1] = s[1] + s[1]*(beta*s[2] - alpha*s[0])
    arr[q,2] = s[2] + s[2]*(gamma*s[0] - beta*s[1])
    #print(arr[q,0,0])
    # Ensure the species concentrations are kept within [0,1]. Normalisation.
    np.clip(arr[q], 0, 1, arr[q])
    return arr

# Initialize the array with uniform random [0,1] amounts of A, B and C.
arr = np.random.random(size=(2, 3, ny, nx))

#All this just image showing headache, cool part ends here
# Set up the image
#fig, ax = plt.subplots(2)
fig, ax = plt.subplots(2,1,gridspec_kw={'height_ratios' : [3,1]})
#fig, ax = plt.subplots()
prism = cm.get_cmap('prism', 256)
vapor = prism(np.linspace(0, 1, 256))

vapor[:,0] = 1 - vapor[:,0]
vapor[:,1] = 1 - vapor[:,1]
vapor[:,2] = 1 - vapor[:,2]
vaporwave = ListedColormap(vapor)
im = ax[0].imshow(arr[0,0], cmap=vaporwave)
ax[0].axis('off')
j=0
yy1 = []
yy2 = []
yy3 = []
yf1 = []
yf2 = []
yf3 = []
xf1 = []
xf2 = []
xf3 = []
def animate(i, arr):
    global yy1, yy2, yy3,yf1,yf2,yf3,xf1,xf2,xf3,j
    global N

   # """Update the image for iteration i of the Matplotlib animation."""
    q = ((i % 2) + 1) % 2
    print(j)

    # Draw x and y lists

    #if(j==250):
     #   plt.savefig("C:/music/bzrxn.svg")
    if(j!=N):
        yy1.append(arr[0,0,100,100])
        yy2.append(arr[0,1,100,100])
        yy3.append(arr[0,2,100,100])
    if(j==N):


        # sample spacing
        T = 1.0 / 2000.0

        yf1 = fft(yy1)
        xf1 = fftfreq(N, T)[:N//2]
        yf2 = fft(yy2)
        xf2 = fftfreq(N, T)[:N//2]
        yf3 = fft(yy3)
        xf3 = fftfreq(N, T)[:N//2]

    if(j>100):
        ydraw = yy1[-100:]

        ax[1].clear()
        ax[1].set_facecolor('#000080')
        ax[1].plot(range(0,100), ydraw, color = '#39ff14')




    arr = update(i % 2, arr)

    im.set_array(arr[i % 2, 0])
    j = j + 1
    return [im]

anim = animation.FuncAnimation(fig, animate, frames=200, interval=5,
                               blit=False, fargs=(arr,))


plt.tight_layout()

plt.show()
#logarithmic plot of fft frequencies vs time. First few frequencies decrease exponentially, so is sloping down line in logarithmic plot
#this is predicted by my math model
plt.plot(xf1, 2.0/N * np.abs(yf1[0:N//2]), color = 'red')
plt.grid()
plt.yscale("log")
plt.show()
plt.plot(xf2, 2.0/N * np.abs(yf2[0:N//2]), color = 'blue')
plt.grid()
plt.yscale("log")
plt.show()
plt.plot(xf3, 2.0/N * np.abs(yf3[0:N//2]), color = 'green')
plt.grid()
plt.yscale("log")
plt.show()
#first few phase differences will be integer multiple of 120 degerees. Also predicted but matching for only first two.
print("First vs Second phase differences for 4 peak harmonics:")
print((np.angle(yf1[120])-np.angle(yf2[120]))*(180/math.pi))
print((np.angle(yf1[240])-np.angle(yf2[240]))*(180/math.pi))
print((np.angle(yf1[360])-np.angle(yf2[360]))*(180/math.pi))
print((np.angle(yf1[480])-np.angle(yf2[480]))*(180/math.pi))

print("Second vs Third phase differences for 4 peak harmonics:")
print((np.angle(yf2[120])-np.angle(yf3[120]))*(180/math.pi))
print((np.angle(yf2[240])-np.angle(yf3[240]))*(180/math.pi))
print((np.angle(yf2[360])-np.angle(yf3[360]))*(180/math.pi))
print((np.angle(yf2[480])-np.angle(yf3[480]))*(180/math.pi))

print("Third vs First phase differences for 4 peak harmonics:")
print((np.angle(yf3[120])-np.angle(yf1[120]))*(180/math.pi))
print((np.angle(yf3[240])-np.angle(yf1[240]))*(180/math.pi))
print((np.angle(yf3[360])-np.angle(yf1[360]))*(180/math.pi))
print((np.angle(yf3[480])-np.angle(yf1[480]))*(180/math.pi))
