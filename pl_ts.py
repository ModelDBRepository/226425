import numpy as np
import matplotlib.pyplot as plt
import struct as st

N=9;

fig1 = plt.figure(1, (9, 8))

area = 1.0; #1.0 - use full figure, 0.9 - 90%

#margins
borderX_l = 0.05; 
borderX_r = 0.03;
borderY_b = 0.06;
borderY_t = 0.03;

#grid numbers
numX = 2;
numY=2;

factorX = 0.9 # factor to have space between panels
factorY = 0.9 # factor to have space between panels
widthX = factorX*(area - (borderX_l+borderX_r))/numX
widthY = factorY*(area - (borderY_t+borderY_b))/numY

X_fig = np.linspace(borderX_l,area-borderX_r-widthX,numX);
Y_fig = np.linspace(borderY_b,area-borderY_t-widthY,numY);

fig1 = plt.figure(1, figsize = (9,8));

ax1 = fig1.add_axes([X_fig[0],Y_fig[1], widthX, widthY])
ax2 = fig1.add_axes([X_fig[1],Y_fig[1], widthX, widthY])
ax3 = fig1.add_axes([X_fig[0],Y_fig[0], widthX, widthY])
ax4 = fig1.add_axes([X_fig[1],Y_fig[0], widthX, widthY])

def sens(x):
	return 3./(1.+np.exp(-1.0*(x-6.)))

f = open('data1');
f1 = open('data2');
f2 = open('data3');

x = np.fromfile(f,dtype=np.float64);
x = np.resize(x,(np.size(x)/N,N));
s = np.fromfile(f2,dtype=np.float64);
s = np.resize(s,(np.size(s)/N,N));
I_ext = np.fromfile(f1,dtype=np.float64);



num = np.linspace(1,N,N);
#ax1 = fig1.add_subplot(323);
#ax1.plot(I_ext[N:2*N],num,'go-',linewidth=1,markersize=8);
#ax1.set_xlim(0.5,9.5)
#ax1.set_ylabel('neuron #')
#ax1.set_xticks([2.,8.])

X = np.linspace(0,12,100)
S = np.linspace(-2,2.,9)
for shift in S:
	ax1.plot(X,sens(X+shift),'b-', linewidth = 1);
ax1.set_xlim(0.0,11);
ax1.set_ylabel('Response',size=24);
ax1.set_xticks([]);
ax1.plot([2.,2.],[0,3],'r--', linewidth = 4)
ax1.plot([6.0,6.0],[0,3],'g--', linewidth = 4)
ax1.plot([9.,9.],[0,3],'b--', linewidth = 4)
ax1.set_xlabel('Odor concentration', size = 24)
ax1.xaxis.set_label_coords(0.5,-0.05)
ax1.set_yticks([])

for I in range(0,N):
	sp = np.where(x[:,I]>20.)[0];
	for i in range(0,np.size(sp)):
		#ax3.plot([sp[i]*0.001*0.01,sp[i]*0.001*0.01],[(1.0+I)/float(N),I/float(N)+0.001],'b-',linewidth=2,rasterized='True');
		ax3.plot([sp[i]*0.001*0.01,sp[i]*0.001*0.01],[(1.0+I)/float(N),I/float(N)+0.001],'b-',linewidth=2);

ax3.set_yticks(np.linspace(0,1,N));
ax3.set_yticklabels(np.arange(1,N+1))
ax3.set_ylim(0,1.0);
ax3.set_xlim(2.,4.)
ax3.set_xticks(np.linspace(2.0,4.0,2));
ax3.set_xticklabels(np.linspace(0.0,2.0,2));

ax3.set_ylabel('LN #', size = 24)

ax3.set_xlabel('time [s]', size = 24)
ax3.xaxis.set_label_coords(0.5,-0.05)

ax2.plot(num,-I_ext[0:N],'go-',linewidth=1,markersize=8);
#ax2.set_yticklabels([]);
ax2.set_ylabel('Input to LNs '+r'$ I^{ext}$', size = 24)
#ax2.set_yticks([1,9])
ax2.set_xlabel('LN #',size = 24)


f = open('data1_1');
f1 = open('data2_1');
f2 = open('data3_1');

x = np.fromfile(f,dtype=np.float64);
x = np.resize(x,(np.size(x)/N,N));
s = np.fromfile(f2,dtype=np.float64);
s = np.resize(s,(np.size(s)/N,N));
I_ext = np.fromfile(f1,dtype=np.float64);



num = np.linspace(1,N,N);
#fig1 = plt.figure(1);
#ax1 = fig1.add_subplot(323);
#ax1.plot(I_ext[N:2*N],num,'rs-',linewidth=1,markersize=8);

#ax3 = fig1.add_subplot(324);
#for I in range(0,N):
#	sp = np.where(x[:,I]>20.)[0];
#	for i in range(0,np.size(sp)):
#		ax3.plot([sp[i]*0.001*0.01,sp[i]*0.001*0.01],[(1.0+I)/float(N),I/float(N)+0.001],'b-',linewidth=2,rasterized='True');
#
#ax3.set_yticks(np.linspace(0,1,N));
#ax3.set_yticklabels(np.linspace(1,N,N).astype(int));
#ax3.set_ylim(0,1.0);
#ax3.set_xlim(2.0,4.0);
#ax3.set_xticks(np.linspace(2.0,4.0,2));
#ax3.set_xticklabels(np.linspace(0.0,2.0,2));
#ax3.set_ylabel('neuron #');


ax2.plot(num,-I_ext[0:N],'rs-',linewidth=1,markersize=8);

f = open('data1_2');
f1 = open('data2_2');
f2 = open('data3_2');

x = np.fromfile(f,dtype=np.float64);
x = np.resize(x,(np.size(x)/N,N));
s = np.fromfile(f2,dtype=np.float64);
s = np.resize(s,(np.size(s)/N,N));
I_ext = np.fromfile(f1,dtype=np.float64);



#num = np.linspace(1,N,N);
#fig1 = plt.figure(1);
#ax1 = fig1.add_subplot(323);
#ax1.plot(I_ext[N:2*N],num,'bv-',linewidth=1,markersize=8);


for I in range(0,N):
	sp = np.where(x[:,I]>20.)[0];
	for i in range(0,np.size(sp)):
		#ax4.plot([sp[i]*0.001*0.01,sp[i]*0.001*0.01],[(1.0+I)/float(N),I/float(N)+0.001],'b-',linewidth=2,rasterized='True');
		ax4.plot([sp[i]*0.001*0.01,sp[i]*0.001*0.01],[(1.0+I)/float(N),I/float(N)+0.001],'b-',linewidth=2);

ax4.set_yticks(np.linspace(0,1,N));
ax4.set_yticklabels(np.arange(1,N+1))
ax4.set_ylim(0,1.0);
ax4.set_xlim(2.,4.)
ax4.set_xticks(np.linspace(2.0,4.0,2));
ax4.set_xticklabels(np.linspace(0.0,2.0,2));

ax4.set_ylabel('LN #', size = 24)

ax4.set_xlabel('time [s]', size = 24)
ax4.xaxis.set_label_coords(0.5,-0.05)


ax2.plot(num,-I_ext[0:N],'bv-',linewidth=1,markersize=8);


#plt.tight_layout()

plt.show();

