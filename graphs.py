import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
import sys
import scipy.stats as stats
import matplotlib.cm as cm

import matplotlib as mpl
mpl.rcParams["errorbar.capsize"] = 2
mpl.rcParams["lines.linewidth"] = 1
mpl.rcParams['pdf.fonttype'] = 42


np.set_printoptions( threshold=999999999999999)


def trevrolls(frates):
	r2=0.
	rs =0.
	n = float(frates.size)
	for i in range(frates.size): # np.nditer(frates):
		r2 += (frates[i]**2)/n
		rs += frates[i]/n
	return 1. - ((rs**2)/r2)
	

def loadspikesdat(filename, tduration):

	ff = open(filename, 'r') 
	fdata = ff.readlines()
	sx = len(fdata)
	sy = tduration;
	raster = np.zeros( (sx, sy) );
	nid=0
	for l in fdata:
		ar = np.fromstring(l, sep=' ' , dtype=int)
		raster[nid, ar] = 1
		raster[nid,0] =0 # XXX bug
		nid += 1

	return raster


def printpairstats(stat, name):
	print( name)
	print( stats.f_oneway( stat[0][0] , stat[1][0] ) ) 
	print( stats.f_oneway( stat[0][1] , stat[1][1] ) )
	print( stats.f_oneway( stat[0][0] , stat[0][1] ) )
	print( stats.f_oneway( stat[1][0] , stat[1][1] ) )

def exportcsv(data, name):
	mycsv = np.array( [data[0,0] , data[1,0], data[0, 1], data[1,1] ]);
	df = pd.DataFrame(mycsv.T)
	df.to_csv(name+".txt")

def label_diff(ax, i,j,text,X,Y):
    x = (X[i]+X[j])/2
    y = 1.1*max(Y[i], Y[j])
    dx = abs(X[i]-X[j])
    props = {'connectionstyle':'bar','arrowstyle':'-',\
                 'shrinkA':20,'shrinkB':20,'linewidth':1}
    ax.annotate(text, xy=(X[i],y-7), zorder=10, transform=ax.transData)
    #ax.text(.5, .5, "text")
    ax.annotate('', xy=(X[i],y), xytext=(X[j],y), arrowprops=props)

NPYRS=400
NINH=100

THRESHOLD=40.


def mkPlots():

	ncases = 2
	NRUNS=10

	saves = np.zeros( (30, NRUNS) );

	trs = np.zeros((ncases, NRUNS))
	corrs = np.zeros((ncases, NRUNS))

	sizeA = np.zeros((ncases, NRUNS))
	sizeB = np.zeros((ncases, NRUNS))
	fratesA = np.zeros((ncases, NRUNS))
	fratesB = np.zeros((ncases, NRUNS))

	trA = np.zeros((ncases, NRUNS))
	trB = np.zeros((ncases, NRUNS))


	pcorr = np.zeros((ncases, 4, NRUNS))

	NLIMIT = NPYRS +NINH

	for run in range(NRUNS):
		spikes = np.loadtxt( './data/control_%d/spikesperpattern.dat'%( run), dtype=float)
		spikes = spikes[:, 0:NLIMIT]



		activepop = np.sum(np.logical_or((spikes[1,:] >THRESHOLD),  ( spikes[0,:] >THRESHOLD) ))

		overlap = 100. * np.sum(np.logical_and((spikes[1,:] >THRESHOLD),  ( spikes[0,:] >THRESHOLD) )) / (activepop+0.0001)
		sizeA[0, run] = 100.* np.sum(spikes[0,:] >THRESHOLD) / NLIMIT   
		sizeB[0, run] = 100.*  np.sum(spikes[1,:] >THRESHOLD) / NLIMIT   
		trs[0, run] = overlap;


		cstart =100
		cend =200

		cf = np.corrcoef(spikes[0, cstart:cend], spikes[1,cstart:cend])
		corrs[0, run] = cf[0,1]



		fratesA[0,run] = np.mean(spikes[0,:])
		fratesB[0,run] = np.mean(spikes[1,:])

		trA[0,run] =  trevrolls(spikes[0,:])
		trB[0,run] =  trevrolls(spikes[1,:])


		spikes = np.loadtxt( './data/blocked_%d/spikesperpattern.dat'%( run), dtype=float)
		spikes = spikes[:, 0:NLIMIT]


		activepop = np.sum(np.logical_or((spikes[1,:] >THRESHOLD),  ( spikes[0,:] >THRESHOLD) ))

		overlap = 100. * np.sum(np.logical_and((spikes[1,:] >THRESHOLD),  ( spikes[0,:] >THRESHOLD) )) / (activepop+0.0001);
		sizeA[1, run] = 100.* np.sum(spikes[0,:] >THRESHOLD) / NLIMIT   
		sizeB[1, run] = 100.*  np.sum(spikes[1,:] >THRESHOLD) / NLIMIT   
		trs[1, run] = overlap;


		cf = np.corrcoef(spikes[0, cstart:cend], spikes[1, cstart:cend])
		corrs[1, run] = cf[0,1]

		fratesA[1,run] = np.mean(spikes[0,:])
		fratesB[1,run] = np.mean(spikes[1,:])

		trA[1,run] =  trevrolls(spikes[0,:])
		trB[1,run] =  trevrolls(spikes[1,:])



	plt.figure()
	plt.subplot(1,2,1)


	

	means = np.mean( trs, axis=1)
	stds = np.std( trs, axis=1)

	plt.ylim((0, 60));
	#plt.boxplot(  (trs[0], trs[1]), notch=True );
	plt.bar([1,2], means, yerr=stds, color=['indigo', 'purple'])
	plt.ylabel('Population Overlap CtxA & CtxB (%)');
	plt.xticks( [1, 2], ['Control', 'LC Block']);

	#np.savetxt("c_overlap.txt", trs, delimiter=',');
	saves[0:2] = np.array(trs)


	"""
	plt.subplot(1,2,2)
	means = np.mean(corrs, axis=1)
	stds = np.std(corrs, axis=1)
	plt.bar([1,2], means, yerr=stds, color=['indigo', 'purple'])
	plt.ylabel('Firing Rate Correlation CtxA / CtxB (%)');
	plt.xticks( [1, 2], ['Control', 'LC Block']);

	saves[3:5] = np.array(corrs)
	"""


	plt.subplot(1,2,2)

	ops = (sizeA[0], sizeB[0], sizeA[1], sizeB[1])
	means = np.mean( ops, axis=1)
	stds = np.std( ops, axis=1)
	plt.bar([1,2,3,4], means, yerr=stds,  color=['indigo', 'indigo', 'purple', 'purple'])



	#plt.ylim((20, 50));
	plt.xticks( [1, 2,3,4], ['CtxA\nControl', 'CtxB\nControl', 'CtxA\nLC Block', 'CtxB\nLC Block']);

	plt.ylabel('Active Population % (ff >10Hz)');

	st = (stats.f_oneway(sizeA[0], sizeB[0])) # , sizeA[1], sizeB[1]) )
	print(st)
	#plt.annotate('1One-way ANOVA %f'%(st[1]), xy = (0.3,0.9), xycoords='figure fraction' )

	saves[6:10] = np.array(ops)


	plt.figure()

	ops = (fratesA[0], fratesB[0], fratesA[1], fratesB[1])
	means = np.mean( ops, axis=1)
	stds = np.std( ops, axis=1)
	plt.bar([1,2,3,4], means, yerr=stds,  color=['indigo', 'indigo', 'purple', 'purple'])

	st = (stats.f_oneway(fratesA[0], fratesB[0], fratesA[1], fratesB[1]))
	print(st)
	plt.annotate('One-way ANOVA p= %g'%(st[1]), xy = (0.3,0.9), xycoords='figure fraction' )


	saves[11:15] = np.array(ops)


	#plt.ylim((15, 30));
	plt.ylabel('Mean Firing Rate (Hz)');
	plt.xticks( [1, 2,3,4], ['MemA\nControl', 'MemB\nControl', 'MemA\nLC Block', 'MemB\nLC Block']);


	plt.figure()
	plt.ylabel('Sparsity');
	ops = (trA[0], trB[0], trA[1], trB[1])
	means = np.mean( ops, axis=1)
	stds = np.std( ops, axis=1)
	plt.bar([1,2,3,4], means, yerr=stds,  color=['indigo', 'indigo', 'purple', 'purple'])
	plt.xticks( [1, 2,3,4], ['MemA\nControl', 'MemB\nControl', 'MemA\nLC Block', 'MemB\nLC Block']);

	saves[18:22] = np.array(ops)

	np.savetxt("saves.txt", saves.T, delimiter=',', fmt='%f');





def mkRampPlot():

	plt.figure()
	plt.ylabel('Mean Firing Rate (Hz)');

	selected = [0,5,6,]

	data = np.loadtxt( './data/ramp_data.txt');
	spikes = data[:,0:50];
	stds = np.std(spikes, axis=1)/np.sqrt(10);
	means = np.mean(spikes, axis=1)


	xlab = [0.01*v for v in range(means.shape[0])];
	plt.errorbar( xlab, means, stds, color='black' );


	spikes = data[:, 200:250]
	stds = np.std(spikes, axis=1)#/np.sqrt(10);
	means = np.mean(spikes, axis=1)
	plt.errorbar( xlab, means, stds,  color='royalblue');
	plt.legend(['Control', 'LC Inhibited']);


	plt.xlabel('Input current (nA)');
	plt.ylim((0, 25));
	#plt.xticks( [1, 2], ['Control', 'LC Block']);

def mkSamples():

	plt.figure()


	myvmax = 17;

	plt.subplot(2,3,1)
	data = np.loadtxt( './data/control_0/spikesperpattern.dat');
	data = data /4;
	spikes = data[0,0:400];
	spikes = spikes.reshape( (20,20))

	plt.imshow(spikes, vmin=0, vmax=myvmax)
	plt.axis('off')


	plt.subplot(2,3,2)
	spikes2 = data[1,0:400];
	spikes2 = spikes2.reshape( (20,20))
	plt.imshow(spikes, vmin=0, vmax=myvmax)
	plt.axis('off')

	plt.subplot(2,3,3)
	act =(np.logical_and( (spikes>10),  (spikes2>10 ) ))

	plt.imshow(act);
	plt.axis('off')


	plt.subplot(2,3,4)
	data = np.loadtxt( './data/blocked_3/spikesperpattern.dat');
	data = data /4;
	spikes = data[0,0:400]
	spikes = spikes.reshape( (20,20))
	plt.imshow(spikes, vmin=0, vmax=myvmax)
	plt.axis('off')

	plt.subplot(2,3,5)
	spikes2 = data[1,0:400];
	spikes2 = spikes2.reshape( (20,20))
	plt.imshow(spikes2, vmin=0, vmax=myvmax)
	plt.axis('off')
	#plt.colorbar()

	plt.subplot(2,3,6)
	act =(np.logical_and( (spikes>10),  (spikes2>10 ) ))
	print(act.sum())
	plt.imshow(act);
	plt.axis('off')

	plt.figure()
	plt.imshow(spikes)
	plt.colorbar();







	#plt.xlabel('Input current (nA)');
	#plt.xticks( [1, 2], ['Control', 'LC Block']);


def mkPairs(cond, run, label):
	raster = np.zeros((500, 8000));

	lines = open( './data/%s_%d/spikes.dat'%(cond,  run),'r').readlines(); 
	ln=0;
	for line in lines:
		cols = [int(n) for n in line.split()]
		raster[ln, cols ] = 1;
		ln +=1;
		if ln>=500: break


	memA = raster[:, 0:4000];
	memB = raster[:, 4000:8000];

	actA = (1*(memA.sum(axis=1)>40))
	actB= (1*(memB.sum(axis=1)>40))
	ov = (((actA + actB)>1))

	print(ov.sum())
	ovA = memA[ov, :]
	ovB = memB[ov, :]

	corrs = []
	for cs in range(40):
		sa= (ovA[:, cs:cs+100].sum(axis=1))
		sb=  (ovB[:, cs:cs+100].sum(axis=1))
		co = np.corrcoef(sa, sb)
		corrs.append(co[0,1])

	plt.figure()
	plt.title(label)
	plt.bar(range(40), corrs)
	plt.xlabel('Time (1/10 sec)')
	plt.ylabel('Corr. Coeff ')
	plt.ylim( (-0.2,0.5) )

def mkVolt():
	volt = np.loadtxt( './data/ramp_voltage.txt');
	plt.figure();
	plt.subplot(2,1,1);
	plt.ylabel("mV")
	plt.plot(volt[:,0]);
	plt.subplot(2,1,2);
	plt.plot(volt[:,1]);
	plt.ylabel("mV")
	plt.xlabel("msec")




mkRampPlot()

mkSamples()

mkPlots()


plt.show()
