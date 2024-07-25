import uproot
import os
from numpy import loadtxt
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import awkward as ak
from array import array
import scipy.optimize as scp
import matplotlib as mpl
import matplotlib.pylab as plt

plot_gaincomp=False
save_multiplicity=False

#print all elements of an array
np.set_printoptions(threshold=sys.maxsize)

def flat_array(akward_array):
  if (akward_array.ndim > 1):
    flatted_arr = ak.flatten(akward_array)
    return flatted_arr 
  else:
    return akward_array

def expo_func(x, D, tau, t):
  return lambda x,D,tau,t: D*(np.exp(-(x-t)/tau))

expo=lambda x,D,tau,t: D*(np.exp(-(x-t)/tau))

#PLOTSPATH="./plots"
#PLOTSPATH="./plots/test_morefiles_run4453_threshold_noQcut"
#PLOTSPATH="./plots/run4496"
#PLOTSPATH="./plots/run4499_cb_lessthan0p3"
#PLOTSPATH="./plots/run4499_trigQcut"
#PLOTSPATH="./plots/run4707"
PLOTSPATH="./plots/run4708"

DATAPATH="./data"
#file = uproot.open("/Users/giancaceresvera/Desktop/ANNIE/AmBe_source/AmBe_DAQ/data_threshold_for_test/files/allfiles.root")
#file = uproot.open("/Users/giancaceresvera/Desktop/ANNIE/AmBe_source/AmBe_DAQ/data_threshold_for_test/run4496/PhaseIITreeMaker_run4496.root")
#file = uproot.open("/Users/giancaceresvera/Desktop/ANNIE/AmBe_source/AmBe_DAQ/data_threshold_for_test/run4499/PhaseIITreeMaker_run4499.root")
#file = uproot.open("/Users/giancaceresvera/Desktop/ANNIE/AmBe_source/AmBe_DAQ/data_threshold_for_test/run4707/PhaseIITreeMaker_run4707.root")
file = uproot.open("/Users/giancaceresvera/Desktop/ANNIE/AmBe_source/AmBe_DAQ/data_threshold_for_test/run4708/PhaseIITreeMaker_run4708.root")

listofkeys=file.keys()
dictofclasses=file.keys()
classnamesoffile=file.classnames()

print("list of keys in file  ",listofkeys)
print("dict of classes  ",dictofclasses)
print("classnames of file",classnamesoffile)
print("arraw of hitQ",file["phaseIITriggerTree/SiPMhitQ"].array())
print(file["phaseIITriggerTree"])

tree = file["phaseIITriggerTree"]
tree_clusters = file["phaseIITankClusterTree"]

print("tree keys ", tree.keys())
print("tree ", tree['SiPMhitQ'])

branches= tree.arrays()
branches_clusters=tree_clusters.arrays()

print('branches_clusters=  ',tree_clusters.arrays())
print('branches_clusters=  ',branches_clusters)
  
for i in tree_clusters.keys():
  if branches_clusters[i].ndim > 1 : 
    print("multidimensional : ", i )
  else :
    print("unidimensional: " , i )

#multidimensional arrays
hitX             = ak.flatten(branches_clusters['hitX'])
hitY             = ak.flatten(branches_clusters['hitY'])
hitZ             = ak.flatten(branches_clusters['hitZ'])
hitT             = ak.flatten(branches_clusters['hitT'])
hitQ             = ak.flatten(branches_clusters['hitQ'])
hitPE            = ak.flatten(branches_clusters['hitPE'])
hitType          = ak.flatten(branches_clusters['hitType'])
hitDetID         = ak.flatten(branches_clusters['hitDetID'])
hitChankey       = ak.flatten(branches_clusters['hitChankey'])
hitChankeyMC     = ak.flatten(branches_clusters['hitChankeyMC'])
SiPMhitQ         = ak.flatten(branches_clusters['SiPMhitQ'])
SiPMhitT         = ak.flatten(branches_clusters['SiPMhitT'])
SiPMhitAmplitude = ak.flatten(branches_clusters['SiPMhitAmplitude'])
SiPMNum          = ak.flatten(branches_clusters['SiPMNum'])
ADCSamples       = ak.flatten(branches_clusters['ADCSamples'])
ADCChankeys      = ak.flatten(branches_clusters['ADCChankeys'])

#unidimensional arrays
runNumber             = flat_array(branches_clusters['runNumber'])             
subrunNumber          = flat_array(branches_clusters['subrunNumber'])  
runType               = flat_array(branches_clusters['runType'])  
startTime             = flat_array(branches_clusters['startTime'])  
eventNumber           = flat_array(branches_clusters['eventNumber'])  
eventTimeTank         = flat_array(branches_clusters['eventTimeTank'])  
clusterNumber         = flat_array(branches_clusters['clusterNumber'])  
clusterTime           = flat_array(branches_clusters['clusterTime'])  
clusterCharge         = flat_array(branches_clusters['clusterCharge'])  
clusterPE             = flat_array(branches_clusters['clusterPE'])  
clusterMaxPE          = flat_array(branches_clusters['clusterMaxPE'])  
clusterChargePointX   = flat_array(branches_clusters['clusterChargePointX'])  
clusterChargePointY   = flat_array(branches_clusters['clusterChargePointY'])  
clusterChargePointZ   = flat_array(branches_clusters['clusterChargePointZ'])  
clusterChargeBalance  = flat_array(branches_clusters['clusterChargeBalance'])  
clusterHits           = flat_array(branches_clusters['clusterHits'])  
trigword              = flat_array(branches_clusters['trigword'])  
TankMRDCoinc          = flat_array(branches_clusters['TankMRDCoinc'])  
NoVeto                = flat_array(branches_clusters['NoVeto'])  
Extended              = flat_array(branches_clusters['Extended'])  
beam_pot              = flat_array(branches_clusters['beam_pot'])  
beam_ok               = flat_array(branches_clusters['beam_ok'])  
SiPM1NPulses          = flat_array(branches_clusters['SiPM1NPulses'])  
SiPM2NPulses          = flat_array(branches_clusters['SiPM2NPulses'])  

print("len SiPMhitT ", len(SiPMhitT)) 
print("len SiPMNum ", len(SiPMNum)) 
print("len SiPM1NPulses ", len(SiPM1NPulses)) 

for i in range(0,2):
  print(" SiPM1NPulses[i] ", (SiPM1NPulses[i])) 
  print(" SiPMhitT[i] ", (SiPMhitT[i])) 

#mask_mediumtrigcharge = (SiPMhitQ>0.1) & (SiPMhitQ < 0.4)
#print("len mask_mediumtrigcharge ", len(mask_mediumtrigcharge)) 
#for i in range (0,5):
#  print(" mask_mediumtrigcharge[i] ", (mask_mediumtrigcharge[i])) 

#exit(1)

#cuts
mask_lowtrigcharge = SiPMhitQ<0.05
mask_mediumtrigcharge = (SiPMhitQ>0.1) & (SiPMhitQ < 0.4)
mask_lowclustertime = clusterTime<1000
#mask_clusterneartrig =(abs(clusterTime - SiPMhitT) < 50) 
mask_prompttrig = SiPMhitT < 4000 
mask_lowclustercb = clusterChargeBalance < 0.3
mask_delaycluster = clusterTime > 2000
mask_clusterPE = ((clusterPE > 15) & (clusterPE < 150)) 
mask_min_clusterPE = clusterPE > 60 
mask_max_clusterPE = clusterPE < 70 
mask_hightrighamp = SiPMhitAmplitude > 0.05
mask_trigword = (trigword ==15)

np_clcb = np.array(clusterChargeBalance)
np_clQ = np.array(clusterCharge)
np_clPE = np.array(clusterPE)
np_clmaxPE = np.array(clusterMaxPE)
np_clT = np.array(clusterTime)
np_clNum = np.array(clusterNumber)


plt.hist(np_clQ, density=False, alpha=0.8, bins=100, range =(0,250) )
plt.xlabel("Cluster Charge (nC)", fontsize=20)
plt.savefig('%s/cluster_charge.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist(np_clQ, density=False, alpha=0.8, bins=100, range =(0,250) )
plt.xlabel("Cluster Charge (nC)", fontsize=20)
plt.yscale('log')
plt.savefig('%s/cluster_charge_log.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist(np_clPE, density=False, alpha=0.8, bins=90 )
plt.xlabel("Total PE", fontsize=20)
#plt.yscale('log')
plt.savefig('%s/PE.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist(np_clPE, density=False, alpha=0.8, bins=90 )
plt.xlabel("Total PE", fontsize=20)
plt.yscale('log')
plt.savefig('%s/PE_log.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist(np_clQ, density=False, alpha=0.8, bins=100, range =(0,250) )
plt.hist(np_clQ[mask_delaycluster], density=False, alpha=0.8, bins=100, range =(0,250), label="delayed cluster" )
plt.xlabel("Cluster Charge (nC)", fontsize=20)
plt.savefig('%s/cluster_charge_and_delaycluster.png'%(PLOTSPATH))
plt.yscale('log')
plt.legend(loc='upper right')
#plt.show()
plt.close()


plt.hist2d( np_clQ, np_clcb, bins=(20,30),cmap=plt.get_cmap('viridis'))
plt.xlabel("Charge (nC)", fontsize=20)
plt.ylabel("Charge Balance", fontsize=20)
plt.colorbar()
plt.savefig('%s/q_cb.png'%(PLOTSPATH))
#plt.show()
plt.close()

#log bar
plt.hist2d( np_clPE, np_clcb, bins=(20,30),norm=mpl.colors.LogNorm(), cmap=plt.get_cmap('viridis'))
plt.xlabel("Total PE", fontsize=20)
plt.ylabel("Charge Balance", fontsize=20)
plt.colorbar()
plt.savefig('%s/PE_cb.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist2d( np_clPE, np_clcb, bins=(20,30), range=[[0,200],[0,1]] , cmap=plt.get_cmap('viridis'))
plt.xlabel("Total PE", fontsize=20)
plt.ylabel("Charge Balance", fontsize=20)
#plt.xlim([0,200]) # only zoom
plt.colorbar()
plt.savefig('%s/PE_cb_zoom.png'%(PLOTSPATH))
#plt.show()
plt.close()


plt.hist2d(np_clT, np_clcb ,bins=(20,30),cmap=plt.get_cmap('viridis'))
plt.xlabel("Cluster Time", fontsize=20)
plt.ylabel("Charge Balance", fontsize=20)
plt.colorbar()
plt.savefig('%s/t_cb.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist2d(np_clT, np_clQ ,bins=(20,30),cmap=plt.get_cmap('viridis'))
plt.xlabel("Cluster Time", fontsize=20)
plt.ylabel("Cluster Charge", fontsize=20)
plt.colorbar()
plt.savefig('%s/t_q.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist(np_clT, density=False, alpha=0.8, bins=100 )
plt.xlabel("Cluster Time (ns)", fontsize=20)
plt.savefig('%s/cluster_time.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist(np_clT, density=False, alpha=0.8, bins=100 )
plt.xlabel("Cluster Time (ns)", fontsize=20)
plt.yscale('log')
plt.savefig('%s/cluster_time_log.png'%(PLOTSPATH))
#plt.show()
plt.close()

###CUTS###
#cut_dummy_init= (mask_delaycluster  & mask_max_clusterPE)
#cut_dummy= (mask_delaycluster  & mask_max_clusterPE & mask_lowclustercb & mask_mediumtrigcharge)
#cut_dummy= (mask_delaycluster  & mask_max_clusterPE & mask_lowclustercb )

cut_dummy= ()
#cut_dummy= (mask_delaycluster  & mask_max_clusterPE & mask_lowclustercb  & mask_trigword)
#cut_dummy= (mask_delaycluster  & mask_max_clusterPE &  mask_trigword  )
#cut_dummy= (mask_delaycluster  & mask_max_clusterPE & mask_lowclustercb & mask_trigword & mask_mediumtrigcharge )
#cut_dummy= (mask_delaycluster  & mask_max_clusterPE & mask_trigword & mask_mediumtrigcharge )


####cut_dummy_init plots

#plt.hist2d( np_clPE[cut_dummy_init], np_clcb[cut_dummy_init], bins=(20,30), cmap=plt.get_cmap('viridis'), label='cuts')
#plt.xlabel("Total PE", fontsize=20)
#plt.ylabel("Charge Balance", fontsize=20)
##plt.title("cuts")
##plt.xlim([0,100])
#plt.colorbar()
#plt.savefig('%s/PE_cb_testcuts_init.png'%(PLOTSPATH))
##plt.show()
#plt.close()
#
#plt.hist2d(np_clT[cut_dummy_init], np_clcb[cut_dummy_init] ,bins=(20,30),cmap=plt.get_cmap('viridis'), label='cuts')
#plt.xlabel("Cluster Time", fontsize=20)
#plt.ylabel("Charge Balance", fontsize=20)
##plt.title("cuts")
#plt.colorbar()
#plt.savefig('%s/t_cb_testcuts_init.png'%(PLOTSPATH))
##plt.show()
#plt.close()
#
#plt.hist(np_clT[cut_dummy_init], density=False, alpha=0.8, bins=100 )
#plt.xlabel("Cluster Time (ns)", fontsize=20)
#plt.savefig('%s/cluster_time_testcuts_init.png'%(PLOTSPATH))
##plt.show()
#plt.close()

########

plt.hist(np_clT[cut_dummy], density=False, alpha=0.8, bins=100 )
plt.xlabel("Cluster Time (ns)", fontsize=20)
plt.savefig('%s/cluster_time_testcuts.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist(np_clT[cut_dummy], density=False, alpha=0.8, bins=100 )
plt.xlabel("Cluster Time (ns)", fontsize=20)
plt.yscale('log')
plt.savefig('%s/cluster_time_testcuts_log.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist(np_clPE[cut_dummy], density=False, alpha=0.8, bins=90 )
plt.xlabel("Total PE", fontsize=20)
#plt.yscale('log')
plt.savefig('%s/PE_testcuts.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist(np_clPE[cut_dummy], density=True, alpha=0.8, bins=90 )
plt.xlabel("Total PE", fontsize=20)
#plt.yscale('log')
plt.savefig('%s/PE_norm_testcuts.png'%(PLOTSPATH))
#plt.show()
plt.close()

###Multiplicity
num_entries=len(np_clNum[cut_dummy])
plt.hist(np_clNum[cut_dummy], density=False, alpha=0.8, bins=6, range =(0,6),label=('entries=%i'%(num_entries)))
plt.xlabel("Number of clusters", fontsize=20)
plt.yscale('log')
plt.legend(loc = 'upper right')
plt.savefig('%s/cluster_num_testcuts.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist(np_clNum[cut_dummy], density=True, alpha=0.8, bins=6, range =(0,6),label=('entries=%i'%(num_entries)))
plt.xlabel("Number of clusters", fontsize=20)
#plt.yscale('log')
plt.legend(loc = 'upper right')
plt.savefig('%s/cluster_num_norm_testcuts.png'%(PLOTSPATH))
#plt.show()
plt.close()



cut_dummy_all= (mask_delaycluster  & mask_max_clusterPE & mask_trigword & mask_lowclustercb & mask_mediumtrigcharge )
cut_dummy_noQ= (mask_delaycluster  & mask_max_clusterPE & mask_trigword & mask_lowclustercb )
cut_dummy_noQ_noCB= (mask_delaycluster  & mask_max_clusterPE & mask_trigword)

if save_multiplicity:
  np.savetxt('%s/multiplicity_all.txt'%(PLOTSPATH),np_clNum[cut_dummy_all],fmt='%.8f',newline=os.linesep)
  np.savetxt('%s/multiplicity_noQ.txt'%(PLOTSPATH),np_clNum[cut_dummy_noQ],fmt='%.8f',newline=os.linesep)
  np.savetxt('%s/multiplicity_noQ_noCB.txt'%(PLOTSPATH),np_clNum[cut_dummy_noQ_noCB],fmt='%.8f',newline=os.linesep)
##########

plt.hist2d( np_clQ[cut_dummy], np_clcb[cut_dummy], bins=(20,30),cmap=plt.get_cmap('viridis'), label='cuts')
plt.xlabel("Cluster Charge (nC)", fontsize=20)
plt.ylabel("Charge Balance", fontsize=20)
#plt.title("cuts")
plt.colorbar()
plt.savefig('%s/q_cb_testcuts.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist2d( np_clPE[cut_dummy], np_clcb[cut_dummy], bins=(20,30), cmap=plt.get_cmap('viridis'), label='cuts')
plt.xlabel("Total PE", fontsize=20)
plt.ylabel("Charge Balance", fontsize=20)
#plt.title("cuts")
#plt.xlim([0,100])
plt.colorbar()
plt.savefig('%s/PE_cb_testcuts.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist2d( np_clPE[cut_dummy], np_clcb[cut_dummy], bins=(20,30), norm=mpl.colors.LogNorm(),cmap=plt.get_cmap('viridis'), label='cuts')
plt.xlabel("Total PE", fontsize=20)
plt.ylabel("Charge Balance", fontsize=20)
#plt.title("cuts")
#plt.xlim([0,100])
plt.colorbar()
plt.savefig('%s/PE_cb_testcuts_log.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist2d( np_clPE[cut_dummy], np_clcb[cut_dummy], bins=(20,30),range=[[0,200],[0,1]], cmap=plt.get_cmap('viridis'), label='cuts')
plt.xlabel("Total PE", fontsize=20)
plt.ylabel("Charge Balance", fontsize=20)
#plt.title("cuts")
#plt.xlim([0,100])
plt.colorbar()
plt.savefig('%s/PE_cb_zoom_testcuts.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist2d(np_clT[cut_dummy], np_clcb[cut_dummy] ,bins=(20,30),cmap=plt.get_cmap('viridis'), label='cuts')
plt.xlabel("Cluster Time", fontsize=20)
plt.ylabel("Charge Balance", fontsize=20)
#plt.title("cuts")
plt.colorbar()
plt.savefig('%s/t_cb_testcuts.png'%(PLOTSPATH))
#plt.show()
plt.close()

plt.hist2d(np_clT[cut_dummy], np_clQ[cut_dummy] ,bins=(20,30),cmap=plt.get_cmap('viridis'), label='cuts')
plt.xlabel("Cluster Time", fontsize=20)
plt.ylabel("Cluster Charge", fontsize=20)
#plt.title("cuts")
plt.colorbar()
plt.savefig('%s/t_q_testcuts.png'%(PLOTSPATH))
#plt.show()
plt.close()

####capture time

neutronlifetime =  np_clT[cut_dummy]  
print ('entries of np_clT = ', len(np_clT))
entries = len(neutronlifetime)
mean = np.mean(neutronlifetime)
std_dev = np.std(neutronlifetime)
n, xedges, patches = plt.hist(neutronlifetime, density=False, alpha=0.8, bins=30 )
bin_centers=(xedges[:-1] + xedges[1:])/2
print('bin_centers = ', bin_centers)
print('n = ', n)
print('len bin_centers = ', len(bin_centers))
print('len n = ', len(n))

#ignoring first bins for fitting

#erasing firt element of array
#n=n[6:]
#bin_centers=bin_centers[6:]

n=n[4:]
bin_centers=bin_centers[4:]
print('after erasing first element')
print('bin_centers = ', bin_centers)
print('n = ', n)

#init_params = [500.0, 500.0, 500.0]
init_params = [1000.0, 1000.0, 1000.0]
#popt, pcov = scp.curve_fit(expo_func, bin_centers, n ,p0=init_params, maxfev=6000)
popt, pcov = scp.curve_fit(expo, bin_centers, n, p0=init_params )
print('popt = ', popt)
print('pcov = ', pcov)
print('entries ', entries)
print('mean ', mean)
print('std dev ', std_dev)
#legend_elements = ['entries= %.2f'%entries, 'lifetime= %.2f'%popt[1] ]
#plt.legend(legend_elements, loc = 'upper right')
plt.plot(bin_centers,expo(bin_centers,popt[0], popt[1], popt[2]), label=(r"$\tau$=%.2f ns"%popt[1]))
plt.legend(loc = 'upper right')
plt.xlabel("Cluster time", fontsize=20)
plt.savefig('%s/capturetime.png'%(PLOTSPATH))
#plt.yscale('log')
#plt.show()
plt.close()

####creating old gains cluster hists
if not plot_gaincomp:
  np.savetxt("oldgains_PEhist.txt",np_clPE,fmt='%.8f',newline=os.linesep)
  np.savetxt("oldgains_Qhist.txt",np_clQ,fmt='%.8f',newline=os.linesep)
  np.savetxt("oldgains_CBhist.txt",np_clcb,fmt='%.8f',newline=os.linesep)
  np.savetxt("oldgains_Thist.txt",np_clT,fmt='%.8f',newline=os.linesep)
  np.savetxt("oldgains_PEhist_cuts.txt",np_clPE[cut_dummy],fmt='%.8f',newline=os.linesep)

if plot_gaincomp:
  ####reading old gains cluster hists
  oldgains_clPE=np.genfromtxt('oldgains_PEhist.txt')
  oldgains_clQ=np.genfromtxt('oldgains_Qhist.txt')
  oldgains_clcb=np.genfromtxt('oldgains_CBhist.txt')
  oldgains_clT=np.genfromtxt('oldgains_Thist.txt')
  oldgains_clPE_cuts=np.genfromtxt('oldgains_PEhist_cuts.txt')
  
  ####plotting comparisons 
  plt.hist(np_clPE, density=True, alpha=0.8, bins=90, label='new gains' )
  plt.hist(oldgains_clPE, density=True, alpha=0.8, histtype='step', bins=90,color='red', label='old gains')
  plt.xlabel("Total PE", fontsize=20)
  plt.yscale('log')
  plt.legend(loc = 'upper right')
  plt.savefig('%s/PE_norm_gaincomp.png'%(PLOTSPATH))
  #plt.show()
  plt.close()
  
  plt.hist(np_clQ, density=True, alpha=0.8, bins=90, label='new gains' )
  plt.hist(oldgains_clQ, density=True, alpha=0.8, histtype='step', bins=90,color='red', label='old gains')
  plt.xlabel("Charge", fontsize=20)
  plt.yscale('log')
  plt.legend(loc = 'upper right')
  plt.savefig('%s/Q_norm_gaincomp.png'%(PLOTSPATH))
  #plt.show()
  plt.close()
  
  plt.hist(np_clcb, density=True, alpha=0.8, bins=90, label='new gains' )
  plt.hist(oldgains_clcb, density=True, alpha=0.8, histtype='step', bins=90,color='red', label='old gains')
  plt.xlabel("Charge Balance", fontsize=20)
  plt.yscale('log')
  plt.legend(loc = 'upper right')
  plt.savefig('%s/cb_norm_gaincomp.png'%(PLOTSPATH))
  #plt.show()
  plt.close()
  
  plt.hist(np_clT, density=True, alpha=0.8, bins=90, label='new gains' )
  plt.hist(oldgains_clT, density=True, alpha=0.8, histtype='step', bins=90,color='red', label='old gains')
  plt.xlabel("Cluster Time", fontsize=20)
  plt.yscale('log')
  plt.legend(loc = 'upper right')
  plt.savefig('%s/T_norm_gaincomp.png'%(PLOTSPATH))
  #plt.show()
  plt.close()
  
  plt.hist(np_clPE[cut_dummy], density=True, alpha=0.8, bins=90, label='new gains' )
  plt.hist(oldgains_clPE_cuts, density=True, alpha=0.8, histtype='step', bins=90,color='red', label='old gains')
  plt.xlabel("Total PE", fontsize=20)
  #plt.yscale('log')
  plt.legend(loc = 'upper right')
  plt.savefig('%s/PE_norm_testcuts_gaincomp.png'%(PLOTSPATH))
  #plt.show()
  plt.close()
########################

