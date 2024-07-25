import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.stats import poisson
import scipy.optimize as opt
from scipy.interpolate import griddata
import scipy.interpolate 
import math

PLOTSPATH="./plots/bkg_models/testing"
#DATAPATH="./data"

#print("random number")
#print(random.uniform(0,1))
plot_model_hists=False
plot_model_hists_2=False
plot_chi2_profile=False

eff_arr =[] 
rate_arr =[] 
chi2_arr =[] 

data_values_2=[2.026e+03, 2.491e+03, 5.620e+02, 6.800e+01, 4.000e+00, 0.000e+00, 0.000e+00]
data_values_error_2raw = [math.sqrt(element) for element in data_values_2]
data_values_error_2 = [0.000001 if elem==0.0 else elem for elem in data_values_error_2raw]

def multiplicity(eff,lambda_parameter):
  corr_draw=random.uniform(0,1)
  bkg_draw=random.uniform(0,1)
  if(corr_draw<=eff):
    corr_n=1
  elif(corr_draw>eff):
    corr_n=0
  y = np.arange(0, 5)
  cumulative_probabilities = [poisson.cdf(k, lambda_parameter) for k in y]
  #cumulative_probabilities = [poisson.pmf(k, lambda_parameter) for k in y]
  dum_list=[]
  for prob in cumulative_probabilities:
    if bkg_draw>prob:
      dum_index=cumulative_probabilities.index(prob)
      dum_list.append(dum_index)
  if dum_list:
    bkg_n=dum_list[-1]+1 
  elif not dum_list:
    bkg_n=0
  M=corr_n + bkg_n
  return M

def chisqfunc(x):
    model =[]
    draws = np.arange(0,1000)
    model = [multiplicity(x[0],x[1]) for a in draws]
    model_v,model_b_e, model_p = plt.hist(model, range = (0,7), bins=7)
    plt.close()
    model_fl = [float(elem) for elem in model_v]
    #print('sums sum(data_values) and sum(model_fl) and ratio')
    #print(sum(data_values))
    #print(sum(model_fl))
    #print(sum(model_v))
    #print(sum(data_values)/sum(model_fl))
    model_norm = [elem * (sum(data_values)/sum(model_fl)) for elem in model_fl]
    #print('elements of data values')
    #print(data_values)
    #print('elements of data values 2')
    #print(data_values_2)
    #print('elements of model')
    #print(model_v)
    #print('elements of model norm')
    #print(model_norm)
    #plt.savefig('%s/model_%s.png'%(PLOTSPATH,i))
    if(plot_model_hists):
      #plt.bar(np.arange(len(model_norm)), model_norm , color='green', ls='dashed')
      #plt.xlabel("Multiplicity",fontsize=20)
      #plt.title('model normalized')
      #plt.show()
      #plt.close()
      #plt.bar(np.arange(len(data_values)),data_values, alpha = 0.3, color ='red', ls='dotted')
      #plt.xlabel("Multiplicity",fontsize=20)
      #plt.title('data')
      #plt.show()
      #plt.close()
      #plt.bar(np.arange(len(model_norm)), model_norm , align='edge',edgecolor='green', color = None, width =1.0)
      plt.step(np.arange(len(model_norm)), model_norm , where="post", color = 'blue' )
      plt.step(np.arange(len(data_values)), data_values,  where="post", color = "green")
      plt.ylim(0.0,5000)
      plt.xlabel("Multiplicity",fontsize=20)
      plt.title('Data and Normalized Model')
      plt.show()
      plt.close()
      #plt.step(np.arange(len(data_values)), data_values,  where="post", color = "green")
      #print(data_values)
      #plt.xlabel("Multiplicity",fontsize=20)
      #plt.title('data')
      #plt.show()
      #plt.close()
    #hisq = np.sum(((data_values_np - model_v_np)/(math.sqrt(data_values_np)))**2)
    chisq = np.sum(((data_values - model_norm)/data_values_error)**2)
    #chisq = np.sum(((data_values - model_norm)/0.002)**2)
    #chisq = np.sum(((data_values - model_v)/0.002)**2)
    #chisq = np.sum(((data_values_2 - model_v)/0.002)**2)
    #print('model_v')
    #print(model_v)
    #print('model_v_norm')
    #print(model_v_norm)
    #print('data_values')
    #print(data_values)
    #print('data_values_error')
    #print(data_values_error)
    #print("difference")
    #print(data_values_norm - model_v_norm)
    #print('ef and rate')
    #print(x[0])
    #print(x[1])
    #print('current X2')
    #print(chisq)
    eff_arr.append(x[0])
    rate_arr.append(x[1])
    chi2_arr.append(chisq)
    return chisq

def chisqfunc_2(x):
    model =[]
    draws = np.arange(0,1000)
    model = [multiplicity(x[0],x[1]) for a in draws]
    model_v,model_b_e, model_p = plt.hist(model, range = (0,7), bins=7)
    plt.close()
    model_fl = [float(elem) for elem in model_v]
    model_norm = [elem * (sum(data_values_2)/sum(model_fl)) for elem in model_fl]
    if(plot_model_hists_2):
      plt.bar(np.arange(len(model_norm)), model_norm , color='green', ls='dashed')
      plt.xlabel("Multiplicity",fontsize=20)
      plt.title('model normalized')
      plt.show()
      plt.close()
      plt.bar(np.arange(len(data_values_2)),data_values_2, alpha = 0.3, color ='red', ls='dotted')
      plt.xlabel("Multiplicity",fontsize=20)
      plt.title('data 2')
      plt.show()
      plt.close()
    model_norm=np.array(model_norm)
    #hisq = np.sum(((data_values_np - model_v_np)/(math.sqrt(data_values_np)))**2)
    #chisq = np.sum(((data_values_2 - model_norm)/data_values_error_2)**2)
    #chisq = np.sum(((data_values_2 - model_norm)/0.002)**2)
    chisq = np.sum(((data_values_2 - model_norm)/data_values_error)**2)
    #chisq = np.sum(((data_values_2 - model_v)/0.002)**2)
    eff_arr.append(x[0])
    rate_arr.append(x[1])
    chi2_arr.append(chisq)
    return chisq

#######data
data_multi_hist=np.genfromtxt('multiplicity_noQ.txt')
data_values, bins_edges,patches=plt.hist(data_multi_hist ,density=False,alpha = 0.4, range=(0,7), color = 'r',bins=7)
data_values_error_raw = [math.sqrt(elem) for elem in data_values]
data_values_error = [0.000001 if elem==0.0 else elem for elem in data_values_error_raw]

plt.title('Multiplicity from Data')
plt.xlabel('number of neutron captures')
plt.ylabel('events')
#plt.yscale('log')
plt.grid(True)
#plt.legend()
plt.savefig('%s/data_multi.png'%(PLOTSPATH))
plt.show()
plt.close()

#######data 2 for testing the Chi2
#x_values_2=np.arange(7)
#plt.bar(x_values_2, data_values_2 , color = 'r', align='edge', width=1.0)
#plt.title('Multiplicity from Data 2')
#plt.xlabel('number of neutron captures')
#plt.ylabel('events')
##plt.yscale('log')
#plt.grid(True)
##plt.legend()
#plt.savefig('%s/data_multi_2.png'%(PLOTSPATH))
#plt.show()
#plt.close()

################

#x0 = np.array([0.64,0.07])
x0 = np.array([0.64,0.07])
#result = opt.minimize(chisqfunc, x0, method = 'BFGS',  callback=store_params)
#result = opt.minimize(chisqfunc, x0, method='BFGS', bounds=((0,1),(0.005,0.09)))
#result = opt.minimize(chisqfunc, x0, method='Newton-CG', bounds=((0,1),(0.005,0.09))) #needs jacobian
#result = opt.minimize(chisqfunc, x0, method='Powell', bounds=((0,1),(0.005,0.09)))   #nan
#result = opt.minimize(chisqfunc, x0, method='CG', bounds=((0,1),(0.005,0.09))) #nan
#result = opt.minimize(chisqfunc, x0, method='Nelder-Mead', bounds=((0,1),(0.004,0.1)))
#result = opt.minimize(chisqfunc, x0, method='Nelder-Mead', bounds=((0,1),(0.0,0.1)))
#result = opt.minimize(chisqfunc, x0, method='Nelder-Mead', bounds=((0,1),(0.0,0.6)))
#print ("result")
#print(result)

result = opt.minimize(chisqfunc, x0, method='Nelder-Mead', bounds=((0,1),(0.0,0.1)))
#result = opt.minimize(chisqfunc, x0, bounds=((0,1),(0.0,0.1)))
print ("result")
print(result)
print(result.x)
print(result.fun)

#result = opt.minimize(chisqfunc, x0, bounds=((0.001,1),(0.001,0.1)))
#print ("result")
#print(result)

#print("len chi2")
#print(len(chi2_arr))
##print("len eff array call back")
##print(len(eff_values))
#print("len eff array inside function")
#print(len(eff_arr))
##print("len rate array call back")
##print(len(lambda_values))
#print("len rate array inside function")
#print(len(rate_arr))

################
#for scatter plot
#using values from the chi2 function
chi2_ = np.asarray(chi2_arr)
eff_  = np.asarray(eff_arr)
rate_  = np.asarray(rate_arr)
print('np lens')
print(len(eff_))
print(len(rate_))
print(len(chi2_))
print('dimensions')
print(eff_.ndim)
print(rate_.ndim)
print(chi2_.ndim)
print('min chi2')
print(min(chi2_))

print('chi2')
print(chi2_)


plt.scatter(eff_,rate_, c=chi2_, s=200, alpha=0.5 ) # s is a size of marker 
plt.scatter(result.x[0], result.x[1], color='red', marker='x', label='Optimal Params')
plt.title("Scatter Plot - optimization values")
plt.legend()
plt.colorbar()
plt.savefig('%s/scatter_plot.png'%(PLOTSPATH))
plt.show()
plt.close()

################

# plotting 2D plot of the chi2 function
# Define range for shift and scale parameters
if plot_chi2_profile:
  shift_range = np.linspace(0, 1, 10)
  scale_range = np.linspace(0.0, 0.1, 10)
  shift_grid, scale_grid = np.meshgrid(shift_range, scale_range)
  
  # Compute chi-square values for each combination of shift and scale
  chi2_values = np.zeros_like(shift_grid)
  for i in range(len(shift_range)):
      for j in range(len(scale_range)):
          params = [shift_grid[i, j], scale_grid[i, j]]
          chi2_values[i, j] = chisqfunc(params)
  
  plt.figure(figsize=(10, 6))
  #contour = plt.contourf(shift_grid, scale_grid, chi2_values, levels=50, cmap='viridis')
  #plt.colorbar(contour, label='Chi-Square Value')
  plt.imshow(chi2_values, extent=(shift_range.min(), shift_range.max(), scale_range.min(), scale_range.max()), aspect='auto', origin='lower')
  plt.colorbar(label='Chi-Square Value')
  plt.scatter(result.x[0], result.x[1], color='red', marker='x', label='Optimal Params')
  plt.xlabel('Efficiency')
  plt.ylabel('Background rate')
  plt.title('Chi-Square Contour Plot - Chi2 function evaluated for a range of values')
  plt.legend()
  plt.savefig('%s/chi2_profile.png'%(PLOTSPATH))
  plt.show()
  plt.close()


####the unique function returns the sorted unique elements of an array
#rate_=np.unique(rate_)
#eff_=np.unique(eff_)
#X,Y=np.meshgrid(rate_,eff_)
#
#Z=chi2_.reshape(len(eff_),len(rate_))
#plt.pcolormesh(X,Y,Z)
#plt.show()


####reshape doens't work because I need to know the 2 axis lenghts in advance before plotting 
#num_bins = len(rate_arr)
#plt.figure(figsize=(8, 6))
#plt.imshow(chi2_.reshape((num_bins, num_bins)), cmap='viridis', origin='lower', extent=(eff_.min(), eff_.max(), rate_.min(), rate_.max()))
#plt.colorbar(label='Intensity')
#plt.xlabel('X-axis')
#plt.ylabel('Y-axis')
#plt.title('Intensity Plot')
#plt.show()

#x0 = np.array([0.64,0.07])
#result = opt.minimize(chisqfunc, x0, method='Nelder-Mead', bounds=([0,0],[1,1]))

#print("result")
#print (result)
#print("len eff array")
#print(len(eff_arr))
#print("len rate array")
#print(len(rate_arr))
#print("len chi2")
#print(len(chi2_arr))
#print("eff array")
#print((eff_arr))
#print("rate array")
#print((rate_arr))
#print("chi2")
#print((chi2_arr))
#
#y=np.array(eff_arr)
#x=np.array(rate_arr)
#z=np.array(chi2_arr)
##print("min numpy chi2")
##print(np.amin(chi2_arr))
#print("result.x")
#print(result.x)
#print("result.fun")
#print(result.fun)
##xi = np.linspace(np.min(x),np.max(x),100)
##yi = np.linspace(np.min(x),np.max(x),100)
#xi = np.linspace(x.min(),x.max(),100)
#yi = np.linspace(y.min(),y.max(),100)
##zi=griddata((x,y),z,(xi,yi),method='linear')
#xi,yi = np.meshgrid(xi,yi)
#zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='linear')
#plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',extent=[x.min(), x.max(), y.min(), y.max()])
#plt.colorbar()
#plt.show()
#plt.close()

#################
#times = np.arange(0,10)
#multi_hist= [multiplicity(0.64,0.07) for a in times]
#
#plt.hist(multi_hist, range = (0,7), bins=7)
#plt.title('Multiplicity distribution model')
#plt.xlabel('number of neutron captures')
#plt.ylabel('events')
#plt.yscale('log')
#plt.grid(True)
#plt.legend()
#plt.savefig('%s/multi.png'%(PLOTSPATH))
#plt.show()
#
#################
## Define the parameters of the Poisson distribution
#lambda_parameter = 0.007  # Adjust this value to your desired rate
#
## Generate a range of values for the x-axis
#x = np.arange(0, 5)
#
## Calculate the cumulative Poisson distribution for each value in x
#cumulative_probabilities_another = [poisson.cdf(k, lambda_parameter) for k in x]
#cumulative_probabilities_greaterthan = [1-(poisson.cdf(k, lambda_parameter)) for k in x]
#pmf_values = [poisson.pmf(k, lambda_parameter) for k in x]

## Plot poisson distribution
#plt.plot(x, pmf_values , marker='o', linestyle='-', color='b', label=f'λ = {lambda_parameter}')
#plt.title('Poisson Distribution')
#plt.xlabel('k')
#plt.ylabel('P(X = k)')
#plt.grid(True)
#plt.legend()
#plt.savefig('%s/pmf.png'%(PLOTSPATH))
#plt.show()
#
## Plot the cumulative distribution
#plt.plot(x, cumulative_probabilities_another, marker='o', linestyle='-', color='b', label=f'λ = {lambda_parameter}')
#plt.title('Cumulative Poisson Distribution')
#plt.xlabel('k')
#plt.ylabel('P(X ≤ k)')
#plt.grid(True)
#plt.legend()
#plt.savefig('%s/cdf.png'%(PLOTSPATH))
#plt.show()
#
## Plot the cumulative distribution
#plt.plot(x, cumulative_probabilities_greaterthan, marker='o', linestyle='-', color='b', label=f'λ = {lambda_parameter}')
#plt.title('Cumulative Poisson Distribution')
#plt.xlabel('k')
#plt.ylabel('P(X >= k)')
#plt.grid(True)
#plt.legend()
#plt.savefig('%s/cdf_gt.png'%(PLOTSPATH))
#plt.show()
