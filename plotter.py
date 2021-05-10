# -*- coding: utf-8 -*-
"""
Created on Sun May  9 15:26:14 2021

@author: mstan
"""
# %%
import numpy as np

import matplotlib.pyplot as plt
import spym
plt.close('all')
central_wavelength=299700.0/374.7
EPP = 0.001
Length=1.0
beam_dia= 12.0
n2=2.5*10**(-20)
beta2=36.0
beta3=0.0
beta4=0.0

data = np.loadtxt(r'initial.txt',skiprows=1,delimiter='\t')
pulse = spym.createPulseSpecNLO(data[:,0], data[:,1],data[:,2],EPP=EPP,central_wavelength=central_wavelength)

fiber = spym.createFiberNLO(central_wavelength=central_wavelength,Length=Length, beam_dia=beam_dia,n2=n2,beta2=beta2,beta3=beta3,beta4=beta4)
_,_,_,pulse_pynlo = spym.runNLOSim(pulse, fiber,Steps=10,Raman=False,Steep=False)

data = np.loadtxt(r'spm.txt',skiprows=1,delimiter='\t')
pulse_spm = spym.createPulseSpecNLO(data[:,0], data[:,1],data[:,2],EPP=0.001)




fig, ax = plt.subplots(ncols=2)

ax[0].plot(pulse.T_ps*1000,np.abs(pulse.AT)**2)
ax[0].plot(pulse_spm.T_ps*1000,np.abs(pulse_spm.AT)**2)
ax[0].plot(pulse_pynlo.T_ps*1000,np.abs(pulse_pynlo.AT)**2,'--')

axP = ax[0].twinx()
axP.plot(pulse.T_ps*1000,pulse.get_phase(domain='time'))
axP.plot(pulse.T_ps*1000,pulse_spm.get_phase(domain='time'))
axP.plot(pulse_pynlo.T_ps*1000,pulse_pynlo.get_phase(domain='time'),'--')


ax[1].plot(pulse.F_THz,np.abs(pulse.AW)**2)
ax[1].plot(pulse_spm.F_THz,np.abs(pulse_spm.AW)**2)
ax[1].plot(pulse_pynlo.F_THz,np.abs(pulse_pynlo.AW)**2,'--')

axP = ax[1].twinx()
axP.plot(pulse.F_THz,pulse.get_phase(domain='freq'))
axP.plot(pulse_spm.F_THz,pulse_spm.get_phase(domain='freq'))
axP.plot(pulse_pynlo.F_THz,pulse_pynlo.get_phase(domain='freq'),'--')

# %%
# data_t = np.abs(np.fft.fftshift(np.fft.fft(np.sqrt(data[:,1]))))**2
# data_time = np.fft.fftshift(np.fft.fftfreq(len(data[:,0]),data[1,0]-data[0,0]))

data_t = np.abs(np.fft.ifft(np.fft.ifftshift(np.sqrt(data[:,1])*np.exp(-1j*data[:,2])*np.exp(1j*data[:,0]*1.025*2*np.pi))))**2
data_time = np.fft.fftshift(np.fft.fftfreq(len(data[:,0]),data[1,0]-data[0,0]))
print(np.max(np.abs(pulse.AT))**2,np.max(data_t),np.max(np.abs(pulse.AT))**2/np.max(data_t))


print(np.max(np.abs(pulse.AW)**2),data[:,1].max())
plt.figure()
plt.plot(data_time,data_t/np.max(data_t))
plt.plot(pulse.T_ps, np.abs(pulse.AT)**2/np.max(np.abs(pulse.AT)**2),'--')


print('Peak Power: ', np.max(np.abs(pulse.AT)**2)/10**10)
# %%

plt.figure()
plt.plot(data[:,0],data[:,1]/np.max(data[:,1]))
plt.plot(data[:,0],data[:,2])
# %%

pulse = spym.createPulseSpecNLO(data[:,0], data[:,1],data[:,2],EPP=0.001)

print(np.max(np.abs(pulse.AW))**2,np.max(data[:,1]))


fig, ax = plt.subplots()
ax.plot(pulse.F_THz,np.abs(pulse.AW)**2/np.max(np.abs(pulse.AW)**2))
axP = ax.twinx()
axP.plot(pulse.F_THz,pulse.get_phase())

print(pulse.fwhm_t())
pulse.add_phase([-1000,-10000])

print(pulse.fwhm_t(),441/15.0)
