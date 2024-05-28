import numpy as np
import matplotlib.pyplot as plt
import datetime
from numpy import nan
import copy
import random
import pandas as pd
from scipy.interpolate import make_interp_spline ##new import
import csv
import os
# Set simulation parameters

# for sim_num in range (10):
Vthresh = -50   #mV, AP threshold voltage
Vreset  = -70   #mV resting potential (going to value after spike)
Vspike  =  -50   #mV
Rm      =  10   #MOhms lower Resistance = lower spikes
tau     =  10   #ms
dt      =  1 #ms, time bins. SHOULD GO DOWN????
counter =  0
i = 0 #counting the number of iterations
Se = 0.01 #strength of spike
counter_poisson = 0 #counting number of spikes
tau_e = 5 #excitatory conductance time constant
Ve = 0 #excitatory reversal potential

#Sine funciton variables
I_max =  1.8 #external current
I_min =  1.7 #external current
I_mid = ((I_max - I_min)/2) + I_min
period = 4000 #in ms = 4 seconds

# ----------------------------------------------------------------------------------------------------------------------#
# CHANGEABLE PARAMETERS (CURRENTLY)
time = 1000*160         # In ms, so do 1000 * desired number of seconds. This is the number of seconds that the simulation will run

r = 30/1000             # The frequency (in #spikes/ms) of the poisson input. 30/1000 means 30hz
                             #FOR JEN: IS THIS RIGHT?

np.random.seed(42)      # Call np.random.seed(int) before any other np calls in the function.
                        # This will ensure that anytime the code is ran, you will get the same sequence of random numbers for the poission input.
                        # The specific integer doesn't matter, 42 is the answer to everything, so we're using it.

using_sine_wave = True
# stimVector defines the oscilating current that modulates the activation/firing of the simulated neurons. Mimics the pattern of an infraslow current.
# it is a sine function. 1 cycle every 4 seconds. Oscilating current.
# if the sine wave is desired, seet using_sine_wave to false, else true

neuron_num = 50 #number of neurons being simulated

file_name = "seeded_neuron_time_matrix_160secs_r=30_yesSine"
directory_path = os.path.join("..", "seeded160SecondData")
groups_filename = file_name+"_groups.csv"
noGroups_filename = file_name+"_noGroups.csv"

folder_path = os.path.join(directory_path, file_name)
os.makedirs(folder_path, exist_ok=True)

#----------------------------------------------------------------------------------------------------------------------#


timeVector = np.arange(0, time, dt) #creating time vector of intervals of size dt.
neuron_arr = np.zeros((neuron_num, len(timeVector))).astype(int) #2D array
neuron_volt = np.zeros((neuron_num, len(timeVector))) # Creates a placeholder for our voltages that is the same size as timeVector

# Creates a placeholder for the external stimulation vector.
# It is also the same size as the time vector.
stimVector = np.zeros(len(timeVector))

if using_sine_wave:
    stimVector[0:] = ((I_max - I_min)/2) * (np.sin((timeVector/(period/(2*np.pi))))) + I_mid #if sine wave
else:
    stimVector[0:] = I_mid



# Euler---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

voltageVector_neurons = np.zeros((neuron_num, len(timeVector))).astype(int) #vector containing the voltages for each neuron at each time bin
voltageVector_neurons_avg = np.zeros(len(timeVector)) #vector containing the average voltage across neurons
#WAYS TO IMPROVE:---------------------------------------------------------------------------------------------------------------------------------------------
#run neurons in parallel (only when not talking to eachother)
#compute external stimulus as you go. Don;t need to hold onto whole vector
#MAYBE: Hold onto time that the neuron spiked.


for n in range(neuron_num): #for each neuron
    print(n)
    counter_poisson = 0 #number of poisson spikes equals 0
    counter = 0
    voltageVector = np.zeros(len(timeVector)) #create a new voltage vector filled with all 0s
    voltageVector[0] = Vreset  # Set the initial voltage to be equal to the resting potential
    geVector = np.zeros(len(timeVector))  # amount spike jumps, excitatory conductance
    spikeVector = np.zeros(len(timeVector)).astype(int)  # holding 1s for spike and 0 for no spike at given interval
    startTime = datetime.datetime.now()

    # This line initiates the loop. "S" counts the number of loops.
    # We are looping for 1 less than the length of the time vector
    # because we have already calculated the voltage for the first
    # iteration.
    for S in range(len(timeVector) - 1):

        # update Vinf using purple sheet. WHAT IS PURPLE SHEET? I HAVE FORGOTTEN
        Vinf = ((Rm * stimVector[S]) / tau) + (Vreset - voltageVector[S]) / tau + (Rm * geVector[S] * (
                    Ve - voltageVector[S])) / tau  # ((E + I*R)/tau) + (-v[i]/tau) #differential equation
        #MAX---------------------------------------------------------------------------------------------------
        #relabel vinf to what?

        #Vinf is the value of the membrane potential at the next time step. It is the sum of three components:

        #The first term represents the influence of the external stimulus (Istim) on the membrane potential. It's essentially the effect of the applied current over time.

        #The second term represents the influence of the difference between the reset potential (Vreset) and the current membrane potential (Vprevious).
        #This term introduces a form of negative feedback, pushing the membrane potential back towards the reset potential.

        #The third term represents the influence of the excitatory conductance (ge) on the membrane potential.
        # It accounts for the effect of excitatory synaptic inputs. Ve is the excitatory reversal potential.

        #Tau: Membrane time constant, representing the time it takes for the membrane potential to change in response to applied current.
        #Ve: Excitatory reversal potential, the equilibrium potential for excitatory inputs.

        voltageVector[S + 1] = Vinf * dt + voltageVector[S]  # v[i+1] = v[i] + dt *Vinf
        #Euler step-- equation of a line

        # upadate conductance decay as if no poisson input spike
        #conductance of external input into the neuron
        #MAX: WHY DO WE DO THIS?
        geVector[S + 1] = geVector[S] * np.exp(-dt / tau_e)

        # This 'if' condition states that if the next voltage is greater than
        # or equal to the threshold, then to run the next section
        if voltageVector[S + 1] >= Vthresh:
            # This states that the next voltage vector will be the Vspike value
            voltageVector[S + 1] = Vspike #neuron spike if above threshold

            # This 'if' statement checks if we are already at Vspike (this is
            # another way we can be above Vthresh)
            if voltageVector[S] == Vspike: #if the current voltage is a spike
                spikeVector[S] = 1  # add a 1 to indicate a spike in spikeVector
                # upadate conductance as if no spike
                #geVector[i] = geVector[i-1] * np.exp(-dt/curTime) #MAX: WHAT IS THE PURPOSE OF THIS?
                counter += 1

                # Set the next voltage equal to the reset value
                voltageVector[S + 1] = Vreset #MAYBE IMPLEMENT A "vRefractory--" less than Vreset

        curTime = timeVector[S]
        while curTime < timeVector[S + 1]: #MAX: I don't understand the purpose of this while loop?
            # responsible for simulating the arrival of excitatory spikes:
            #Generates a random waiting time according to a Poisson process

            curTime = -np.log(np.random.uniform(0, 1)) / r + timeVector[S]  # r in ms Calculate wait time. PROPOSED SPIKE TIME
            # generates a random number from an exponential distribution with a mean of 1/r. This random number represents the time until the next spike
            #If the generated curTime is within the current time interval (timeVector[S] to timeVector[S + 1]),
            # it updates the excitatory conductance geVector by adding the strength of a spike Se.
            # This simulates the effect of an excitatory spike arriving at the neuron during that time interval

            if curTime <= timeVector[S + 1]:  # see if wait time is in interval
                # if so, update conductance
                geVector[S + 1] = geVector[S] + Se
                # spike will occur
                # If the generated curTime is within the current time interval (timeVector[S] to timeVector[S + 1]),
                # it updates the excitatory conductance geVector by adding the strength of a spike Se.
                # This simulates the effect of an excitatory spike arriving at the neuron during that time interval
        #INEFFICENT
        neuron_arr[n] = spikeVector
        neuron_volt[n] = voltageVector
        voltageVector_neurons[n, ] = voltageVector

#END LOOP

for i in range(len(voltageVector)):
    voltageVector_neurons_avg[i, ] = voltageVector_neurons[:, i].mean()

spike_time_matrix = []
# Vector of specific times of spikes
for neuron in range(neuron_num):
    spikeVector = neuron_arr[neuron]
    spike_time_Vector = np.zeros(sum(spikeVector), dtype=int)
    z = 0  # initialize z to begin at index #MAX: Rename? What is z?
    for y in range(len(spikeVector)):  # grab index of spike and compare to timeVector

        if spikeVector[y] == 1:
            spike_time_Vector[z] = timeVector[y]  # add spike times in this vector
            z = z + 1
    spike_time_matrix.append(list(spike_time_Vector))

spike_time_matrix_seconds = [[value / 1000.0 for value in sublist] for sublist in spike_time_matrix]

# Find the length of the longest row
max_row_length = max(len(row) for row in spike_time_matrix_seconds)

# Pad shorter rows with NaN to make them the same length as the longest row
for row in spike_time_matrix_seconds:
    row.extend([nan] * (max_row_length - len(row)))

#Convert spike_time_matrix_seconds to a NumPy array
spike_time_matrix_padded = np.array(spike_time_matrix_seconds)

noGroups_filePath = os.path.join(folder_path, noGroups_filename)
with open(noGroups_filePath, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(spike_time_matrix_padded)


spike_time_matrix_seconds_groups = copy.deepcopy(spike_time_matrix_seconds) #make a deep copy so that a matrix with assemblies can be generated

last_row = spike_time_matrix_seconds_groups[-1]  # Get the last row
for i in range(1, 4):  # Add three rows
    new_row = [round(value + (0.025 * i), 4) if not np.isnan(value) else np.nan for value in last_row]
    spike_time_matrix_seconds_groups.append(new_row)

#Convert spike_time_matrix_seconds to a NumPy array
spike_time_matrix_groups_padded = np.array(spike_time_matrix_seconds_groups)

groups_filePath = os.path.join(folder_path, groups_filename)
with open(groups_filePath, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(spike_time_matrix_groups_padded)

