import numpy as np
import matplotlib.pyplot as plt
import datetime
from numpy import nan
import copy
import random
import pandas as pd
from scipy.interpolate import make_interp_spline
import csv
import os


# THIS CODE TAKES AN ARRAY OF DURATIONS AND AN ARRAY OF PARAMETER COMBINATIONS AND GENERATES A DATASET FOR EACH PARAMETER COMBINATION AND TIME DURATION
# IT IS NOT SEEDED, SO THE DATSETS WILL NOT BE CONSISTENT ACROSS DURATIONS

# For each duration and parameter combination, it generates two datasets: one with the suffix "noGroups" and one with the suffix "groups"
# These 2 datasets are the same, except the "groups" dataset has 3 more neurons (currently neurons 51, 52, 53).
# These additional neurons represent a contrived assembly, with neuron 51 always firing 0.025 seconds after neuron 50, 52 always firing 0.025 seconds after 51, etc.



# CHANGEABLE PARAMETERS--------------------------------------------------------------------------------------------------------------------------------------

#Durations (seconds)
durations = [10, 40, 80, 120, 160, 320, 640]



# Define parameter combinations
# r: firing rate (hz) of the poission spike train input to the simulation
# using_sine_wave influences stimVector, which defines the oscilating current that modulates the activation/firing of the simulated neurons. Mimics the pattern of an infraslow current.
# it is a sine function. 1 cycle every 4 seconds. Oscilating current.
# if the sine wave is desired, set using_sine_wave to true, else false

parameter_combinations = [
    {"r": 30 / 1000, "using_sine_wave": True},
    {"r": 60 / 1000, "using_sine_wave": True},
    {"r": 30 / 1000, "using_sine_wave": False},
    {"r": 60 / 1000, "using_sine_wave": False}
]

#the number of neurons in the simulation.
neuron_num = 50

#IF DESIRED, YOU CAN ALSO CHANGE THE FILE NAMES AND FOLDER ON LINES 79/80

# CHANGEABLE PARAMETERS--------------------------------------------------------------------------------------------------------------------------------------

for duration in durations:
    for params in parameter_combinations:
        # Neuron model parameters
        Vthresh = -50  # mV, AP (action potential) threshold voltage
        Vreset = -70  # mV, resting potential (value after spike)
        Vspike = -50  # mV, spike potential
        Rm = 10  # MOhms, membrane resistance
        tau = 10  # ms, membrane time constant (Rm * Cm, where Cm is the membrane capacitance)
        dt = 1  # ms, time step size
        counter = 0
        i = 0  # counting the number of iterations
        Se = 0.01  # synaptic conductance increment after a spike
        counter_poisson = 0  # counting number of spikes
        tau_e = 5  # ms, excitatory synaptic conductance time constant
        Ve = 0  # mV, excitatory reversal potential

        # The sine wave current
        I_max = 1.8  # nA, maximum external current
        I_min = 1.7  # nA, minimum external current
        I_mid = ((I_max - I_min) / 2) + I_min

        period = 4000  # ms, period of the sine wave current

        # Convert duration to milliseconds
        time = 1000 * duration

        r = params["r"]
        using_sine_wave = params["using_sine_wave"]

        # File naming based on parameters
        file_name = f"neuron_time_matrix_{duration}secs_r={int(r * 1000)}_{'yes' if using_sine_wave else 'no'}Sine"
        directory_path = os.path.join(r"../SimulatedDataSets/singleDatasetsUnseeded", f"{duration}SecondData")
        groups_filename = file_name + "_groups.csv"
        noGroups_filename = file_name + "_noGroups.csv"

        folder_path = os.path.join(directory_path, file_name)
        os.makedirs(folder_path, exist_ok=True)

        timeVector = np.arange(0, time, dt) #creating time vector of intervals of size dt.
        neuron_arr = np.zeros((neuron_num, len(timeVector))).astype(int) #2D array
        neuron_volt = np.zeros((neuron_num, len(timeVector))) # Creates a placeholder for our voltages that is the same size as timeVector

        # Creates a placeholder for the external stimulation vector.
        # It is also the same size as the time vector.
        stimVector = np.zeros(len(timeVector))

        # Generate the stimulation vector based on sine wave or constant input
        if using_sine_wave:
            stimVector[0:] = ((I_max - I_min) / 2) * (np.sin((timeVector / (period / (2 * np.pi))))) + I_mid
        else:
            stimVector[0:] = I_mid

        voltageVector_neurons = np.zeros((neuron_num, len(timeVector))).astype(int) #vector containing the voltages for each neuron at each time bin
        voltageVector_neurons_avg = np.zeros(len(timeVector)) #vector containing the average voltage across neurons

        for n in range(neuron_num):
            print(params, ",", n)
            counter_poisson = 0
            counter = 0
            voltageVector = np.zeros(len(timeVector))
            voltageVector[0] = Vreset
            geVector = np.zeros(len(timeVector))
            spikeVector = np.zeros(len(timeVector)).astype(int)

            for S in range(len(timeVector) - 1):
                Vinf = ((Rm * stimVector[S]) / tau) + (Vreset - voltageVector[S]) / tau + (Rm * geVector[S] * (Ve - voltageVector[S])) / tau
                voltageVector[S + 1] = Vinf * dt + voltageVector[S]
                geVector[S + 1] = geVector[S] * np.exp(-dt / tau_e)

                if voltageVector[S + 1] >= Vthresh:
                    voltageVector[S + 1] = Vspike

                    if voltageVector[S] == Vspike:
                        spikeVector[S] = 1
                        counter += 1
                        voltageVector[S + 1] = Vreset

                curTime = timeVector[S]
                while curTime < timeVector[S + 1]:
                    curTime = -np.log(np.random.uniform(0, 1)) / r + timeVector[S]

                    if curTime <= timeVector[S + 1]:
                        geVector[S + 1] = geVector[S] + Se

            neuron_arr[n] = spikeVector
            neuron_volt[n] = voltageVector
            voltageVector_neurons[n, ] = voltageVector

        for i in range(len(voltageVector)):
            voltageVector_neurons_avg[i, ] = voltageVector_neurons[:, i].mean()

        spike_time_matrix = []

        for neuron in range(neuron_num):
            spikeVector = neuron_arr[neuron]
            spike_time_Vector = np.zeros(sum(spikeVector), dtype=int)
            z = 0
            for y in range(len(spikeVector)):
                if spikeVector[y] == 1:
                    spike_time_Vector[z] = timeVector[y]
                    z = z + 1
            spike_time_matrix.append(list(spike_time_Vector))

        spike_time_matrix_seconds = [[value / 1000.0 for value in sublist] for sublist in spike_time_matrix]

        max_row_length = max(len(row) for row in spike_time_matrix_seconds)

        for row in spike_time_matrix_seconds:
            row.extend([nan] * (max_row_length - len(row)))

        spike_time_matrix_padded = np.array(spike_time_matrix_seconds)

        noGroups_filePath = os.path.join(folder_path, noGroups_filename)
        with open(noGroups_filePath, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(spike_time_matrix_padded)

        spike_time_matrix_seconds_groups = copy.deepcopy(spike_time_matrix_seconds)

        last_row = spike_time_matrix_seconds_groups[-1]
        for i in range(1, 4):
            new_row = [round(value + (0.025 * i), 4) if not np.isnan(value) else np.nan for value in last_row]
            spike_time_matrix_seconds_groups.append(new_row)

        spike_time_matrix_groups_padded = np.array(spike_time_matrix_seconds_groups)

        groups_filePath = os.path.join(folder_path, groups_filename)
        with open(groups_filePath, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(spike_time_matrix_groups_padded)