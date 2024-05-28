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


# THIS CODE TAKES AN ARRAY OF DURATIONS AND AN ARRAY OF PARAMETER COMBINATIONS AND GENERATES MULTIPLE DATASETS FOR EACH PARAMETER COMBINATION AND TIME DURATION

# For each duration and parameter combination, it generates two datasets: one with the suffix "noGroups" and one with the suffix "groups"
# These 2 datasets are the same, except the "groups" dataset has 3 more neurons (currently neurons 51, 52, 53).
# These additional neurons represent a contrived assembly, with neuron 51 always firing 0.025 seconds after neuron 50, 52 always firing 0.025 seconds after 51, etc.



# CHANGEABLE PARAMETERS--------------------------------------------------------------------------------------------------------------------------------------

#Durations (seconds)
durations = [10, 40, 80, 120, 160, 320, 640]


# Define parameter combinations
# r: firing rate (hz) of the poission spike train input to the simulation
# using_sine_wave: Determines if the simulation has the modulating sine wave input representing the infraslow oscilation. If it is set to false, then the sine input is just a flat line
# these are the only parameters we talked about looking at.

parameter_combinations = [
    {"r": 30 / 1000, "using_sine_wave": True},
    {"r": 60 / 1000, "using_sine_wave": True},
    {"r": 30 / 1000, "using_sine_wave": False},
    {"r": 60 / 1000, "using_sine_wave": False}
]

# the number of neurons in the simulation.
neuron_num = 50

# how many times each dataset is generated. The seed changes each time. This is useful for determining if certain seeds have unpredicted results.
num_simulations = 10

#IF DESIRED, YOU CAN ALSO CHANGE THE FILE NAMES AND FOLDER ON LINES 79/80

# CHANGEABLE PARAMETERS--------------------------------------------------------------------------------------------------------------------------------------

for duration in durations:
    for param_index, params in enumerate(parameter_combinations):
        name_prefix = 0
        for simulation_num in range(1, num_simulations + 1):
            Vthresh = -50   # mV, AP threshold voltage
            Vreset = -70   # mV resting potential (going to value after spike)
            Vspike = -50   # mV
            Rm = 10   # MOhms lower Resistance = lower spikes
            tau = 10   # ms
            dt = 1  # ms, time bins. SHOULD GO DOWN????
            counter = 0
            i = 0  # counting the number of iterations
            Se = 0.01  # strength of spike
            counter_poisson = 0  # counting number of spikes
            tau_e = 5  # excitatory conductance time constant
            Ve = 0  # excitatory reversal potential

            # the sine wave current
            I_max = 1.8  # external current
            I_min = 1.7  # external current
            I_mid = ((I_max - I_min) / 2) + I_min


            period = 4000  # in ms = 4 seconds

            time = 1000 * duration  # Convert duration to milliseconds

            r = params["r"]
            using_sine_wave = params["using_sine_wave"]

            param_identifier = f"r={int(r * 1000)}_{'yes' if using_sine_wave else 'no'}Sine"
            file_name = f"sim_{name_prefix}_neuron_time_matrix_{duration}secs_r={int(r * 1000)}_{'yes' if using_sine_wave else 'no'}Sine"
            directory_path = os.path.join(r"../SimulatedDataSets/multipleDatasetsUnseeded", f"unseeded{duration}SecondData", param_identifier)
            groups_filename = file_name + "_groups.csv"
            noGroups_filename = file_name + "_noGroups.csv"
            name_prefix += 1

            folder_path = os.path.join(directory_path)
            os.makedirs(folder_path, exist_ok=True)

            timeVector = np.arange(0, time, dt)
            neuron_arr = np.zeros((neuron_num, len(timeVector))).astype(int)
            neuron_volt = np.zeros((neuron_num, len(timeVector)))
            stimVector = np.zeros(len(timeVector))

            if using_sine_wave:
                stimVector[0:] = ((I_max - I_min) / 2) * (np.sin((timeVector / (period / (2 * np.pi))))) + I_mid
            else:
                stimVector[0:] = I_mid

            voltageVector_neurons = np.zeros((neuron_num, len(timeVector))).astype(int)
            voltageVector_neurons_avg = np.zeros(len(timeVector))

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
