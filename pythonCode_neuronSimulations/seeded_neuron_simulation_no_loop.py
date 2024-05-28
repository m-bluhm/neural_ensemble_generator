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


# THIS CODE TAKES AN ARRAY OF DURATIONS AND AN ARRAY OF PARAMETER COMBINATIONS AND GENERATES A SEEDED DATASET FOR EACH PARAMETER COMBINATION AND TIME DURATION

#it is seeded in these sense that the datasets for different durations will be the same until the shorter dataset ends.
#For example, the first 40 seconds of a 160 second dataset will be identical to the 40 second dataset

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

# Initial seed number to generate consistent datasets.
seed_num = 42

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
        file_name = f"seeded_neuron_time_matrix_{duration}secs_r={int(r * 1000)}_{'yes' if using_sine_wave else 'no'}Sine"
        directory_path = os.path.join(r"../SimulatedDataSets/singleDatasetsSeeded", f"seeded{duration}SecondData")
        groups_filename = file_name + "_groups.csv"
        noGroups_filename = file_name + "_noGroups.csv"

        folder_path = os.path.join(directory_path)
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

        voltageVector_neurons = np.zeros((neuron_num, len(timeVector))).astype(int)
        voltageVector_neurons_avg = np.zeros(len(timeVector))

        #Loop through each neuron
        for n in range(neuron_num):
            # Seed the random number generator for reproducibility of results
            np.random.seed(seed_num + n)
            print(params, ",", n)  # Print current parameters and neuron index for tracking progress

            # Initialize counters and vectors for simulation
            counter_poisson = 0  # Counter for Poisson spikes (not explicitly used in this block)
            counter = 0  # Counter for action potentials
            voltageVector = np.zeros(len(timeVector))  # Vector to store membrane potential over time
            voltageVector[0] = Vreset  # Set initial membrane potential to resting potential
            geVector = np.zeros(len(timeVector))  # Vector to store excitatory synaptic conductance over time. # amount spike jumps, excitatory conductance
            spikeVector = np.zeros(len(timeVector)).astype(int)  # Vector to record spikes (1 if spike occurs, 0 otherwise)

            # Simulate neuron dynamics over time
            # This line initiates the loop. "S" counts the number of loops.
            # We are looping for 1 less than the length of the time vector
            # because we have already calculated the voltage for the first iteration.
            for S in range(len(timeVector) - 1):
                # Calculate the steady-state membrane potential (Vinf) for the next time step
                Vinf = ((Rm * stimVector[S]) / tau) + (Vreset - voltageVector[S]) / tau + (
                            Rm * geVector[S] * (Ve - voltageVector[S])) / tau
                # ((E + I*R)/tau) + (-v[i]/tau) #differential equation

                # Vinf is the value of the membrane potential at the next time step. It is the sum of three components:

                # The first term represents the influence of the external stimulus (Istim) on the membrane potential. It's essentially the effect of the applied current over time.

                # The second term represents the influence of the difference between the reset potential (Vreset) and the current membrane potential (Vprevious).
                # This term introduces a form of negative feedback, pushing the membrane potential back towards the reset potential.

                # The third term represents the influence of the excitatory conductance (ge) on the membrane potential.
                # It accounts for the effect of excitatory synaptic inputs. Ve is the excitatory reversal potential.

                # Tau: Membrane time constant, representing the time it takes for the membrane potential to change in response to applied current.
                # Ve: Excitatory reversal potential, the equilibrium potential for excitatory inputs.

                # Update the membrane potential for the next time step using Euler's method
                voltageVector[S + 1] = Vinf * dt + voltageVector[S]

                # Update the synaptic conductance, decaying exponentially
                # upadate conductance decay as if no poisson input spike
                # conductance of external input into the neuron
                geVector[S + 1] = geVector[S] * np.exp(-dt / tau_e)

                # Check if the membrane potential exceeds the action potential threshold
                if voltageVector[S + 1] >= Vthresh:
                    voltageVector[S + 1] = Vspike  # Set membrane potential to spike value

                    # This 'if' statement checks if we are already at Vspike (this is another way we can be above Vthresh)
                    if voltageVector[S] == Vspike: # if the current voltage is a spike
                        spikeVector[S] = 1  # Record the spike occurrence
                        counter += 1  # Increment spike counter
                        voltageVector[S + 1] = Vreset  # Reset membrane potential to resting potential

                # Implement the Poisson process for generating synaptic input spikes
                curTime = timeVector[S]
                while curTime < timeVector[S + 1]:
                    # Generate the next spike time using the inverse transform method for the exponential distribution
                    curTime = -np.log(np.random.uniform(0, 1)) / r + timeVector[S]

                    # r in ms Calculate wait time. PROPOSED SPIKE TIME
                    # generates a random number from an exponential distribution with a mean of 1/r. This random number represents the time until the next spike
                    # If the generated curTime is within the current time interval (timeVector[S] to timeVector[S + 1]),
                    # it updates the excitatory conductance geVector by adding the strength of a spike Se.
                    # This simulates the effect of an excitatory spike arriving at the neuron during that time interval

                    # If the generated spike time is within the current time step, update the conductance
                    if curTime <= timeVector[S + 1]:
                        geVector[S + 1] = geVector[S] + Se  # Increase synaptic conductance by the strength of the spike

            # Store the spike train and voltage trace for the current neuron
            # INEFFICIENT
            neuron_arr[n] = spikeVector  # Record the spike train for the current neuron
            neuron_volt[n] = voltageVector  # Record the voltage trace for the current neuron
            voltageVector_neurons[n,] = voltageVector  # Store the voltage trace in the matrix for averaging later

        for i in range(len(voltageVector)):
            voltageVector_neurons_avg[i, ] = voltageVector_neurons[:, i].mean()

        spike_time_matrix = []
        # Vector of specific times of spikes

        for neuron in range(neuron_num):
            spikeVector = neuron_arr[neuron]
            spike_time_Vector = np.zeros(sum(spikeVector), dtype=int)
            z = 0
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

        # Convert spike_time_matrix_seconds to a NumPy array
        spike_time_matrix_padded = np.array(spike_time_matrix_seconds)

        # Write "noGroups" dataset to CSV
        noGroups_filePath = os.path.join(folder_path, noGroups_filename)
        with open(noGroups_filePath, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(spike_time_matrix_padded)

        # this code block copies the simulated dataset to make and write the "groups" dataset.
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