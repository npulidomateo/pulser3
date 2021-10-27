    
    # Pulser3 - Create 2-qubit sequences for cycle benchmarking.
    # Copyright (C) 2021  The Ninety-niners.

    # This program is free software: you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation, either version 3 of the License, or
    # (at your option) any later version.

    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.

    # You should have received a copy of the GNU General Public License
    # along with this program.  If not, see <https://www.gnu.org/licenses/>.


import qiskit as qk
from qiskit.circuit.library.standard_gates import RGate, IGate, PhaseGate, RXXGate
from qiskit.quantum_info import partial_trace
from qiskit.visualization import plot_bloch_multivector, plot_histogram
import numpy as np
import matplotlib.pyplot as plt
from qiskit.visualization.state_visualization import plot_bloch_multivector

# Alias numpy pi for convenience
pi = np.pi

# Define custom gates
# https://quantumcomputing.stackexchange.com/a/13160
# Use via e.g. `qc.append(ms_gate, [0, 1])`

ms_qc = qk.QuantumCircuit(2, name='MS')
ms_qc.h(0)
ms_qc.cx(0, 1)
ms_gate = ms_qc.to_instruction()

B_x_qc = qk.QuantumCircuit(1, name='B x')
B_x_qc.append(RGate(theta=pi/2, phi=0), [0])
B_x = B_x_qc.to_instruction()

B_y_qc = qk.QuantumCircuit(1, name='B y')
B_y_qc.append(RGate(theta=pi/2, phi=pi/2), [0])
B_y = B_y_qc.to_instruction()

B_z_qc = qk.QuantumCircuit(1, name='B z')
B_z_qc.append(PhaseGate(theta=pi/2), [0])
B_z = B_z_qc.to_instruction()

R_x_qc = qk.QuantumCircuit(1, name='R x')
R_x_qc.append(RGate(theta=pi, phi=0), [0])
R_x = R_x_qc.to_instruction()

R_y_qc = qk.QuantumCircuit(1, name='R y')
R_y_qc.append(RGate(theta=pi, phi=pi/2), [0])
R_y = R_y_qc.to_instruction()

R_z_qc = qk.QuantumCircuit(1, name='R z')
R_z_qc.append(PhaseGate(theta=pi), [0])
R_z = R_z_qc.to_instruction()

R_mz_qc = qk.QuantumCircuit(1, name='R -z')
R_mz_qc.append(PhaseGate(theta=-pi), [0])
R_mz = R_mz_qc.to_instruction()

R_mx_qc = qk.QuantumCircuit(1, name='R -x')
R_mx_qc.append(RGate(theta=pi, phi=0+pi), [0])
R_mx = R_mx_qc.to_instruction()

R_my_qc = qk.QuantumCircuit(1, name='R -y')
R_my_qc.append(RGate(theta=pi, phi=pi/2+pi), [0])
R_my = R_my_qc.to_instruction()

I_qc = qk.QuantumCircuit(1, name='I')
I_qc.append(IGate(), [0])
I = I_qc.to_instruction()

# Good ol' Pauli matrices
sx = np.array([[0, 1],[1, 0]])
sy = np.array([[0, -1j],[1j, 0]])
sz = np.array([[1, 0],[0, -1]])

# https://quantumcomputing.stackexchange.com/a/16921
# Style definition for qc.draw('mpl', style=style)
style={'displaycolor': {
        'B x': ('#00ff00', '#000000'),
        'B y': ('#00ff00', '#000000'),
        'B z': ('#00ff00', '#000000'),
        'MS': ('#ff5599', '#000000'),
        }}

# Make the generated sequences reproducible


class Pulser:
    """Create, store and compile quantum circuits"""

    svsim = qk.Aer.get_backend('statevector_simulator')

    def __init__(self, seed=999):
        self.circuits = []
        self.legends_l_k_m = []
        self.final_states = []
        
        self.rng = np.random.default_rng(seed)

        

    @staticmethod
    def cartesian_to_spherical(vector_car):
    # https://en.wikipedia.org/wiki/Spherical_coordinate_system#Cartesian_coordinates
        r = np.linalg.norm(vector_car)
        x, y, z = vector_car
        theta = np.arccos(z / r)
        if x > 0:
            phi = np.arctan(y/ x)
        elif x < 0:
            phi = np.arctan(y/x) + np.pi
        else:
            phi = 0
        return np.array([r, theta, phi]).transpose()


    def get_statevector(self, qc):
        return self.svsim.run(qk.transpile(qc, self.svsim)).result().get_statevector()

    def get_counts(self, qc):
        return self.svsim.run(qk.transpile(qc, self.svsim)).result().get_counts()

    def get_subsystems(self, statevector):
        return partial_trace(statevector, [1]), partial_trace(statevector, [0])

    def back_to_z(self, qc):
        """ Go back to north or south poles.

        params: qc
        returns: list of engineered_gates (two)
        """
        final_state = self.rng.choice([0, 1])
        if not final_state:
            dest_vector = np.array([0, 0, -1])
        else:
            dest_vector = np.array([0, 0, 1])
        
        engineered_gates = [None for _ in range(2)]
        subsystems = self.get_subsystems(self.get_statevector(qc))
        sigmas = [sx, sy, sz]
        
        for qubit, subsystem in enumerate(subsystems):
            bloch_vector = np.ndarray(shape=(3,), dtype=float)

            for axis, sigma in enumerate(sigmas):
                proj = np.real(np.trace(np.matmul(subsystem, sigma)))
                bloch_vector[axis] = proj
                
            rotation_axis = np.cross(bloch_vector, dest_vector)
            if np.linalg.norm(rotation_axis) <= 1e-5 and np.inner(bloch_vector, dest_vector) < 0:
                phi = 0
                rotation_angle = pi
                print()
                print('bloch_vector', bloch_vector)
                print('dest_vector', dest_vector)
                print('rotation_axis', rotation_axis)
                print('anti-parallel --> pi flip')
            else:
                print()
                print('bloch_vector', bloch_vector)
                print('dest_vector', dest_vector)
                print('rotation_axis', rotation_axis)
                rotation_angle = np.arcsin(np.linalg.norm(rotation_axis) / np.linalg.norm(bloch_vector) / np.linalg.norm(dest_vector))
                r, theta, phi = self.cartesian_to_spherical(rotation_axis)
                print('rotation_angle', rotation_angle)
            if np.abs(rotation_angle) < 1e-5:
                gate = IGate()
            else:
                gate = RGate(rotation_angle, phi)
            engineered_gates[qubit] = gate
        print(self.get_counts(qc))
        qc2 = qc.copy()
        for i, gate in enumerate(engineered_gates):
            qc2.append(gate, [i])
        print(self.get_counts(qc2))
        return engineered_gates


    def cycle_benchmark(self, m_list=[4, 8], L=5, final_pulses=True):
        print(__name__)

        # Prepare basis change combinations
        B_operators = [B_x, B_y, B_z, I]
        B_combinations = []
        for Bi in B_operators:
            for Bj in B_operators:
                B_combinations.append((Bi, Bj))
        # B_combinations = B_combinations[:-1]

        # Paulis to choose from:
        paulis = [R_x, R_mx, R_y, R_my, R_z, R_mz, I]

        # Assemble and store circuits
        for l in range(L):
            for k, (Bi, Bj) in enumerate(B_combinations):
                for m in m_list:
                    qc = qk.QuantumCircuit(2)
                    qc.append(Bi, [0])
                    qc.append(Bj, [1])
                    qc.append(self.rng.choice(paulis), [0])
                    qc.append(self.rng.choice(paulis), [1])
                    for _ in range(m):
                        qc.append(RXXGate(pi/2), [0, 1])
                        qc.append(self.rng.choice(paulis), [0])
                        qc.append(self.rng.choice(paulis), [1])
                    # calculate and append the engineered gates
                    if final_pulses:
                        engineered_gates = self.back_to_z(qc)
                        qc.append(engineered_gates[0], [0])
                        qc.append(engineered_gates[1], [1])
                    self.circuits.append(qc)
                    self.legends_l_k_m.append((l, k, m))


def main_cycle_benchmark():
    pulser = Pulser()
    pulser.cycle_benchmark(L=1, m_list=[8])
     # Print circuits
    for i, (legend, qc) in enumerate(zip(pulser.legends_l_k_m, pulser.circuits)):
        print(i, '(l, k, m):', legend)
        # print(qc)
        print(pulser.get_counts(qc))


def main_debug_back_to_z():
    pulser = Pulser()
    pulser.cycle_benchmark(L=1, m_list=[2], final_pulses=False)
    for i, qc in enumerate(pulser.circuits):
        if not(i == 15 or True):
            continue
        print('\n\ncircuit', i)
        engineered_gates = pulser.back_to_z(qc)
        qc.append(engineered_gates[0], [0])
        qc.append(engineered_gates[1], [1])
        print('circuits from main')
        print(qc)
        counts = pulser.get_counts(qc)
        print(qc.qasm())
        print('counts from main', counts)
        # plot_bloch_multivector(qc)
        # plt.show()

def main_nonconsistent_epulses():
    my_pulser = Pulser()
    your_pulser = Pulser()

    my_pulser.cycle_benchmark(L=1, m_list=[4], final_pulses=False)
    your_pulser.cycle_benchmark(L=1, m_list=[4], final_pulses=True)
    epulses = your_pulser.back_to_z(your_pulser.circuits[0])

    qc = qk.QuantumCircuit(2)
    qc.append(epulses[0],[0])
    qc.append(epulses[1],[1])
    print(qc)

    print(my_pulser.circuits[0])
    print(your_pulser.circuits[0])


if __name__ == '__main__':
    main_debug_back_to_z()

