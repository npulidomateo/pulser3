    # Pulser3 - Create pulse sequences for cycle benchmarking
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


from qiskit import *
from qiskit.circuit.library.standard_gates import XGate, x, RGate, IGate, PhaseGate, YGate, ZGate
from qiskit.circuit import instruction
import numpy as np
import matplotlib.pyplot as plt

pi = np.pi

# Define custom gates
# https://quantumcomputing.stackexchange.com/a/13160
# Use via `qc.append(ms_gate, [0, 1])`
ms_qc = QuantumCircuit(2, name='MS')
ms_qc.h(0)
ms_qc.cx(0, 1)
ms_gate = ms_qc.to_instruction()

B_x_qc = QuantumCircuit(1, name='B x')
B_x_qc.append(RGate(theta=pi/2, phi=0), [0])
B_x = B_x_qc.to_instruction()

B_y_qc = QuantumCircuit(1, name='B y')
B_y_qc.append(RGate(theta=pi/2, phi=pi/2), [0])
B_y = B_y_qc.to_instruction()

B_z_qc = QuantumCircuit(1, name='B z')
B_z_qc.append(PhaseGate(theta=pi/2), [0])
B_z = B_z_qc.to_instruction()

R_x_qc = QuantumCircuit(1, name='R x')
R_x_qc.append(RGate(theta=pi, phi=0), [0])
R_x = R_x_qc.to_instruction()

R_y_qc = QuantumCircuit(1, name='R y')
R_y_qc.append(RGate(theta=pi, phi=pi/2), [0])
R_y = R_y_qc.to_instruction()

R_z_qc = QuantumCircuit(1, name='R z')
R_z_qc.append(PhaseGate(theta=pi), [0])
R_z = R_z_qc.to_instruction()

R_mx_qc = QuantumCircuit(1, name='R -x')
R_mx_qc.append(RGate(theta=pi, phi=0+pi), [0])
R_mx = R_mx_qc.to_instruction()

R_my_qc = QuantumCircuit(1, name='R -y')
R_my_qc.append(RGate(theta=pi, phi=pi/2+pi), [0])
R_my = R_my_qc.to_instruction()

I_qc = QuantumCircuit(1, name='I')
I_qc.append(IGate(), [0])
I = I_qc.to_instruction()

style={'displaycolor': {
        'B x': ('#00ff00', '#000000'),
        'B y': ('#00ff00', '#000000'),
        'B x': ('#00ff00', '#000000'),
        'MS': ('#ff5599', '#000000'),
        }}

# Make the generated sequences reproducible
np.random.seed(999)


class Pulser:
    def __init__(self):
        self.circuits = []

    def get_state(qc):
        pass

    def cycle_benchmark(self, m_list=[4, 8], L=5):

        # Prepare basis change combinations
        B_operators = [B_x, B_y, B_z, I]
        B_combinations = []
        for Bi in B_operators:
            for Bj in B_operators:
                B_combinations.append((Bi, Bj))

        # Paulis to choose from:
        paulis = [R_x, R_mx, R_y, R_my, R_z, I]

        # Assemble and store circuits
        for l in range(L):
            for k, (Bi, Bj) in enumerate(B_combinations):
                for m in m_list:
                    qc = QuantumCircuit(2)
                    qc.append(Bi, [0])
                    qc.append(Bj, [1])
                    qc.append(np.random.choice(paulis), [0])
                    qc.append(np.random.choice(paulis), [1])
                    for _ in range(m):
                        qc.append(ms_gate, [0, 1])
                        qc.append(np.random.choice(paulis), [0])
                        qc.append(np.random.choice(paulis), [1])
                    self.circuits.append(qc)
        
        # Print circuits
        for qc in self.circuits:
            print(qc)
            print()
                    

def main_custom_gates():
    qc = QuantumCircuit(2)
    qc.append(B_x, [0])
    qc.append(B_y, [1])
    qc.append(R_z, [0])
    qc.append(R_my, [1])
    qc.append(R_z, [1])
    qc.append(ms_gate, [0, 1])
    qc.append(I, [0])
    print(qc)
    # https://quantumcomputing.stackexchange.com/a/16921
    
    qc.draw('mpl', style=style)
    plt.show()


def main_cycle_benchmark():
    pulser = Pulser()
    pulser.cycle_benchmark()

if __name__ == '__main__':
    main_cycle_benchmark()
    


            

