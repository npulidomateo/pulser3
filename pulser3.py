    
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
from qiskit.circuit.library.standard_gates import RGate, IGate, RXGate, RYGate, RZGate, RXXGate
from qiskit.quantum_info import partial_trace
from qiskit.visualization import plot_bloch_multivector, plot_histogram
import numpy as np
from pathlib import Path
import os
import warnings
from multiprocessing import Pool

# Alias numpy pi for convenience
pi = np.pi

# Good ol' Pauli matrices
sx = np.array([[0, 1],[1, 0]])
sy = np.array([[0, -1j],[1j, 0]])
sz = np.array([[1, 0],[0, -1]])

# https://quantumcomputing.stackexchange.com/a/16921
# Style definition for qc.draw('mpl', style=style)
style={'displaycolor': {
        'Rx': ('#00ff00', '#000000'),
        'Ry': ('#00ff00', '#000000'),
        'Rz': ('#00ff00', '#000000'),
        'Rxx': ('#ff5599', '#000000'),
        }}


class Pulser:
    """Create, store, transpile and compile 2-qubit quantum circuits"""

    svsim = qk.Aer.get_backend('statevector_simulator')

    def __init__(self, seed=999):
        self.circuits = []
        self.legends_l_k_m = []
        self.final_states = []
        self.hfgui_sequences = []
        self.hfgui_circuits = []

        # Make the generated sequences reproducible
        self.rng = np.random.default_rng(seed)


    @staticmethod
    def cartesian_to_spherical(vector_car):
    # https://en.wikipedia.org/wiki/Spherical_coordinate_system#Cartesian_coordinates
        r = np.linalg.norm(vector_car)
        x, y, z = vector_car
        if r < 0.00001:
            return np.zeros(shape=(3,))
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

    @staticmethod
    def get_bloch_vector(rho):
        bloch_vector = np.ndarray(shape=(3,), dtype=float)
        sigmas = [sx, sy, sz]
        for axis, sigma in enumerate(sigmas):
            proj = np.real(np.trace(np.matmul(rho, sigma)))
            bloch_vector[axis] = proj
        return bloch_vector


    def back_to_z(self, qc, verbose=False):
        """ Go back to north or south poles.

        Pululates `self.final_states` with random choices 
        0 = |00>, 1 = |11>

        params: qc (two qubit quantum circuit)
        returns: list of engineered_gates (two)
        """
        final_state = self.rng.choice([0, 1])
        if not final_state:
            dest_vector = np.array([0, 0, -1])
        else:
            dest_vector = np.array([0, 0, 1])
        self.final_states.append(final_state)
        
        engineered_gates = [None for _ in range(2)]
        subsystems = self.get_subsystems(self.get_statevector(qc))

        for qubit, subsystem in enumerate(subsystems):
            bloch_vector = self.get_bloch_vector(subsystem)
            rotation_axis = np.cross(bloch_vector, dest_vector)
            if np.linalg.norm(rotation_axis) <= 1e-5 and np.inner(bloch_vector, dest_vector) < 0:
                phi = 0
                rotation_angle = pi
                if verbose:
                    print()
                    print('bloch_vector', bloch_vector)
                    print('dest_vector', dest_vector)
                    print('rotation_axis', rotation_axis)
                    print('anti-parallel --> pi flip')
            else:
                rotation_angle = np.arcsin(np.linalg.norm(rotation_axis) / np.linalg.norm(bloch_vector) / np.linalg.norm(dest_vector))
                r, theta, phi = self.cartesian_to_spherical(rotation_axis)
                if verbose:
                    print()
                    print('bloch_vector', bloch_vector)
                    print('dest_vector', dest_vector)
                    print('rotation_axis', rotation_axis)
                    print('rotation_angle', rotation_angle)
            if np.abs(rotation_angle) < 1e-5:
                gate = IGate()
            else:
                gate = RGate(rotation_angle, phi)
            engineered_gates[qubit] = gate
        
        if verbose:
            print(self.get_counts(qc))
            qc2 = qc.copy()
            for i, gate in enumerate(engineered_gates):
                qc2.append(gate, [i])
            print(self.get_counts(qc2))

        return engineered_gates


    def cycle_benchmark(self, m_list=[4, 8], L=5, final_pulses=True, do_ms_gate=True):
        self.circuits = []  # Reset circuits

        # Prepare basis change combinations
        B_operators = [RXGate(pi/2), RYGate(pi/2), RZGate(pi/2), IGate()]
        B_combinations = []
        for Bi in B_operators:
            for Bj in B_operators:
                B_combinations.append((Bi, Bj))
        # Don't use last combination (I, I)
        B_combinations = B_combinations[:-1]

        # Paulis to choose from:
        paulis = [RXGate(pi), RXGate(-pi), RYGate(pi), RYGate(-pi), RZGate(pi), RZGate(-pi), IGate()]

        # Assemble and store circuits and legends
        for l in range(L):
            for k, (Bi, Bj) in enumerate(B_combinations):
                for m in m_list:
                    qc = qk.QuantumCircuit(2)
                    qc.append(Bi, [0])
                    qc.append(Bj, [1])
                    i = self.rng.choice([idx for idx in range(len(paulis))])
                    j = self.rng.choice([idx for idx in range(len(paulis))])
                    Ri, Rj = paulis[i], paulis[j]
                    qc.append(Ri, [0])
                    qc.append(Rj, [1])
                    for _ in range(m):
                        i = self.rng.choice([idx for idx in range(len(paulis))])
                        j = self.rng.choice([idx for idx in range(len(paulis))])
                        Ri, Rj = paulis[i], paulis[j]
                        if do_ms_gate:
                            qc.append(RXXGate(pi/2), [0, 1])
                        else:
                            qc.append(IGate(), [0, 1])
                        qc.append(Ri, [0])
                        qc.append(Rj, [1])
                    # calculate and append the engineered gates
                    if final_pulses:
                        engineered_gates = self.back_to_z(qc)
                        qc.append(engineered_gates[0], [0])
                        qc.append(engineered_gates[1], [1])
                    self.circuits.append(qc)
                    self.legends_l_k_m.append((l, k, m))


    @staticmethod
    def SIA_pitime(word, ion):
        """Converts angles to SIA pitime for a given ion"""

        # get angle from expression
        exec('expression = ' + word, globals())
        angle = expression
        
        # Don't allow negative angles
        if angle < 0:
            angle = 2*pi + angle

        angle_in_pi_units = angle / pi
        new_word = '%f*scan_pi_%d' % (angle_in_pi_units, ion)
        return new_word


    @staticmethod
    def put_the_dot(word='pi_0 * 23.456/4 * (pi_1/pi_0) / 4. / *(42.0 + 34) * pi3'):
        """Add a dot to the end of floats in expressions to make hfgui smile"""
        
        separators = ' +-/*()'
        digits = '.0123456789'
        tainted = False
        new_word = ''
        number = ''
        
        # Parse the word
        for c in word:
            if c in separators and number != '':
                if '.' not in number:
                    number += '.'
                new_word += number
                new_word += c
                number = ''
            elif c in separators and number == '':
                tainted = False
                new_word += c
            elif c in digits and not tainted:
                number += c
            else:
                new_word += c
                tainted = True
        # Catch the last number
        if number != '':
            if '.' not in number:
                number += '.'
            new_word += number

        return new_word
            

    @staticmethod
    def get_from_parenthesis(text):
        sentence = []
        word = ''
        rec = False
        for character in text:
            if character == '(':
                rec = True
            elif character == ',':
                sentence.append(word)
                word = ''
            elif character == ')':
                if word == '':
                    sentence = ['', '']
                else:
                    sentence.append(word)
                break
            elif rec:
                word += character

        # remove the aczs compensation phase
        word = ''
        first_arg = sentence[0]
        for c in first_arg:
            if c in [' ', '+']:
                sentence[0] = word
            else:
                    word += c
        return sentence

    @staticmethod
    def adapt_phase_for_negative_rotations(theta, phi):
        if '-' in theta:
            print(theta)
            negative_theta = ''
            for c in theta:
                if c != '-':
                    negative_theta += c
            exec('phi_expr = %s' % phi, globals())
            adapted_phi = str((phi_expr + pi) % (2*pi)) 
            
            return negative_theta, adapted_phi
        
        else:
            return theta, phi

            

    def compile_circuit(self, qc, basis_gates=['rxx', 'r', 'id'], do_ms_gate=True, aczs_comp=True, transpile=True, transpile_optimization=0):
       
        # Transpile to something compatible with our system
        if transpile:
            # https://stackoverflow.com/a/17654868
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=UserWarning)

                tqc = qk.transpile(qc, basis_gates=basis_gates, optimization_level=transpile_optimization)
        
        # Get qasm instructions
        text = tqc.copy().qasm()  # .copy because .qasm() changes the gate names
        ignore = ['\n']
        qasm_instructions = []
        sentence = ''
        for c in text:
            if c == ';':
                qasm_instructions.append(sentence)
                sentence = ''
            elif c not in ignore:
                sentence += c
        if sentence != '':
            qasm_instructions.append(sentence)

        # Get instructions after `qreg`
        high_level_index = None
        for i, raw_instruction in enumerate(qasm_instructions):
            if 'qreg' in raw_instruction:
                high_level_index = i + 1
                break
        high_level_instructions = qasm_instructions[high_level_index:]

        # Parse high_level_instructions, reorder, transport, translate to pulses
        hfgui_pulses = []
        qubit_0_pulses = []
        qubit_1_pulses = []
        
        # Translate instructions to hfgui pulses
        n_ms_gates = 0
        for instruction in high_level_instructions:
            
            # MS Gates
            if 'q[0]' in instruction and 'q[1]' in instruction:

                if 'rxx' not in instruction:
                    continue
                
                if qubit_0_pulses != []:
                    hfgui_pulses.append('inline ion_1_potential();')
                    for pulse in qubit_0_pulses:
                        hfgui_pulses.append(pulse)
                    qubit_0_pulses = []
                
                if qubit_1_pulses != []:
                    hfgui_pulses.append('inline ion_2_potential();')
                    for pulse in qubit_1_pulses:
                        hfgui_pulses.append(pulse)
                    qubit_1_pulses = []
                
                if do_ms_gate:
                    hfgui_pulses.append('inline ms_potential();')
                    hfgui_pulses.append('inline ms_gate(0.);')
                    n_ms_gates+=1
                else:
                    hfgui_pulses.append('inline identity();')
            
            # Identities
            elif 'id' in instruction:
                continue

            # SIA qubit 0
            elif 'q[0]' in instruction:
                theta, phi = self.get_from_parenthesis(instruction)
                theta, phi = self.adapt_phase_for_negative_rotations(theta, phi)
                theta = self.SIA_pitime(theta, ion=0)
                theta = self.put_the_dot(theta)
                if aczs_comp and n_ms_gates > 0:
                    phi += ' + ACZS*gate_time*360.*%d.' % n_ms_gates
                phi = self.put_the_dot(phi)
                pulse = 'inline rot_1(%s, %s);' % (phi, theta)
                qubit_0_pulses.append(pulse)
            
            # SIA qubit 1
            elif 'q[1]' in instruction:
                theta, phi = self.get_from_parenthesis(instruction)
                theta, phi = self.adapt_phase_for_negative_rotations(theta, phi)
                theta = self.SIA_pitime(theta, ion=1)
                theta = self.put_the_dot(theta)
                if aczs_comp and n_ms_gates > 0:
                    phi += ' + ACZS*gate_time*360.*%d.' % n_ms_gates
                phi = self.put_the_dot(phi)
                pulse = 'inline rot_2(%s, %s);' % (phi, theta)
                qubit_1_pulses.append(pulse)
        
        # Catch the pulses after the last two-qubit gate
        if qubit_0_pulses != []:
            hfgui_pulses.append('inline ion_1_potential();')
            for pulse in qubit_0_pulses:
                hfgui_pulses.append(pulse)
            qubit_0_pulses = []
        if qubit_1_pulses != []:
            hfgui_pulses.append('inline ion_2_potential();')
            for pulse in qubit_1_pulses:
                hfgui_pulses.append(pulse)
            qubit_1_pulses = []
        
        return hfgui_pulses, tqc


    def compile(self, basis_gates=['rxx', 'r', 'id'], do_ms_gate=True, aczs_comp=True, transpile=True, transpile_optimization=0, n_cores=4):

        # Get rid of (SymPy) deprecation warnings
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)

            # Use multiprocessing
            if n_cores > 0:
                compile_circuit_args = [(circuit, basis_gates, do_ms_gate, aczs_comp, transpile, transpile_optimization) for circuit in self.circuits]
                with Pool(n_cores) as p:
                    print(compile_circuit_args)
                    result = p.starmap(self.compile_circuit, compile_circuit_args)
                for sub_result in result:
                    self.hfgui_sequences.append(sub_result[0])
                    self.hfgui_circuits.append(sub_result[1])

            # Do not use multiprocessing
            else:
                for circuit in self.circuits:
                    hfgui_pulses, tqc = self.compile_circuit(circuit, basis_gates, do_ms_gate, aczs_comp, transpile, transpile_optimization)
                    self.hfgui_sequences.append(hfgui_pulses)
                    self.hfgui_circuits.append(tqc)


    def decompile_sequence(self, sequence):
        qc = qk.QuantumCircuit(2)
        for pulse in sequence:
            # Get qubit_idx
            if 'ms_gate' in pulse:
                qubit_idx = 2
            elif 'rot_1' in pulse:
                qubit_idx = 0
            elif 'rot_2' in pulse:
                qubit_idx = 1
            else:
                continue

            if qubit_idx == 2:
                qc.rxx(pi/2, 0, 1)
            else:
                phi, theta = self.get_from_parenthesis(pulse)
                scan_pi_0 = pi
                scan_pi_1 = pi
                exec('theta_expr = ' + theta ,locals(), globals())
                exec('phi_expr = ' + phi, locals() ,globals())
                theta = theta_expr
                phi = phi_expr
                qc.r(theta, phi, qubit_idx)
        return qc

    def decompile(self, n_cores=4):
        with Pool(n_cores) as p:
            self.hfgui_circuits = p.map(self.decompile_sequence, self.hfgui_sequences)
            

    def write_files(self, out_folder='sequences', remove_old_files=False, prefix='seq', filenames=None, write_circuits=True):
        
        # Allow empty prefix to be nice
        if prefix != '':
            prefix += '_'
       
        # Create output folder 
        out_folder = Path(out_folder)
        if not os.path.isdir(out_folder):
            os.mkdir(out_folder)
        
        # Remove file inside output folder
        elif remove_old_files:
            for fname in os.listdir(out_folder):
                fname = Path(fname)
                if not os.path.isdir(fname):
                    os.remove(out_folder/fname)
    
        for i, sequence in enumerate(self.hfgui_sequences):
            
            # Check if cycle benchmarking was done
            if self.legends_l_k_m != [] and filenames is None:
                l, k, m = self.legends_l_k_m[i]
                final_state = self.final_states[i]
                filename = Path('seq_k_%02d_m_%d_l_%02d_f_%d.dc' % (k, m, l, final_state))
            
            # Assemble filenames using prefix and automatic counter
            elif filenames is None: 
                if self.final_states != []:
                    final_state = self.final_states[i]
                    filename = Path('%s%03d_f_%d.dc' % (prefix, i, final_state))
                else:
                    filename = Path('%s%03d.dc' % (prefix, i))
            
            # Use specified filename list
            else:
                filename = filenames[i]
                
            # Write files
            with open(out_folder/filename, 'w') as f:

                # Write circuit as comment
                if write_circuits:
                    circuit_text = str(self.hfgui_circuits[i].draw(output='text'))
                    line = ''
                    for c in circuit_text:
                        if c in ['\n', '\r']:
                            print('//' + line, file=f)
                            line = ''
                        else:
                            line += c
                    if line != '':
                        print('//' + line, file=f)

                # Write pulses
                for pulse in sequence:
                    print(pulse, file=f)


def main():
    pulser = Pulser()
    pulser.cycle_benchmark(m_list=[4, 8], L=5)
    pulser.compile()
    pulser.write_files()


def main_debug_decompile():
    pulser = Pulser()
    pulser.cycle_benchmark(m_list=[2], L=1)
    pulser.compile()
    pulser.decompile()
    final_states_alt = []
    for qc in pulser.circuits:
        state = pulser.svsim.run(qk.transpile(qc, pulser.svsim)).result().get_statevector()
        final_state_alt = int(np.around(np.abs(state[0]), 0))
        final_states_alt.append(final_state_alt)
    print(final_states_alt)
    print(pulser.final_states)
    


if __name__ == '__main__':
    main()

