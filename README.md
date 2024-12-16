This project delves into the world of computational quantum physics, focusing on calculating the ground state energy of a **hydrogen molecule (H2)**. It leverages the power of the **Born-Oppenheimer approximation**, a fundamental concept in molecular physics, and employs ** Quantum Monte Carlo (QMC)** techniques to achieve its goal. This GitHub repository houses the Python code developed for this project. It is a compilation of the different python files used, and would have to be separated into different files to emulate the results.

**Theoretical Foundation and Methodology**

At the heart of the project lies the application of Monte Carlo methods for evaluating complex multi-dimensional integrals that arise in quantum mechanical calculations. The code explores and compares the effectiveness of different Monte Carlo techniques, including: 

* **Basic Monte Carlo Integration:** This method relies on random sampling to estimate the value of an integral.
* **Importance Sampling:**  This technique enhances the efficiency of Monte Carlo integration by focusing the sampling in regions where the integrand is most significant.
* **Metropolis Method:** This approach generates a sequence of random samples (electron positions in our case) that are distributed according to a desired probability distribution, which is determined by the trial wavefunction.

The core methodology consists of the following key steps:

* **Choosing a Trial Wavefunction:**  The first step involves selecting a physically reasonable trial wavefunction to represent the H2 molecule's ground state. This trial function, a product of individual atomic orbitals and a correlation term,  incorporates parameters (β and proton separation, S) that are subsequently optimized. 
* **Variational Method for Optimization:** The variational method is applied to determine the optimal values of β and S that minimize the energy expectation value calculated using the trial wavefunction. This optimization ensures that the trial wavefunction provides the best possible approximation to the true ground state within its functional form.
* **Imaginary Time Evolution:**  The project implements imaginary time evolution as a powerful refinement technique. By treating the time-dependent Schrödinger equation as a diffusion equation, the code iteratively updates the wavefunction to converge towards the true ground state.  This step significantly improves the accuracy of the energy calculation.
* **Extension to Helium Atom:** To highlight the versatility of the implemented methods, the code includes a simulation of a helium atom. This is achieved by modifying the proton separation parameter, effectively treating the two protons of the helium nucleus as a single entity. 

**Key Results and Insights**

The project successfully calculates the ground state energy of the H2 molecule and achieves the following notable results:

* **Binding Energy (Eb) for H2:** The project obtains a binding energy of **Eb = 4.171 ± 0.088eV** for the hydrogen molecule. This result aligns well with values obtained from other studies employing the Born-Oppenheimer approximation, which typically yield a binding energy of approximately 4.25 eV. This agreement confirms the effectiveness and accuracy of the implemented QMC methods. 
* **Proton Separation (S) in H2:** The calculated proton separation, **S = 0.77 ± 0.02 Å**,  demonstrates good agreement with the experimentally accepted value of 0.74 Å.  This accurate determination of the equilibrium bond length further validates the reliability of the computational approach. 
* **Local Energy (E) for Helium:**  The simulation of the helium atom results in a local energy of **E = −78.32 ± 0.11 eV**. While this value is close to the experimental ground state energy of -78.9 eV, it's important to note that this simulation does not include the imaginary time evolution step. Applying imaginary time evolution to the helium atom case would likely lead to an energy value that is even closer to the experimental ground state energy.
