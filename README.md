# OpenWME

OpenWME is the first implementation of a general-purpose electromagnetic field solver (EM Solver), that is, software for simulation of the physics of electromagnetic fields, based on Weber-Maxwell electrodynamics. By using the modern Weber-Maxwell electrodynamics, OpenWME is not only orders of magnitude faster than any other EM solver that presently exists, but at the same time provides more accurate results with fewer numerical errors. The aforementioned advantages result from the fact that in Weber-Maxwell electrodynamics it is not necessary to solve Maxwell's equations numerically. Complex electromagnetic fields, waves and their reflections can therefore usually be calculated in a few seconds, while conventional solvers need many hours or even days for similar tasks.

------

## Weber-Maxwell electrodynamics

Key element of Weber-Maxwell electrodynamics is a force formula which describes the electromagnetic force between two point charges and which satisfies Newton's third law, i.e. *actio* = *reactio*. Unlike the well-known Coulomb force, however, the Weber-Maxwell force is universal and can be also applied to moving and arbitrarily accelerated point charges. Thus, compared to the Coulomb force, it contains all electromagnetic effects, e.g. magnetism or induction. 

However, Weber-Maxwell electrodynamics is also excellently suited for the modeling and analysis of electromagnetic waves. In this respect it goes beyond the limits of Weber electrodynamics, which is included as a subset. The relation to Maxwell electrodynamics is that the Weber-Maxwell force is the general solution of Maxwell's equations for point charges, provided that one interprets Maxwell's equations in an adequate way. 

The theoretical basics are described here: https://www.techrxiv.org/articles/preprint/24087840

![](examples/02_electromagnetic_waves/03_interference_at_a_double_slit/interference.gif)

OpenWME is currently in the alpha stage.

## Usage

OpenWME consists of a C++ library and a number of example applications that include the library at source code level.

For most of the examples there is a video demonstrating the result of the simulation:

- Magnetism: [Force between two long straight wires](examples/01_quasistatics/01_magnetic_force_between_wires/magnetic_force_between_wires.webm?raw=true)
- Electromagnetic induction: [Generation of a field near a moving conductor loop with direct current](examples/01_quasistatics/02_moving_current_loop_dc_current/moving_current_loop_dc_current.webm?raw=true)
- Lorentz force: [Moving point charge in the field of a Helmholtz coil](examples/01_quasistatics/03_point_particle_in_a_helmholtz_coil/point_particle_in_a_helmholtz_coil.webm?raw=true)
- Electromagnetic induction: [Generation of a field in the vicinity of a stationary conductor loop with alternating current](examples/01_quasistatics/04_current_loop_with_ac_current/current_loop_with_ac_current.webm?raw=true)
- Electromagnetic waves: [Field of a very fast moving electromagnetic transmitter](examples/02_electromagnetic_waves/01_moving_hertzian_dipole/moving_hertzian_dipole.webm?raw=true)
- Reflection of electromagnetic waves: [Electromagnetic pulse reflected back and forth between two point charges](examples/02_electromagnetic_waves/02_reflection_two_resting_hertzian_dipols/reflection_two_resting_hertzian_dipols.webm?raw=true)
- Interference: [Interference of an electromagnetic wave at a double slit](examples/02_electromagnetic_waves/03_interference_at_a_double_slit/interference_at_a_double_slit.webm?raw=true)
- Special relativity: [Einstein's light clock](examples/03_special_relativity/01_light_clock/light_clock.webm?raw=true)
- Relativity of simultaneity: [Each intertial frame has its own fields](examples/03_special_relativity/02_reflection_moving_transmitter_resting_receiver/reflection_moving_transmitter_resting_receiver.webm?raw=true)
- Quantum forces: [Field of the force of a Hertzian dipole on itself when it is in front of a double slit as a function of the location](examples/04_quantum_mechanics/01_quantum_forces_at_a_double_slit/quantum_forces_at_a_double_slit.png?raw=true)

For the visualization of the results, [Cairo](https://www.cairographics.org/) is currently used. However, the library can also be used without graphical elements.

## Building from source

The development is done on Linux [Ubuntu](https://ubuntu.com/) or [Debian](https://packages.debian.org/), using [GNU Make](https://www.gnu.org/software/make/) to compile the applications. To build OpenWME, the packages `build-essential` and `libcairo2-dev` are needed.

## Licensing

OpenWME is a free, open source software. The source code is available under the GPL v3 license.
