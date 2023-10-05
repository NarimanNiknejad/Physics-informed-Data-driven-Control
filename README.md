# Physics-informed Data-driven Safe and Optimal Control Design

# Safe Physics Informed Data-Driven Control in MATLAB

This repository contains MATLAB code for implementing safe physics-informed data-driven control. The code includes the following components:

## Code Sections

### 1. Define general parameters

- `nStates`, `nInputs`, `nRealization`, `nSample`, `lambda`, `P_safe`, and `x0` are general parameters defined for the control problem.

### 2. Time sequence

- `T0`, `Ts`, `Tf`, and `t` define the time sequence for simulation.

### 3. Noise generation

- Noise with specific statistics is generated for use in the simulation.

### 4. Define the parameters related to the source system

- `As` and `Bs` define the system dynamics.

### 5. The difference between the actual system and the target system

- `ep_AB`, `delta_A`, and `delta_B` represent differences between the actual and target systems.

### 6. Define the parameters of the target system

- `At` and `Bt` represent the target system dynamics.

### 7. Define parameters for uncertain system

- Parameters related to the uncertain system are defined here.

### 8. Define parameters for collected data

- Data is collected from the target system and used for further calculations.

### 9. Define parameters for data-based uncertainty

- Parameters for data-based uncertainty are defined in this section.

### 10. Define parameters for lambda_contracted ellipsoid

- Parameters for the lambda-contracted ellipsoid are defined here.

### 11. Main loop

- The main simulation loop is defined here, including the control law.

### 12. Plot results

- Code for plotting simulation results is provided in this section.

## Usage

To run the MATLAB code, follow these steps:

1. Clone this repository to your local machine.

2. Open the MATLAB script in your MATLAB environment.

3. Make any necessary parameter adjustments in the script.

4. Run the script to simulate and visualize the results.

## Dependencies

This code relies on MATLAB.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


## Conclusion

This MATLAB code provides an implementation of safe physics-informed data-driven control. You can customize the parameters and use them for your specific control system applications.

If you have any questions or encounter issues, please feel free to [contact me](mailto:niknejad@msu.edu).

Enjoy exploring and using the code!



Watch the video of the safe hovering of the drone here:

[![Watch the video](https://img.youtube.com/vi/LdfYQQp4STU/maxresdefault.jpg)](https://www.youtube.com/watch?v=LdfYQQp4STU)



