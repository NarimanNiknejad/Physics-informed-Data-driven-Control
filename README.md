# Physics-informed Data-driven Safe and Optimal Control Design
This repository contains MATLAB codes for implementing safe and optimal physics-informed data-driven controls for "Physics-informed Data-driven Safe and Optimal Control Design" submitted to IEEE Control Systems Letters.
[![DOI](https://img.shields.io/badge/DOI-10.1234/abcd-blue)](https://doi.org/10.1109/LCSYS.2023.3333257)


# Abstract
This paper introduces a novel approach that is based on a physics-informed data-driven approach to design a robust controller for discrete-time linear time-invariant systems. The goal is to enhance both the performance and feasibility of safe and optimal control design, all without the need for explicit system modeling. The idea is to start from a robust control design when no data samples are available using physics information and progressively move towards a fully adaptive controller as more data becomes available. This will enhance the feasibility of designing a safe control system and elevate the performance of an optimal control system. It achieves this by integrating safety and performance specifications for systems that lie at the intersection of two information sets: the physics-informed set of possible system models and the data-conformity set of models. This intersection set is non-empty if the prior knowledge includes the actual system model. Besides, it is smaller than physics-informed and data-conformity sets. Linear Matrix inequality conditions are provided to robustly satisfy the safety and performance of the systems that fall at the intersection set. Two applications are presented to verify the theoretical results: the safe hovering of a quadcopter and the quadratic stabilization of a Lithium-ion battery.


## Usage

To run the MATLAB code, follow these steps:

1. Clone this repository to your local machine.

2. Open the MATLAB script in your MATLAB environment.

3. Make any necessary parameter adjustments in the script.

4. Run the script to simulate and visualize the results.

## Dependencies

This code relies on MATLAB.

## Obtaining and Licensing MOSEK for MATLAB

MOSEK is a powerful optimization solver that can be used with MATLAB to solve various mathematical optimization problems. To get started with MOSEK for MATLAB, follow these steps:

### 1. Download MOSEK

1. Visit the MOSEK download page: [MOSEK Download](https://www.mosek.com/downloads/).

2. Select the appropriate version of MOSEK for your operating system. MOSEK provides versions for Windows, Linux, and macOS.

3. Download the MOSEK installation package.

### 2. Install MOSEK

Follow the installation instructions provided by MOSEK to install the software on your system.

### 3. Obtain a License

1. MOSEK requires a license to use. You can request a free academic license, a trial license, or purchase a commercial license.

2. To request an academic license or a trial license, visit the MOSEK License Request page.

3. Follow the steps on the license request page to obtain your license file. This file will be used to activate MOSEK on your machine.

4. If you decide to purchase a commercial license, contact MOSEK directly through their website for more information on pricing and licensing options.

### 4. Configure MOSEK for MATLAB

Once you have installed MOSEK and obtained a valid license, you need to configure it for MATLAB:

1. Locate the MOSEK installation directory on your system.

2. In MATLAB, set up the path to include the MOSEK directory. 


## License and Contact Info

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. You can customize the parameters and use them for your specific control system applications.

If you have any questions or encounter issues, please feel free to [contact me](mailto:niknejad@msu.edu).

Enjoy exploring and using the code!



Watch the video of the safe hovering of the drone here:

[![Watch the video](https://img.youtube.com/vi/LdfYQQp4STU/maxresdefault.jpg)](https://www.youtube.com/watch?v=LdfYQQp4STU)



