# Physics-informed Data-driven Safe and Optimal Control Design
In this project, we introduce a novel paradigm known as physics-informed data-driven control, presenting an approach aimed at steadily enhancing the performance and feasibility of controlled systems. Our methodology revolves around finite-adaptive control, which initiates with robust control design in the absence of any data samples and progressively evolves the controller into a fully adaptive form as more data becomes available. The fundamental premise is to ensure the satisfaction of critical system properties, such as safety and stability, within the intersection of two pivotal sets: the physics-informed set encompassing potential system models informed by prior knowledge and the data-conformed set representing models aligned with observed data. Importantly, this intersection is non-empty only when the prior knowledge incorporates the true system model, and it remains smaller in size compared to both the physics-informed and data-conformed sets, thereby enhancing the feasibility and performance of the controlled system. To facilitate these objectives, we provide Linear Matrix Inequality conditions that robustly guarantee the safety and performance of the controlled system within the intersection set, offering a powerful framework for the design of physics-informed data-driven control systems.


This repository contains MATLAB codes for implementing safe and optimal physics-informed data-driven control.

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

2. To request an academic license or a trial license, visit the MOSEK License Request page: [MOSEK License Request](https://www.mosek.com/license/).

3. Follow the steps on the license request page to obtain your license file. This file will be used to activate MOSEK on your machine.

4. If you decide to purchase a commercial license, contact MOSEK directly through their website for more information on pricing and licensing options.

### 4. Configure MOSEK for MATLAB

Once you have installed MOSEK and obtained a valid license, you need to configure it for MATLAB:

1. Locate the MOSEK installation directory on your system.

2. In MATLAB, set up the path to include the MOSEK directory. You can do this by running the following command in the MATLAB command window, replacing `/path/to/mosek` with the actual path to your MOSEK installation directory:

   ```matlab
   addpath('/path/to/mosek');


## License and Contact Info

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. You can customize the parameters and use them for your specific control system applications.

If you have any questions or encounter issues, please feel free to [contact me](mailto:niknejad@msu.edu).

Enjoy exploring and using the code!



Watch the video of the safe hovering of the drone here:

[![Watch the video](https://img.youtube.com/vi/LdfYQQp4STU/maxresdefault.jpg)](https://www.youtube.com/watch?v=LdfYQQp4STU)



