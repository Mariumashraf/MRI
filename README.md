
<img src = "/images/f6.PNG" width = "50%">


# 1- Bloch Equation: Bloch equations are a set of macroscopic equations that are used to calculate the nuclear magnetization M = (Mx, My, Mz) as a function of time when relaxation times T1 and T2 are present.
# 2- Function that rotates the bulk magnetization vector:
# The net magnetization vector in MRI is the summation of all the magnetic moments of the individual hydrogen nuclei and the bulk magnetization vector M, defined as the volume density of the spin magnetic moments, is aligned and proportional to the magnetic field. Its intensity characterizes the imbalance between the two spin populations.
# -We started rotation and we get Propagation Matrix first then we started magnetization with Keeping track of magnetization at all-time points then we simulate the Decay and the graph below shows the Result.

<img src = "/images/f1.PNG" width = "50%">


# In this graph we made a simple sequence simply consists of 60-degree excitation RF pulses about the y-axis, spaced TR apart and we plot the magnetization and decay of x, y and z

<img src = "/images/f2.PNG" width = "50%">


# In the following graphs we can load matlab's default set of MR image or browse image then we compute the Fourier transform for each image and then we allocate memory for k-space domain as K-space is the Fourier transform of the MR images.

<img src = "/images/f3.PNG" width = "50%">


# This graph shows a function that simulates the uniformity effect of magnetic field where its effect is imposed in the axial (z), but no radial (x or y) dependence on the magnetic field strength, the graph shows a plot of Bz(z).

<img src = "/images/f4.PNG" width = "50%">


# This graph shows the difference in angular frequency using heat map depending on the change in the Static Magnetic Field that was already calculated

<img src = "/images/f5.PNG" width = "50%">


# In the following graphs we considered Phase and Frequency Encoding to make a phantom made of oil and water and costructed in Kspace with gradient function and then we reconstructed the Image using inverse fourir transform.

<img src = "/images/f7.PNG" width = "50%">


<img src = "/images/f8.PNG" width = "50%">

