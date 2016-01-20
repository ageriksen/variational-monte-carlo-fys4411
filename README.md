# variational-monte-carlo-fys4411

Example class structure for the first VMC project of FYS4411 (spring 2016). You may, if you wish, fork this repository and make it the basis of your project. If you choose to do this, you will have to implement a lot of the crucial functions yourself. The relevant functions you need to implement are spread throughout the project, but they are all commented with a note saying what each function should do.

Please note that this is only a start, and when you have implemented all of these functions you will only have completed the first exercise. However, once this is done, you will have a very good basis for further work, and adding functionality will be easier with a good class structure.

If you want to write your own code from scratch, you are of course welcome to do so, and feel free to use this code as inspiration for your own class structure.

- If you choose to use this code as a basis for your work, the first thing you should do is fork it, pull it down to your computer, and make sure it compiles and runs. After this you should spend at least 10 minutes looking at the structure and familiarizing yourself with how the classes interact with eachother. 
- A good way to do this may be to simply start at the top of the main.cpp file, and go through all the calls to the System class functions. Consider also the base classes WaveFunction, Hamiltonian, and InitialState and see which functions are virtual (which functions NEED to be implemented by any sub-class).
- You can skip over the output function in the Sampler class and the entire Random class.
- Next, you should start implementing the Gaussian wave function class.

