# Adjoint Optimization software for Lumerical from UC Berkeley

## Intro

This is a matlab implementation of a Continuous Adjoint Optimization method for Lumerical FDTD. This code was developped in the EECS department of UC Berkeley between 2008 and 2013 by sudents in Prof. Eli Yablonovitch's group. It is distributed under the MIT License. 

Presentations of the code here:
http://optoelectronics.eecs.berkeley.edu/PhotonicInverseDesign/index.html

Graduate Students having worked on this code:
  - Owen Miller
  - Samarth Bhargava
  - Vidya Ganapati
  - Christopher Lalau-Keraly
  
 ## Some Resources
 
 ### Thesis:
 
 	- Owen Miller: Photonic Design: From Fundamental Solar Cell Physics to Computational Inverse Design 	
 		http://optoelectronics.eecs.berkeley.edu/ThesisOwenMiller.pdf
	- Samarth Bhargava: Heat-Assisted Magnetic Recording: Fundamental Limits to Inverse Electromagnetic Design
		 http://optoelectronics.eecs.berkeley.edu/BhargavaDissertation.pdf
 	- Vidya Ganapati: Optical Design Considerations for High Conversion E fficiency in Photovoltaics
		http://optoelectronics.eecs.berkeley.edu/GanapatiDissertation.pdf
	- Christopher Lalau-Keraly: Optimizing Nanophotonics: from Photoreceivers to Waveguides	
		https://www2.eecs.berkeley.edu/Pubs/TechRpts/2017/EECS-2017-20.pdf
### Papers: 
 	- Adjoint shape optimization applied to electromagnetic design 
 	https://www.osapublishing.org/oe/abstract.cfm?uri=oe-21-18-21693
 	- Shape optimization of nanophotonic devices using the adjoint method
 	http://ieeexplore.ieee.org/document/6989970/
  
The examples provided have been tested and work with Matlab 2016a and Lumerical, on Mac OS High Sierra, and Red Hat 7. It should also work with Windows with a little tweaking.
 
 ## Installation 

Optimization code download:
In a terminal, go to the folder where you would like to install and run:
  `git clone https://github.com/chriskeraly/BerkeleyPhotonicInverseDesign.git`

Matlab and Lumerical must be installed.

The only tricky part in the installation is to get Lumerical and Matlab to get to talk to each other correctly. Both Matlab and Lumerical need to be able to launch each other. 
 
 An AWS AMI will be provided very shortly, so that you can get the software working almost 'out of the box'.
 
 ### Lumerical --> Matlab
 
 Follow these instructions: https://kb.lumerical.com/en/install_matlab_integration.html
 
 For Mac it is a little tricky to find the executable. After a simple installation on High Sierra, it is here: `/Applications/Lumerical/FDTD Solutions/FDTD Solutions.app/Contents/MacOS`
 For the Red Hat AWS AMI, if you follow the instructions given in the following section, it is here: `/usr/local/bin`
 
 ### Matlab --> Lumerical
 
 In `Optimization_software/runOpt_params.m` you must matlab where your Lumerical FDTD executable is. In linux, `which fdtd-solutions` in a terminal should give you the path. For a normal Mac install: `/Applications/Lumerical/FDTD\ Solutions/FDTD\ Solutions.app/Contents/MacOS/fdtd-solutions`
 
## First Optimization

Once you have set everything up correctly, open Matlab, go to `/Optimization_software` and run `runOpt.m`. A window should pop out, asking you to select a 'setupFile.m'. This is a matlab file containing all the information about the type of optimization we wish to run (Figure of Merit, Geometry representation, etc...).  Navigate to `/Optimization_software/examples/EZSiPh_compactsplitter` and choose `setup_EZSiPh_params.m`. Then you are asked to select a 'baseFile.fsp'. This is a Lumerical file which contains all the geometry information about your simulation. Navigate to `/Optimization_software/examples/EZSiPh_compactsplitter` and choose `base_compactsplitter.fsp`. Then, sit back and enjoy!

## AWS install (still work in progress)

In order to facilitate the setup and use of the software, an Amazon AMI machine is available. In a few steps, you can launch a powerfull machine on AWS with Lumerical and Matlab (almost...) pre-installed and preconfigured, which will be able to run these simulations for you. This is a good way to run optimizations, since you can simply get them started and then let them run without the need for your personal computer ot be always on. 

### Launching and connecting to the AMI. 

You will need an AWS account. Go to the EC2 dashboard and start a machine, in Oregon (this is where the AMI is stored). c4.large is a great starting point, but you can choose much bigger machines if you want. You will have to save the ssh key that is given to you in a safe location (.ssh/keys/ is a nice place), and change the access rights to it: `chmod 400 mykey.pem`. (All this info is given to you if you click `Connect` in the EC2 dashboard, although aws suggests using root, when you should really use ec2-user)

Once your machine is up, I strongly recommend setting up your ssh config file (located at ~/.ssh/config) and adding these lines to it (obviously replacing with the right values for your HostName and Identityfile)
```
Host BerkeleyOpt
	Hostname ec2-xxx-xxx-xxx-xxx.us-west-2.compute.amazonaws.com
	User ec2-user
	IdentityFile ~/.ssh/keys/mykey.pem
```
 
At this point check if your ssh access is working by typing `ssh BerkeleyOpt` in a terminal. If you connect: Congrats! You are in your new machine.

Now connect to it again using these commands instead: `ssh -v -C -L 5901:localhost:5901 BerkeleyOpt`. This will do port forwarding to your machine, so that you can use vnc. A vncserver should already be running on the machine, but if it doesnt, just type `vncserver` in your ssh session and that should boot it up. 

Using your favorite VNC client, connect to localhost:5901. (In macOs, there is a native client called ScreenSharing, you can even just type localhost:5901 in the address bar of safari and that should boot it up for you). If you see your linux Desktop, you can go to the next step

### Finishing the installation of matlab and Lumerical

#### Matlab
While the installation files for matlab are on the AMI, it is not installed yet: start a terminal in your vnc session, and type `sudo matlab/matlab_R2016a_glnxa64/install`, then follow the instructions and activate matlab with an appropriate license. 

#### Lumerical
Lumerical is pre-installed, but you must still get your license in. In a terminal, type `fdtd-solutions`, and then enter the correct license info for lumerical with your license server information. If you have a temporary license dedicated for this machine, then you can follow the instuctions here: https://kb.lumerical.com/en/install_linux_instructions_fdtd.html 

WARNING: Do not install a fixed license on this machine! If you do so, you will loose your license if you disconnect from it! You can start a small machine (t2.micro) with the same AMI, and create a license server there, although you will need to open a few ports in aws to make your licenses available to the outside world. 

#### Clone the repo

The branch with the settings for the aws install is appropriatly named `aws`
```bash
cd
git clone https://github.com/chriskeraly/BerkeleyPhotonicInverseDesign.git
git checkout aws
```

At this point your machine should be all set, and you can run your first optimization!
