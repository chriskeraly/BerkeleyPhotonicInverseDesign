### Adjoint Optimization software for Lumerical from UC Berkeley

# Intro

This is a matlab implementation of a Continuous Adjoint Optimization method for Lumerical FDTD. This code was developped in the EECS department of UC Berkeley between 2008 and 2013 by sudents in Prof. ELi Yablonovitch's group. It is distributed under the MIT License. 

Presentations of the code here:
http://optoelectronics.eecs.berkeley.edu/PhotonicInverseDesign/index.html

Graduate Students having worked on this code:
  -Owen Miller
  -Samarth Bhargava
  -Vidya Ganapati
  -Christopher Lalau-Keraly
  
The examples provided have been tested and work with Matlab 2016a and Lumerical, on Mac OS High Sierra, and Red Hat 7. It should also work with Windows with a little tweaking.
 
 # Installation 

Optimization code download:
In a terminal, go to the folder where you would like to install and run:
  `git clone https://github.com/chriskeraly/BerkeleyPhotonicInverseDesign.git`

Matlab and Lumerical must be installed.

The only tricky part in the installation is to get Lumerical and Matlab to get to talk to each other correctly. Both Matlab and Lumerical need to be able to launch each other. 
 
 An AWS AMI will be provided very shortly, so that you can get the software working almost 'out of the box'.
 
 ## Lumerical --> Matlab
 
 Follow these instructions: https://kb.lumerical.com/en/install_matlab_integration.html
 
 For Mac it is a little tricky to find the executable. After a simple installation on High Sierra, it was here: `/Applications/Lumerical/FDTD Solutions/FDTD Solutions.app/Contents/MacOS`
 For the Red Hat AWS AMI, if you follow the instructions given in the following section, it is here: `/usr/local/bin`
 
 ## Matlab --> Lumerical
 
 In `Optimization_software/runOpt_params.m` you want to tell matlab where you Lumerical FDTD executable is. The paths already there typical for Mac High Sierra or the Red Hat AWS AMI provided.
 
# First Optimization

sudo yum install git

## AWS install




# Connecting with aws

ssh into the machine, and start the vnc server: `vncserver`

SSH with 5901 port forwarding and start screen sharing

Open terminal on machine

cd /home/ec2-user/matlab/matlab_R2016a_glnxa64

sudo ./install

install only base matlab

