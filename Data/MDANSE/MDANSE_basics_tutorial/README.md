# MDANSE Basics

This tutorial requires us to analyse simple trajectories.
We are going to generate those by running scripts.

## Requirements

Please install the following Python packages:

* `mdanse[all]`
* `mdanse_gui[all]`
* `ase`
* `asap3`

## Basic MDANSE tutorial

We will follow the slides which explain different tasks to be
completed within this tutorial.

### Scripts

These scripts will use ASE and the ASAP calculator to generate
trajectories of a binary alloy of copper and gold.

* run_isotropic_trajectory.py
* run_anisotropic_trajectory.py

### Answers

Some of the results are stored in the "answers" directory.
You are welcome to have a look at them if you already completed the
tutorial tasks, but want to verify the answers.

Since many tasks in this tutorial deal with different functions of
the graphical interface, it is not possible to save the outcome of
many of them as answers.

## Trajectory stability exercise

Once the previous part has been completed, you can move on to working
on the trajectory called `tutorial1_trajectory.h5`.
It is a trajectory in H5MD format, and your task is to find out what
this file contains and what it can be used for, or not.

You are encouraged to run different analysis types on limited subsets
of trajectory frames and plot the results together to check if the
results changed within the simulation.
