# Unscented Kalman Filter Project
Self-Driving Car Engineer Nanodegree Program

---

## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt [mode]`

   mode is optional 0, 1, or 2. 0=Use both laser and radar (default), 1=Laser only, 2=Radar only
   You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt output.txt`
    - eg. `./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt output.txt 1`
    - eg. `./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt output.txt 2`

## Project Instructions and Rubric

This information is only accessible by people who are already enrolled in Term 2
of CarND. If you are enrolled, see [the project page](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/c3eb3583-17b2-4d83-abf7-d852ae1b9fff/concepts/4d0420af-0527-4c9f-a5cd-56ee0fe4f09e)
for instructions and the project rubric.
