# Bridge Game

A little project to play around with the Finite Element Method.

## Prerequisites
### MacOS
This is the platform I've been testing on, so it should "just work". Try:
```
brew install bazel
bazel build :main
```
### Linux
This should mostly "just work". Make sure to install `xorg-dev" though for xll
```
sudo apt-get install bazel
sudo apt-get install xorg-dev
bazel build :main
```

## Control
To run:
```
bazel run :main
```

The UI should launch and will allow you to:
 - Click to draw a brick
 - Control click to erase a brick
 - Press `z` to increase the cursor size
 - Press `x` to reset the cursor size

When you press the space-bar, the simulation will start. Bricks that are under too much stress will be destroyed. Pressing the space-bar again will go back to building mode.

## Notes
 - The game is a work in progress without too much testing, so your mileage may very
 - I haven't tested on other platforms besides my Mac

## Sources
[1] http://www.unm.edu/~bgreen/ME360/2D%20Triangular%20Elements.pdf

[2] http://web.mit.edu/kjb/www/Books/FEP_2nd_Edition_4th_Printing.pdf
