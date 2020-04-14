OpenCV_Rectangle_Detection
==========================

Alters the default rectangle detection algorithm in OpenCV to filter out found rectangles.

### Installing opencv on Ubuntu:
Simple way:
sudo apt-get install libopencv-dev

### Compile on Ubuntu in project folder:
g++ squares.cpp -o squares `pkg-config --cflags --libs opencv`

### Help Example:
A program using pyramid scaling, Canny, contours, contour simpification and
memory storage (it's got it all folks) to find
squares in a list of images pic1-6.png
Returns sequence of squares detected on the image.
the sequence is stored in the specified memory storage
Call:
./squares
Using OpenCV version 3.2.0

