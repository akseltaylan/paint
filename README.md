# paint
Drawing Program in C++ with OpenGL

0 - 7: Change color (black, red, green, yellow, blue, magenta, cyan, white)
L: Line drawing mode
P: Polygon drawing mode
    - left click to add a point, right click to indicate the final point of the polygon
F: Free-hand drawing mode
Q: Boundary fill mode
    - make sure current color matches boundary you are filling in
O: Flood fill mode
C: Clear canvas
R: increment R of current color by 0.05 (max of 1)
G: increment G of current color by 0.05 (max of 1)
B: increment B of current color by 0.05 (max of 1)

Known bugs:
-line drawing erases other colors in the canvas
-Boundary fill has a chance to stack overflow
