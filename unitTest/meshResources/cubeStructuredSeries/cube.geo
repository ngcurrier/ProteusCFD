//+
h = 0.5;
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {0, 1, 0, h};
Point(4) = {1, 1, 0, h};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1} = {3, 1, 2, 4};
//+
//Recombine Surface {1};
//+
Extrude {0, 0, 1} {
  Surface{1}; Layers{(1.0/h)}; Recombine;
}
