
Mesh.Algorithm = 6;

//h = 1.5;   /// approx 1k tet : 1054 tetrahera and 526 triangles

h = 100;   /// approx 1k tet : 10589 tetrahera and 2884 triangles

//h = 0.24;   /// approx 1k tet : 101524 tetrahera and 13436 triangles




// Dimensions: x is left to right, y is front to back and z is top to bottom
bottom_siz = 15 * h;
fracture_left_siz = 12 * h;
size_between_layers_left = 25 * h;
left_siz = bottom_siz;
right_siz = 35 * h;
size_between_layers_right = 12 * h;
fracture_right_siz = 12 * h;
inlet_siz = 12 * h;

// Bounding box points
Point(1) = {0.0, 0.0, 0.0, bottom_siz};
Point(2) = {100.0, 0.0, 0.0, bottom_siz};
Point(3) = {100.0, 100.0, 0.0, bottom_siz};
Point(4) = {0.0, 100.0, 0.0, bottom_siz};
Point(5) = {0.0, 0.0, 100.0, left_siz};
Point(6) = {100.0, 0.0, 100.0, right_siz};
Point(7) = {100.0, 100.0, 100.0, right_siz};
Point(8) = {0.0, 100.0, 100.0, left_siz};


// Lower layer outlet points
Point(9) = {0.0, 0.0, 10.0, size_between_layers_left / 2};
Point(10) = {100.0, 0.0, 10.0, size_between_layers_right};
Point(11) = {100.0, 100.0, 10.0, size_between_layers_right};
Point(12) = {0.0, 100.0, 10.0, size_between_layers_left};

// Fault boundary points
Point(13) = {0.0, 0.0, 80.0, fracture_left_siz};
Point(14) = {100.0, 0.0, 20.0, fracture_right_siz};
Point(15) = {100.0, 100.0, 20.0, fracture_right_siz};
Point(16) = {0.0, 100.0, 80.0, fracture_left_siz};

// Upper layer inlet points
Point(17) = {0.0, 0.0, 90.0, inlet_siz};
Point(18) = {0.0, 100.0, 90.0, inlet_siz};

//  layer inlet points JOSE
Point(19) = {100.0, 0.0, 90.0, inlet_siz};
Point(20) = {100.0, 100.0, 90.0, inlet_siz};

Point(21) = {100.0, 0.0, 80.0, fracture_left_siz};
Point(22) = {100.0, 100.0, 80.0, fracture_left_siz};
Point(23) = {0.0, 0.0, 20.0, fracture_right_siz};
Point(24) = {0.0, 100.0, 20.0, fracture_right_siz};

// Lines layer one vertical discretization
Line(1) = {1, 9};
Line(2) = {4, 12};
Line(3) = {2, 10};
Line(4) = {3, 11};

// Lines layer two vertical discretization
Line(5) = {10, 14};
Line(6) = {11, 15};
Line(7) = {16, 24};
Line(8) = {13, 23};

Line(9) = {5, 17};
Line(109) = {17, 13};
Line(10) = {8, 18};
Line(110) = {18, 16};
Line(11) = {14, 21};
Line(12) = {15, 22};

// Discretization in x-direction
Line(13) = {1, 2};
Line(14) = {9, 10};
Line(15) = {13, 14};
Line(16) = {5, 6};
Line(17) = {4, 3};
Line(18) = {12, 11};
Line(19) = {16, 15};
Line(20) = {8, 7};

// Discretization in y-direction
Line(21) = {2, 3};
Line(22) = {10, 11};
Line(23) = {14, 15};
Line(24) = {6, 7};
Line(25) = {1, 4};
Line(26) = {9, 12};
Line(27) = {13, 16};
Line(28) = {5, 8};
Line(112) = {17, 18};
Line(113) = {19, 20};
Line(114) = {19, 6};
Line(115) = {7, 20};
Line(116) = {20, 18};
Line(117) = {17, 19};


//JOSE
Line(118) = {21, 13};
Line(119) = {21, 22};
Line(120) = {22, 16};
Line(121) = {19, 21};
Line(122) = {20, 22};

Line(123) = {15, 24};
Line(124) = {23, 24};
Line(125) = {23, 14};

Line(126) = {9, 23};
Line(127) = {12, 24};//+
Curve Loop(1) = {20, -24, -16, 28};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {112, -116, -113, -117};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {118, 27, -120, -119};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {125, 23, 123, -124};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {14, 22, -18, -26};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {13, 21, -17, -25};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {15, 23, -19, -27};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {24, 115, -113, 114};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {121, 119, -122, -113};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {11, 119, -12, -23};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {5, 23, -6, -22};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {3, 22, -4, -21};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {112, 110, -27, -109};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {28, 10, -112, -9};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {8, 124, -7, -27};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {124, -127, -26, 126};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {26, -2, -25, 1};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {20, 115, 116, -10};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {116, 110, -120, -122};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {120, 19, 12};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {19, 123, -7};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {123, -127, 18, 6};
//+
Plane Surface(22) = {22};
//+
Curve Loop(23) = {18, -4, -17, 2};
//+
Plane Surface(23) = {23};
//+
Curve Loop(24) = {16, -114, -117, -9};
//+
Plane Surface(24) = {24};
//+
Curve Loop(25) = {117, 121, 118, -109};
//+
Plane Surface(25) = {25};
//+
Curve Loop(26) = {118, 15, 11};
//+
Plane Surface(26) = {26};
//+
Curve Loop(27) = {15, -125, -8};
//+
Plane Surface(27) = {27};
//+
Curve Loop(28) = {125, -5, -14, 126};
//+
Plane Surface(28) = {28};
//+
Curve Loop(29) = {14, -3, -13, 1};
//+
Plane Surface(29) = {29};
//+
Surface Loop(1) = {6, 29, 12, 23, 17, 5};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {5, 11, 28, 16, 22, 4};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {4, 7, 27, 15, 21};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {7, 26, 10, 20, 3};
//+
Volume(4) = {4};
//+
Surface Loop(5) = {19, 13, 25, 9, 3, 2};
//+
Volume(5) = {5};
//+
Surface Loop(6) = {2, 1, 18, 8, 24, 14};
//+
Volume(6) = {6};

//+
Physical Volume("RockMatrix_1") = {6, 5, 3, 2,4};

//+
Physical Volume("RockMatrix_2") = {1};


Coherence;

Transfinite Line "*" = 2 Using Bump 1;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";//+



