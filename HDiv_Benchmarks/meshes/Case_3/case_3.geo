Geometry.Tolerance = 3.33333333333e-05;

s = 1.0;

//h_domain = 0.125; /// 30k
//ng = 12; /// 30k

h_domain = 0.05; /// 30k
ng = 25; /// 30k


nf_0 = ng;
nf_1 = ng;
nf_2 = ng;
nf_3 = ng;
nf_4 = ng;
nf_5 = ng;
nf_6 = ng;
nf_7 = ng;




// Fracture 1, left (low x) and right ends
h_1_left = s*0.1;
h_1_right = s*0.1;

// Fracture 2
h_2 = s*0.1;

// Intersection of 1 and 2
h_1_2 = s*0.05;

// Intersection of fractures 1 and 3
h_1_3 = s*0.04;
// Other points on fracture 3
h_3 = s*0.1;

// Points on fracture 4, close to f1 and far away
h_4_close = s*0.04;
h_4_away = s*0.1;

// Endpoints of fractures 5 and 6. Not intersection (below)
h_5_6 = s*0.04;
// Intersection between 5 and 6
h_5_6_isect = s*0.03;

// Intersection of 1 with 5 and 6
h_1_5_6 = s*0.03;

// Fracture 7
h_7 = s*0.06;
// Fracture 8
h_8 = s*0.06;
// Intersection between 1 and 7 and 8
h_1_7 = s*0.06;
h_1_8 = s*0.06;

ymax = 2.25;

// Fracture 1 points
p0 = newp; Point(p0) = {0.05, 0.25, 0.5, 0.1 , h_1_left };
p1 = newp; Point(p1) = {0.95, 0.25, 0.5, 0.1 , h_1_left };
p2 = newp; Point(p2) = {0.95, 2.0, 0.5, 0.1 , h_1_right };
p3 = newp; Point(p3) = {0.05, 2.0, 0.5, 0.1 , h_1_right };

// Fracture 2 points
p4 = newp; Point(p4) = {0.5, 0.05, 0.95, 0.1 , h_2 };
p5 = newp; Point(p5) = {0.5, 0.05, 0.05, 0.1 , h_2 };
p6 = newp; Point(p6) = {0.5, 0.3, 0.05, 0.1 , h_2 };
p7 = newp; Point(p7) = {0.5, 0.3, 0.95, 0.1 , h_2 };

// Intersection of fracture 1 and fracture 3 points
p8 = newp; Point(p8) = {0.05, 1.0, 0.5, 0.04 , h_1_3 };
p9 = newp; Point(p9) = {0.95, 1.0, 0.5, 0.04 , h_1_3 };

// Other points of fracture 3
p10 = newp; Point(p10) = {0.95, 2.2, 0.85, 0.1 , h_3 };
p11 = newp; Point(p11) = {0.05, 2.2, 0.85, 0.1 , h_3 };

// Fracture 4 points
p12 = newp; Point(p12) = {0.05, 1.0, 0.48, 0.04 , h_4_close };
p13 = newp; Point(p13) = {0.95, 1.0, 0.48, 0.04 , h_4_close };
p14 = newp; Point(p14) = {0.95, 2.2, 0.14, 0.1 , h_4_away };
p15 = newp; Point(p15) = {0.05, 2.2, 0.14, 0.1 , h_4_away };

// Fractures 5 and 6 points
p16 = newp; Point(p16) = {0.23, 1.9, 0.3, 0.04, h_5_6 };
p17 = newp; Point(p17) = {0.23, 1.9, 0.7, 0.04, h_5_6 };
p18 = newp; Point(p18) = {0.17, 2.2, 0.7, 0.04, h_5_6 };
p19 = newp; Point(p19) = {0.17, 2.2, 0.3, 0.04, h_5_6 };
p20 = newp; Point(p20) = {0.17, 1.9, 0.3, 0.04, h_5_6 };
p21 = newp; Point(p21) = {0.17, 1.9, 0.7, 0.04, h_5_6 };
p22 = newp; Point(p22) = {0.23, 2.2, 0.7, 0.04, h_5_6 };
p23 = newp; Point(p23) = {0.23, 2.2, 0.3, 0.04, h_5_6 };

// Fracture 7 points
p24 = newp; Point(p24) = {.77, 1.9, 0.3, 0.06 , h_7 };
p25 = newp; Point(p25) = {.77, 1.9, 0.7, 0.06 , h_7 };
p26 = newp; Point(p26) = {.77, 2.2, 0.7, 0.06, h_7 };
p27 = newp; Point(p27) = {.77, 2.2, 0.3, 0.06, h_7 };

// Fracture 8 points
p28 = newp; Point(p28) = {0.83, 1.9, 0.3, 0.06, h_8 };
p29 = newp; Point(p29) = {0.83, 1.9, 0.7, 0.06, h_8 };
p30 = newp; Point(p30) = {0.83, 2.2, 0.7, 0.06, h_8 };
p31 = newp; Point(p31) = {0.83, 2.2, 0.3, 0.06 , h_8 };

// Domain corner points
p32 = newp; Point(p32) = {0.0, 0.0, 1.0, h_domain };
p33 = newp; Point(p33) = {0.0, 0.0, 0.0, h_domain };
p34 = newp; Point(p34) = {0.0, ymax, 0.0, h_domain };
p35 = newp; Point(p35) = {0.0, ymax, 1.0, h_domain };
p36 = newp; Point(p36) = {1.0, 0.0, 1.0, h_domain };
p37 = newp; Point(p37) = {1.0, 0.0, 0.0, h_domain };
p38 = newp; Point(p38) = {1.0, ymax, 0.0, h_domain };
p39 = newp; Point(p39) = {1.0, ymax, 1.0, h_domain };

// Intersection of Fractures 1 and 2 points
p40 = newp; Point(p40) = {0.5, 0.3, 0.5, 0.05 , h_1_2 };
p41 = newp; Point(p41) = {0.5, 0.25, 0.5, 0.05 , h_1_2 };

// Intersections between fractures 1 and 5 and 6 points
p42 = newp; Point(p42) = {0.23, 1.9, 0.5, 0.03 , h_1_5_6 };
p43 = newp; Point(p43) = {0.21, 2.0, 0.5, 0.03 , h_1_5_6 };
p44 = newp; Point(p44) = {0.17, 1.9, 0.5, 0.03 , h_1_5_6 };
p45 = newp; Point(p45) = {0.19, 2.0, 0.5, 0.03 , h_1_5_6 };

// Intersections between fractures 1 and 7 and 8 points
p46 = newp; Point(p46) = {0.77, 1.9, 0.5, 0.06, h_1_7 };
p47 = newp; Point(p47) = {0.77, 2.0, 0.5, 0.06, h_1_7 };
p48 = newp; Point(p48) = {0.83, 1.9, 0.5, 0.06, h_1_8 };
p49 = newp; Point(p49) = {0.83, 2.0, 0.5, 0.06, h_1_8 };

// Intersection of fractures 5 and 6 points
p50 = newp; Point(p50) = {0.2, 2.05, 0.7, 0.03 , h_5_6_isect };
p51 = newp; Point(p51) = {0.2, 2.05, 0.3, 0.03 , h_5_6_isect };

// Domain boundaries conditions points
pin00 = newp; Point(pin00) = {0., 0., 0.333333, h_domain};
pin01 = newp; Point(pin01) = {0., 0., 0.666667, h_domain};
pin11 = newp; Point(pin11) = {1., 0., 0.666667, h_domain};
pin10 = newp; Point(pin10) = {1., 0., 0.333333, h_domain};
pout00 = newp; Point(pout00) = {0., ymax, 0.333333, h_domain};
pout01 = newp; Point(pout01) = {0., ymax, 0.666667, h_domain};
pout11 = newp; Point(pout11) = {1., ymax, 0.666667, h_domain};
pout10 = newp; Point(pout10) = {1., ymax, 0.333333, h_domain};

// End of point specification

// Define lines 
frac_line_0= newl; Line(frac_line_0) = {p0, p8};
frac_line_1= newl; Line(frac_line_1) = {p0, p41};
frac_line_2= newl; Line(frac_line_2) = {p1, p9};
frac_line_3= newl; Line(frac_line_3) = {p1, p41};
frac_line_4= newl; Line(frac_line_4) = {p2, p9};
frac_line_5= newl; Line(frac_line_5) = {p2, p49};
frac_line_6= newl; Line(frac_line_6) = {p3, p8};
frac_line_7= newl; Line(frac_line_7) = {p3, p45};
frac_line_8= newl; Line(frac_line_8) = {p4, p5};
frac_line_9= newl; Line(frac_line_9) = {p4, p7};
frac_line_10= newl; Line(frac_line_10) = {p5, p6};
frac_line_11= newl; Line(frac_line_11) = {p6, p40};
frac_line_12= newl; Line(frac_line_12) = {p7, p40};
frac_line_13= newl; Line(frac_line_13) = {p8, p9};
frac_line_14= newl; Line(frac_line_14) = {p8, p11};
frac_line_15= newl; Line(frac_line_15) = {p9, p10};
frac_line_16= newl; Line(frac_line_16) = {p10, p11};
frac_line_17= newl; Line(frac_line_17) = {p12, p13};
frac_line_18= newl; Line(frac_line_18) = {p12, p15};
frac_line_19= newl; Line(frac_line_19) = {p13, p14};
frac_line_20= newl; Line(frac_line_20) = {p14, p15};
frac_line_21= newl; Line(frac_line_21) = {p16, p42};
frac_line_22= newl; Line(frac_line_22) = {p16, p51};
frac_line_23= newl; Line(frac_line_23) = {p17, p42};
frac_line_24= newl; Line(frac_line_24) = {p17, p50};
frac_line_25= newl; Line(frac_line_25) = {p18, p19};
frac_line_26= newl; Line(frac_line_26) = {p18, p50};
frac_line_27= newl; Line(frac_line_27) = {p19, p51};
frac_line_28= newl; Line(frac_line_28) = {p20, p44};
frac_line_29= newl; Line(frac_line_29) = {p20, p51};
frac_line_30= newl; Line(frac_line_30) = {p21, p44};
frac_line_31= newl; Line(frac_line_31) = {p21, p50};
frac_line_32= newl; Line(frac_line_32) = {p22, p23};
frac_line_33= newl; Line(frac_line_33) = {p22, p50};
frac_line_34= newl; Line(frac_line_34) = {p23, p51};
frac_line_35= newl; Line(frac_line_35) = {p24, p27};
frac_line_36= newl; Line(frac_line_36) = {p24, p46};
frac_line_37= newl; Line(frac_line_37) = {p25, p26};
frac_line_38= newl; Line(frac_line_38) = {p25, p46};
frac_line_39= newl; Line(frac_line_39) = {p26, p27};
frac_line_40= newl; Line(frac_line_40) = {p28, p31};
frac_line_41= newl; Line(frac_line_41) = {p28, p48};
frac_line_42= newl; Line(frac_line_42) = {p29, p30};
frac_line_43= newl; Line(frac_line_43) = {p29, p48};
frac_line_44= newl; Line(frac_line_44) = {p30, p31};
frac_line_45= newl; Line(frac_line_45) = {p32, pin01};
frac_line_45_1= newl; Line(frac_line_45_1) = {pin01, pin00};
frac_line_45_2= newl; Line(frac_line_45_2) = {pin00, p33};
frac_line_46= newl; Line(frac_line_46) = {p32, p35};
frac_line_47= newl; Line(frac_line_47) = {p32, p36};
frac_line_48= newl; Line(frac_line_48) = {p33, p34};
frac_line_49= newl; Line(frac_line_49) = {p33, p37};
frac_line_50= newl; Line(frac_line_50) = {p34, pout00};
frac_line_50_1= newl; Line(frac_line_50_1) = {pout00, pout01};
frac_line_50_2= newl; Line(frac_line_50_2) = {pout01, p35};
frac_line_51= newl; Line(frac_line_51) = {p34, p38};
frac_line_52= newl; Line(frac_line_52) = {p35, p39};
frac_line_53= newl; Line(frac_line_53) = {p36, pin11};
frac_line_53_1= newl; Line(frac_line_53_1) = {pin11, pin10};
frac_line_53_2= newl; Line(frac_line_53_2) = {pin10, p37};
frac_line_54= newl; Line(frac_line_54) = {p36, p39};
frac_line_55= newl; Line(frac_line_55) = {p37, p38};
frac_line_56= newl; Line(frac_line_56) = {p38, pout10};
frac_line_56_1= newl; Line(frac_line_56_1) = {pout10, pout11};
frac_line_56_2= newl; Line(frac_line_56_2) = {pout11, p39};
frac_line_57= newl; Line(frac_line_57) = {p40, p41};
frac_line_58= newl; Line(frac_line_58) = {p42, p43};
frac_line_59= newl; Line(frac_line_59) = {p43, p45};
frac_line_60= newl; Line(frac_line_60) = {p43, p47};
frac_line_61= newl; Line(frac_line_61) = {p44, p45};
frac_line_62= newl; Line(frac_line_62) = {p46, p47};
frac_line_63= newl; Line(frac_line_63) = {p47, p49};
frac_line_64= newl; Line(frac_line_64) = {p48, p49};
frac_line_65= newl; Line(frac_line_65) = {p50, p51};

in_line_low = newl; Line(in_line_low) = {pin00, pin10};

in_line_hi = newl; Line(in_line_hi) = {pin01, pin11};

out_line_low = newl; Line(out_line_low) = {pout00, pout10};

out_line_hi = newl; Line(out_line_hi) = {pout01, pout11};


f0_l[]={1,2,3,4,5,6,7,8,14,69};
f1_l[]={9,10,11,12,13};
f2_l[]={14,15,16,17};
f3_l[]={18,19,20,21};
f4_l[]={22,23,24,25,26,27,28};
f5_l[]={29,30,31,32,33,35};
f6_l[]={36,37,38,39,40};
f7_l[]={41,42,43,44,45};

Transfinite Line {f0_l[]} = nf_0;
Transfinite Line {f1_l[]} = nf_1;
Transfinite Line {f2_l[]} = nf_2;
Transfinite Line {f3_l[]} = nf_3;
Transfinite Line {f4_l[]} = nf_4;
Transfinite Line {f5_l[]} = nf_5;
Transfinite Line {f6_l[]} = nf_6;
Transfinite Line {f7_l[]} = nf_7;

// End of line specification

//Physical bodies 
Curve Loop(1) = {9, 11, 12, -13, -10};
Plane Surface(1) = {1};

Line{66} In Surface{1};

Curve Loop(2) = {2, -4, 3, -14, -1};
Plane Surface(2) = {2};
Curve Loop(3) = {5, -14, -7, 8, -68, 69, 72, -6};
Plane Surface(3) = {3};

Line{71} In Surface{3};
Line{73} In Surface{3};
Line{70} In Surface{3};
Line{67} In Surface{3};

Curve Loop(4) = {15, -17, -16, -14};
Plane Surface(4) = {4};

Curve Loop(5) = {20, 21, -19, 18};
Plane Surface(5) = {5};

Curve Loop(6) = {44, -42, 41, -45, -43};
Plane Surface(6) = {6};

Line{73} In Surface{6};

Curve Loop(7) = {38, 40, -36, 37, -39};
Plane Surface(7) = {7};

Line{71} In Surface{7};

Curve Loop(8) = {25, -27, 26, 28, -23, 22, -24};
Plane Surface(8) = {8};

Line{67} In Surface{8};

Curve Loop(9) = {33, 35, -30, 29, -31, 32, -34};
Plane Surface(9) = {9};

Line{70} In Surface{9};

Curve Loop(10) = {76, 59, -75, -47};
Plane Surface(10) = {10};


Curve Loop(11) = {65, -57, -55, 78};
Plane Surface(11) = {11};
Curve Loop(12) = {77, -63, -56, 53};
Plane Surface(12) = {12};


Curve Loop(13) = {61, -65, -64, -63, -62, -60, -59, -58};
Plane Surface(13) = {13};
Curve Loop(14) = {64, -78, -54, 77};
Plane Surface(14) = {14};
Curve Loop(15) = {49, -55, -54, -53, -51, -48, -47, -46};
Plane Surface(15) = {15};
Curve Loop(16) = {58, -76, -46, 50};
Plane Surface(16) = {16};
Curve Loop(17) = {51, 56, -62, -52};
Plane Surface(17) = {17};
Curve Loop(18) = {49, 57, -61, -50};
Plane Surface(18) = {18};
Curve Loop(19) = {48, 52, -60, -75};
Plane Surface(19) = {19};

Line{74} In Surface{8};
Line{74} In Surface{9};
Line{14} In Surface{4};
Line{14} In Surface{3};
Line{14} In Surface{2};
Line{66} In Surface{2};


//Volume body
Surface Loop(1) = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
Volume(1) = {1};
Surface {1,2,3,4,5,6,7,8,9} In Volume{1};
Coherence;


Physical Volume("RockMatrix_1") = {1};
Physical Surface("BCInlet") = {10};
Physical Surface("BCOutlet") = {11, 12};
Physical Surface("BCImpervious") = {13, 14, 15, 16, 17, 18, 19};
Physical Surface("FRACTURE_0") = {2,3};
Physical Surface("FRACTURE_1") = {1};
Physical Surface("FRACTURE_2") = {4};
Physical Surface("FRACTURE_3") = {5};
Physical Surface("FRACTURE_4") = {8};
Physical Surface("FRACTURE_5") = {9};
Physical Surface("FRACTURE_6") = {7};
Physical Surface("FRACTURE_7") = {6};


Physical Curve("FracturesIntersections") = {14, 67, 70, 71, 73, 66, 74};
