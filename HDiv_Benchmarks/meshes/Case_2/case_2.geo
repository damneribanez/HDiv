

// Gmsh project created on Sun Mar 31 19:13:00 2019
SetFactory("OpenCASCADE");

//h = 1.0; // 923 tet

h = 0.14; // 3907 tet

//h = 0.062; // 38742 tet

//Small fractures defintion 
p1 = newp; Point(p1) = {0.5, 0.75, 0.625,h};
p2 = newp; Point(p2) = {0.625, 0.75, 0.625,h};
p3 = newp; Point(p3) = {0.75, 0.75, 0.625,h};
p4 = newp; Point(p4) = {0.5, 0.625, 0.625,h};
p5 = newp; Point(p5) = {0.75, 0.625, 0.625,h};
p6 = newp; Point(p6) = {0.5, 0.5, 0.625,h};
p7 = newp; Point(p7) = {0.625, 0.5, 0.625,h};
p8 = newp; Point(p8) = {0.75, 0.5, 0.625,h};

p9 = newp; Point(p9) = {0.625, 0.75, 0.75,h};
p10 = newp; Point(p10) = {0.625, 0.75, 0.5,h};
p11 = newp; Point(p11) = {0.625, 0.625, 0.75,h};
p12 = newp; Point(p12) = {0.625, 0.625, 0.5,h};
p13 = newp; Point(p13) = {0.625, 0.5, 0.75,h};
p14 = newp; Point(p14) = {0.625, 0.5, 0.5,h};

p15 = newp; Point(p15) = {0.5, 0.625, 0.75,h};
p16 = newp; Point(p16) = {0.5, 0.625, 0.5,h};
p17 = newp; Point(p17) = {0.625, 0.625, 0.75,h};
p18 = newp; Point(p18) = {0.625, 0.625, 0.5,h};
p19 = newp; Point(p19) = {0.75, 0.625, 0.75,h};
p20 = newp; Point(p20) = {0.75, 0.625, 0.5,h};

//Center point small fractures
p21 = newp; Point(p21) = {0.625, 0.625, 0.625,h};

l1 = newl; Curve(l1) = {p1,p2};
l2 = newl; Line(l2) = {p2,p3};
l3 = newl; Line(l3) = {p4,p21};
l4 = newl; Line(l4) = {p21,p5};
l5 = newl; Line(l5) = {p6,p7};
l6 = newl; Line(l6) = {p7,p8};

l7 = newl; Line(l7) = {p9,p2};
l8 = newl; Line(l8) = {p2,p10};
l9 = newl; Line(l9) = {p11,p21};
l10 = newl; Line(l10) = {p21,p12};
l11 = newl; Line(l11) = {p7,p13};
l12 = newl; Line(l12) = {p7,p14};

l13 = newl; Line(l13) = {p4,p15};
l14 = newl; Line(l14) = {p4,p16};
l15 = newl; Line(l15) = {p5,p19};
l16 = newl; Line(l16) = {p5,p20};

l17 = newl; Line(l17) = {p1,p4};
l18 = newl; Line(l18) = {p4,p6};
l17 = newl; Line(l17) = {p2,p21};
l18 = newl; Line(l18) = {p21,p7};
l19 = newl; Line(l19) = {p3,p5};
l20 = newl; Line(l20) = {p5,p8};

l21 = newl; Line(l21) = {p9,p11};
l22 = newl; Line(l22) = {p11,p13};
l23 = newl; Line(l23) = {p10,p12};
l24 = newl; Line(l24) = {p12,p14};

l25 = newl; Line(l25) = {p11,p15};
l26 = newl; Line(l26) = {p11,p19};
l27 = newl; Line(l27) = {p12,p16};
l28 = newl; Line(l28) = {p12,p20};

//Small fractures surfaces 

frac_1_1 = newll; Line Loop(frac_1_1) = {1, 19, -3, -17};
frac_1_2 = newll; Line Loop(frac_1_2) = {2, 21, -4, -19};
frac_1_3 = newll; Line Loop(frac_1_3) = {3, 20, -5, -18};
frac_1_4 = newll; Line Loop(frac_1_4) = {20, 6, -22, -4};

su_frac_1_1 = news; Plane Surface(su_frac_1_1) = {frac_1_1};
su_frac_1_2 = news; Plane Surface(su_frac_1_2) = {frac_1_2};
su_frac_1_3 = news; Plane Surface(su_frac_1_3) = {frac_1_3};
su_frac_1_4 = news; Plane Surface(su_frac_1_4) = {frac_1_4};

frac_2_1 = newll; Line Loop(frac_2_1) = {7, 19, -9, -23};
frac_2_2 = newll; Line Loop(frac_2_2) = {8, 25, -10, -19};
frac_2_3 = newll; Line Loop(frac_2_3) = {9, 20, 11, -24};
frac_2_4 = newll; Line Loop(frac_2_4) = {10, 26, -12, -20};

su_frac_2_1 = news; Plane Surface(su_frac_2_1) = {frac_2_1};
su_frac_2_2 = news; Plane Surface(su_frac_2_2) = {frac_2_2};
su_frac_2_3 = news; Plane Surface(su_frac_2_3) = {frac_2_3};
su_frac_2_4 = news; Plane Surface(su_frac_2_4) = {frac_2_4};

frac_3_1 = newll; Line Loop(frac_3_1) = {13, -27, 9, -3};
frac_3_2 = newll; Line Loop(frac_3_2) = {14, -29, -10, -3};
frac_3_3 = newll; Line Loop(frac_3_3) = {9, 4, 15, -28};
frac_3_4 = newll; Line Loop(frac_3_4) = {10, 30, -16, -4};

su_frac_3_1 = news; Plane Surface(su_frac_3_1) = {frac_3_1};
su_frac_3_2 = news; Plane Surface(su_frac_3_2) = {frac_3_2};
su_frac_3_3 = news; Plane Surface(su_frac_3_3) = {frac_3_3};
su_frac_3_4 = news; Plane Surface(su_frac_3_4) = {frac_3_4};

//Medium fractures definition

p22 = newp; Point(p22) = {1.0, 1.0, 0.75,h};
p23 = newp; Point(p23) = {0.75, 1.0, 0.75,h};
p24 = newp; Point(p24) = {0.5, 1.0, 0.75,h};
p25 = newp; Point(p25) = {1.0, 0.75, 0.75,h};
p26 = newp; Point(p26) = {0.5, 0.75, 0.75,h};
p27 = newp; Point(p27) = {1.0, 0.5, 0.75,h};
p28 = newp; Point(p28) = {0.75, 0.5, 0.75,h};
p29 = newp; Point(p29) = {0.5, 0.5, 0.75,h};

p30 = newp; Point(p30) = {0.75, 1.0, 0.5,h};
p31 = newp; Point(p31) = {0.75, 1.0, 1.0,h};
p32 = newp; Point(p32) = {0.75, 0.75, 0.5,h};
p33 = newp; Point(p33) = {0.75, 0.75, 1.0,h};
p34 = newp; Point(p34) = {0.75, 0.5, 0.5,h};
p35 = newp; Point(p35) = {0.75, 0.5, 1.0,h};

p36 = newp; Point(p36) = {1.0, 0.75, 0.5,h};
p37 = newp; Point(p37) = {1.0, 0.75, 1.0,h};
p38 = newp; Point(p38) = {0.5, 0.75, 0.5,h};
p39 = newp; Point(p39) = {0.5, 0.75, 1.0,h};

//Center point medium fractures
p40 = newp; Point(p40) = {0.75, 0.75, 0.75,h};

p41 = newp; Point(p41) = {0.75, 1., 0.875,h};
p42 = newp; Point(p42) = {0.75, 0.875, 1.,h};
p43 = newp; Point(p43) = {0.875, 0.75, 1.,h};
p44 = newp; Point(p44) = {1., 0.875, 0.75,h};
p45 = newp; Point(p45) = {1., 0.75, 0.875,h};
p46 = newp; Point(p46) = {0.875, 1., 0.75,h};

l55 = newl; Line(l55) = {p23,p24};
l56 = newl; Line(l56) = {p23,p46};
l57 = newl; Line(l57) = {p46,p22};
l58 = newl; Line(l58) = {p9,p26};
l59 = newl; Line(l59) = {p9,p40};
l60 = newl; Line(l60) = {p25,p40};

l61 = newl; Line(l61) = {p13,p29};
l62 = newl; Line(l62) = {p13,p28};
l63 = newl; Line(l63) = {p27, p28};

l64 = newl; Line(l64) = {p23,p30};
l65 = newl; Line(l65) = {p33,p40};
l66 = newl; Line(l66) = {p3, p40};
l67 = newl; Line(l67) = {p3,p32};
l68 = newl; Line(l68) = {p28,p35};
l69 = newl; Line(l69) = {p8,p28};
l70 = newl; Line(l70) = {p8,p34};
l71 = newl; Line(l71) = {p23,p41};
l72 = newl; Line(l72) = {p41,p31};

l73 = newl; Line(l73) = {p26,p39};
l74 = newl; Line(l74) = {p1,p26};
l75 = newl; Line(l75) = {p1, p38};
l76 = newl; Line(l76) = {p45,p37};
l77 = newl; Line(l77) = {p25,p45};
l78 = newl; Line(l78) = {p25,p36};

l79 = newl; Line(l79) = {p24,p26};
l80 = newl; Line(l80) = {p15,p26};
l81 = newl; Line(l81) = {p15,p29};
l82 = newl; Line(l82) = {p23,p40};
l83 = newl; Line(l83) = {p40, p19};
l84 = newl; Line(l84) = {p19, p28};
l85 = newl; Line(l85) = {p27,p25};
l86 = newl; Line(l86) = {p25, p44};
l87 = newl; Line(l87) = {p44, p22};

l88 = newl; Line(l88) = {p30,p32};
l89 = newl; Line(l89) = {p32,p20};
l90 = newl; Line(l90) = {p31,p42};
l91 = newl; Line(l91) = {p42,p33};
l92 = newl; Line(l92) = {p33,p35};

l93 = newl; Line(l93) = {p10,p38};
l94 = newl; Line(l94) = {p10,p32};
l95 = newl; Line(l95) = {p32,p36};

l96 = newl; Line(l96) = {p39,p33};
l97 = newl; Line(l97) = {p33,p43};
l98 = newl; Line(l98) = {p43,p37};
l99 = newl; Line(l99) = {p20,p34};





//Medium fractures surfaces 
frac_5_1 = newll; Line Loop(frac_5_1) = {l85, l60, l83, l84, l63};
frac_5_2 = newll; Line Loop(frac_5_2) = {l86, l87, l57, l56, l82, l60};
frac_5_3_1 = newll; Line Loop(frac_5_3_1) = {l84, l62, 24, 28};
frac_5_3_2 = newll; Line Loop(frac_5_3_2) = {l83, 28, 23, l59};
frac_5_3_3 = newll; Line Loop(frac_5_3_3) = {24, l61, l81, 27};
frac_5_3_4 = newll; Line Loop(frac_5_3_4) = {23, 27, l80, l58};
frac_5_4 = newll; Line Loop(frac_5_4) = {l82, l59, l58, l79, l55};

su_frac_5_1 = news; Plane Surface(su_frac_5_1) = {frac_5_1};
su_frac_5_2 = news; Plane Surface(su_frac_5_2) = {frac_5_2};
su_frac_5_3_1 = news; Plane Surface(su_frac_5_3_1) = {frac_5_3_1};
su_frac_5_3_2 = news; Plane Surface(su_frac_5_3_2) = {frac_5_3_2};
su_frac_5_3_3 = news; Plane Surface(su_frac_5_3_3) = {frac_5_3_3};
su_frac_5_3_4 = news; Plane Surface(su_frac_5_3_4) = {frac_5_3_4};
su_frac_5_4 = news; Plane Surface(su_frac_5_4) = {frac_5_4};


frac_6_1 = newll; Line Loop(frac_6_1) = {l97, l98, l76, l77, l60, l65};
frac_6_2 = newll; Line Loop(frac_6_2) = {l96, l65, l59, l58, l73};
frac_6_3 = newll; Line Loop(frac_6_3) = {l60, l78, l95, l67, l66};
frac_6_4_1 = newll; Line Loop(frac_6_4_1) = {l58, l74, 1, 7};
frac_6_4_2 = newll; Line Loop(frac_6_4_2) = {l59, 7, 2, l66};
frac_6_4_3 = newll; Line Loop(frac_6_4_3) = {1, l75, l93, 8};
frac_6_4_4 = newll; Line Loop(frac_6_4_4) = {2, 8, l94, l67};


su_frac_6_1 = news; Plane Surface(su_frac_6_1) = {frac_6_1};
su_frac_6_2 = news; Plane Surface(su_frac_6_2) = {frac_6_2};
su_frac_6_3 = news; Plane Surface(su_frac_6_3) = {frac_6_3};su_frac_6_4_1 = news; Plane Surface(su_frac_6_4_1) = {frac_6_4_1};
su_frac_6_4_2 = news; Plane Surface(su_frac_6_4_2) = {frac_6_4_2};
su_frac_6_4_3 = news; Plane Surface(su_frac_6_4_3) = {frac_6_4_3};
su_frac_6_4_4 = news; Plane Surface(su_frac_6_4_4) = {frac_6_4_4};

frac_7_1 = newll; Line Loop(frac_7_1) = {l91, l90, l72, l71, l82, l65};
frac_7_2_1 = newll; Line Loop(frac_7_2_1) = {l83, l66, 21, 15};
frac_7_2_2 = newll; Line Loop(frac_7_2_2) = {15, l84, l69, 22};
frac_7_2_3 = newll; Line Loop(frac_7_2_3) = {l67, 21, 16, l89};
frac_7_2_4 = newll; Line Loop(frac_7_2_4) = {16, 22, l70, l99};
frac_7_3 = newll; Line Loop(frac_7_3) = {l92, l65, l83, l84, l68};
frac_7_4 = newll; Line Loop(frac_7_4) = {l82, l64, l88, l67, l66};

su_frac_7_1 = news; Plane Surface(su_frac_7_1) = {frac_7_1};
su_frac_7_2_1 = news; Plane Surface(su_frac_7_2_1) =
{frac_7_2_1};
su_frac_7_2_2 = news; Plane Surface(su_frac_7_2_2) =
{frac_7_2_2};
su_frac_7_2_3 = news; Plane Surface(su_frac_7_2_3) =
{frac_7_2_3};
su_frac_7_2_4 = news; Plane Surface(su_frac_7_2_4) =
{frac_7_2_4};
su_frac_7_3 = news; Plane Surface(su_frac_7_3) = {frac_7_3};
su_frac_7_4 = news; Plane Surface(su_frac_7_4) = {frac_7_4};
//Big fractures defintion

p47 = newp; Point(p47) = {1.0, 1.0, 0.5,h};
p48 = newp; Point(p48) = {0.5, 1.0, 0.5,h};
p49 = newp; Point(p49) = {0.0, 1.0, 0.5,h};
p50 = newp; Point(p50) = {1.0, 0.5, 0.5,h};
p51 = newp; Point(p51) = {0.0, 0.5, 0.5,h};
p52 = newp; Point(p52) = {1.0, 0.0, 0.5,h};
p53 = newp; Point(p53) = {0.5, 0.0, 0.5,h};
p54 = newp; Point(p54) = {0.0, 0.0, 0.5,h};

p55 = newp; Point(p55) = {0.5, 1.0, 0.0,h};
p56 = newp; Point(p56) = {0.5, 1.0, 1.0,h};
p57 = newp; Point(p57) = {0.5, 0.5, 0.0,h};
p58 = newp; Point(p58) = {0.5, 0.5, 1.0,h};
p59 = newp; Point(p59) = {0.5, 0.0, 0.0,h};
p60 = newp; Point(p60) = {0.5, 0.0, 1.0,h};

p61 = newp; Point(p61) = {0.5, 1.0, 0.0,h};
p62 = newp; Point(p62) = {0.5, 1.0, 1.0,h};
p63 = newp; Point(p63) = {0.5, 0.5, 0.0,h};
p64 = newp; Point(p64) = {0.5, 0.5, 1.0,h};
p65 = newp; Point(p65)= {0.5, 0.0, 0.0,h};
p66 = newp; Point(p66) = {0.5, 0.0, 1.0,h};

p67 = newp; Point(p67) = {1.0, 0.5, 0.0,h};
p68 = newp; Point(p68) = {1.0, 0.5, 1.0,h};
p69 = newp; Point(p69) = {0.0, 0.5, 0.0,h};
p70 = newp; Point(p70) = {0.0, 0.5, 1.0,h};

//Center point big fractures
p71 = newp; Point(p71) = {0.5, 0.5, 0.5,h};

p600 = newp; Point(p600) = {0.25, 0.5, 0.,h};
p601 = newp; Point(p601) = {0.5, 0.25, 0.,h};
p602 = newp; Point(p602) = {0.25, 0., 0.5,h};
//p603 = newp; Point(p603) = {0.25, 0., 0.5,h};
p604 = newp; Point(p604) = {0.5, 0., 0.25,h};
p605 = newp; Point(p605) = {0., 0.25, 0.5,h};
p606 = newp; Point(p606) = {0., 0.5, 0.25,h};

l142 = newl; Line(l142) = {p47,p30};
l143 = newl; Line(l143) = {p30,p48};
l144 = newl; Line(l144) = {p48,p49};
l145 = newl; Line(l145) = {p50,p34};
l146 = newl; Line(l146) = {p34,p14};
l147 = newl; Line(l147) = {p14,p71};
l148 = newl; Line(l148) = {p71,p51};
l149 = newl; Line(l149) = {p52,p53};
l150_1 = newl; Line(l150_1) = {p53,p602};
l150_2 = newl; Line(l150_2) = {p602,p54};

l151 = newl; Line(l151) = {p47,p36};
l152 = newl; Line(l152) = {p36,p50};
l153 = newl; Line(l153) = {p50,p52};
l154 = newl; Line(l154) = {p48,p38};
l155 = newl; Line(l155) = {p38,p16};
l156 = newl; Line(l156) = {p16,p71};
l157 = newl; Line(l157) = {p71,p53};
l158 = newl; Line(l158) = {p49,p51};
l159_1 = newl; Line(l159_1) = {p51,p605};
l159_2 = newl; Line(l159_2) = {p605,p54};

l160 = newl; Line(l160) = {p70,p58};
l161 = newl; Line(l161) = {p58,p35};
l162 = newl; Line(l162) = {p35,p68};
l163_1 = newl; Line(l163_1) = {p69,p600};
l163_2 = newl; Line(l163_2) = {p600,p57};
l164 = newl; Line(l164) = {p57,p67};

l165 = newl; Line(l165) = {p60,p58};
l166 = newl; Line(l166) = {p58,p39};
l167 = newl; Line(l167) = {p39,p56};
l168_1 = newl; Line(l168_1) = {p59,p601};
l168_2 = newl; Line(l168_2) = {p601,p57};
l169 = newl; Line(l169) = {p57,p55};

l170 = newl; Line(l170) = {p55,p48};
l171 = newl; Line(l171) = {p48,p24};
l172 = newl; Line(l172) = {p24,p56};
l173 = newl; Line(l173) = {p57,p71};
l174 = newl; Line(l174) = {p71,p6};
l175 = newl; Line(l175) = {p6,p29};
l176 = newl; Line(l176) = {p29,p58};
l177_1 = newl; Line(l177_1) = {p59,p604};
l177_2 = newl; Line(l177_2) = {p604,p53};
l178 = newl; Line(l178) = {p53,p60};

l179 = newl; Line(l179) = {p67,p50};
l180 = newl; Line(l180) = {p50,p27};
l181 = newl; Line(l181) = {p27,p68};
l182_1 = newl; Line(l182_1) = {p69,p606};
l182_2 = newl; Line(l182_2) = {p606,p51};
l183 = newl; Line(l183) = {p51,p70};

frac_8_1 = newll; Line Loop(frac_8_1) = {l160, l176, l175, l174, l148, l183};
frac_8_2_1 = newll; Line Loop(frac_8_2_1) = {l161, l68, 62, l61, l176};
frac_8_2_2 = newll; Line Loop(frac_8_2_2) = {l162, l181, l63, l68};
frac_8_2_3_1 = newll; Line Loop(frac_8_2_3_1) = {l61, 11, 5, l175};
frac_8_2_3_2 = newll; Line Loop(frac_8_2_3_2) = {l62, l69, 6, 11};
frac_8_2_3_3 = newll; Line Loop(frac_8_2_3_3) = {5, 12, l147, l174};
frac_8_2_3_4 = newll; Line Loop(frac_8_2_3_4) = {6, l70, l146, 12};
frac_8_2_4 = newll; Line Loop(frac_8_2_4) = {l63, l180, l145, l70, l69};
frac_8_3 = newll; Line Loop(frac_8_3) = {l145, l146, l147, l173, l164, l179};
frac_8_4 = newll; Line Loop(frac_8_4) = {l148, l182_2, l182_1, l163_1, l163_2, l173};

su_frac_8_1 = news; Plane Surface(su_frac_8_1) = {frac_8_1};
su_frac_8_2_1 = news; Plane Surface(su_frac_8_2_1) = {frac_8_2_1};
su_frac_8_2_2 = news; Plane Surface(su_frac_8_2_2) = {frac_8_2_2};
su_frac_8_2_3_1 = news; Plane Surface(su_frac_8_2_3_1) =
{frac_8_2_3_1};
su_frac_8_2_3_2 = news; Plane Surface(su_frac_8_2_3_2) =
{frac_8_2_3_2};
su_frac_8_2_3_3 = news; Plane Surface(su_frac_8_2_3_3) =
{frac_8_2_3_3};
su_frac_8_2_3_4 = news; Plane Surface(su_frac_8_2_3_4) =
{frac_8_2_3_4};
su_frac_8_2_4 = news; Plane Surface(su_frac_8_2_4) = {frac_8_2_4};
su_frac_8_3 = news; Plane Surface(su_frac_8_3) = {frac_8_3};
su_frac_8_4 = news; Plane Surface(su_frac_8_4) = {frac_8_4};

frac_9_1 = newll; Line Loop(frac_9_1) = {l170, l154, l155, l156, l173, l169};
frac_9_2_1 = newll; Line Loop(frac_9_2_1) = {l171, l79, l74, l75, l154};
frac_9_2_2 = newll; Line Loop(frac_9_2_2) = {l172, l167, l73, l79};
frac_9_2_3 = newll; Line Loop(frac_9_2_3) = {l73, l166, l176, l81, l80};
frac_9_2_4_1 = newll; Line Loop(frac_9_2_4_1) = {l75, 17, 14, l155};
frac_9_2_4_2 = newll; Line Loop(frac_9_2_4_2) = {l74, l80, 13, 17};
frac_9_2_4_3 = newll; Line Loop(frac_9_2_4_3) = {14, 18, l174, l156};
frac_9_2_4_4 = newll; Line Loop(frac_9_2_4_4) = {13, l81, l175, 18};
frac_9_3 = newll; Line Loop(frac_9_3) = {l173, l157, l177_2, l177_1, l168_1, l168_2};
frac_9_4 = newll; Line Loop(frac_9_4) = {l174, l175, l176, l165, l178, l157};

su_frac_9_1 = news; Plane Surface(su_frac_9_1) = {frac_9_1};
su_frac_9_2_1 = news; Plane Surface(su_frac_9_2_1) =
{frac_9_2_1};
su_frac_9_2_2 = news; Plane Surface(su_frac_9_2_2) =
{frac_9_2_2};
su_frac_9_2_3 = news; Plane Surface(su_frac_9_2_3) =
{frac_9_2_3};
su_frac_9_2_4_1 = news; Plane Surface(su_frac_9_2_4_1) =
{frac_9_2_4_1};
su_frac_9_2_4_2 = news; Plane Surface(su_frac_9_2_4_2) =
{frac_9_2_4_2};
su_frac_9_2_4_3 = news; Plane Surface(su_frac_9_2_4_3) =
{frac_9_2_4_3};
su_frac_9_2_4_4 = news; Plane Surface(su_frac_9_2_4_4) =
{frac_9_2_4_4};
su_frac_9_3 = news; Plane Surface(su_frac_9_3) = {frac_9_3};
su_frac_9_4 = news; Plane Surface(su_frac_9_4) = {frac_9_4};

frac_10_1 = newll; Line Loop(frac_10_1) = {l153, l145, l146, l147, l157, l149};
frac_10_2_1 = newll; Line Loop(frac_10_2_1) = {l152, l95, l89, l99, l145};
frac_10_2_2 = newll; Line Loop(frac_10_2_2) = {l151, l142, l88, l95};
frac_10_2_3_1 = newll; Line Loop(frac_10_2_3_1) = {l99, 30, 26, l146};
frac_10_2_3_2 = newll; Line Loop(frac_10_2_3_2) = {l89, l94, 25, 30};
frac_10_2_3_3 = newll; Line Loop(frac_10_2_3_3) = {26, 29, l156, l147};
frac_10_2_3_4 = newll; Line Loop(frac_10_2_3_4) = {25, l93, l155, 29};
frac_10_2_4 = newll; Line Loop(frac_10_2_4) = {l88, l143, l154, l93, l94};
frac_10_3 = newll; Line Loop(frac_10_3) = {l157, l148, l159_1, l159_2, l150_1, l150_2};
frac_10_4 = newll; Line Loop(frac_10_4) = {l156, l155, l154, l144, l158, l148};

su_frac_10_1 = news; Plane Surface(su_frac_10_1) = {frac_10_1};
su_frac_10_2_1 = news; Plane Surface(su_frac_10_2_1) = {frac_10_2_1};
su_frac_10_2_2 = news; Plane Surface(su_frac_10_2_2) = {frac_10_2_2};
su_frac_10_2_3_1 = news; Plane Surface(su_frac_10_2_3_1) = {frac_10_2_3_1};
su_frac_10_2_3_2 = news; Plane Surface(su_frac_10_2_3_2) = {frac_10_2_3_2};
su_frac_10_2_3_3 = news; Plane Surface(su_frac_10_2_3_3) = {frac_10_2_3_3};
su_frac_10_2_3_4 = news; Plane Surface(su_frac_10_2_3_4) = {frac_10_2_3_4};
su_frac_10_2_4 = news; Plane Surface(su_frac_10_2_4) = {frac_10_2_4};
su_frac_10_3 = news; Plane Surface(su_frac_10_3) = {frac_10_3};
su_frac_10_4 = news; Plane Surface(su_frac_10_4) = {frac_10_4};

//Volume small block 
Surface Loop(1) = {su_frac_7_2_2, su_frac_8_2_3_2, su_frac_5_3_1, su_frac_2_3, su_frac_1_4, su_frac_3_3};
Volume(1) = {1};

//Volume small block 2
Surface Loop(2) = {su_frac_7_2_4, su_frac_8_2_3_4, su_frac_3_4, su_frac_10_2_3_1, su_frac_2_4, su_frac_1_4};
Volume(2) = {2};

//Volume small block 3
Surface Loop(3) = {su_frac_1_1, su_frac_9_2_4_1, su_frac_3_2, su_frac_6_4_3, su_frac_10_2_3_4, su_frac_2_2};
Volume(3) = {3};

//Volume small block 4
Surface Loop(4) = {su_frac_2_4, su_frac_8_2_3_3, su_frac_1_3, su_frac_9_2_4_3, su_frac_10_2_3_3, su_frac_3_2};
Volume(4) = {4};

//Volume small block 5
Surface Loop(5) = {su_frac_6_4_2, su_frac_1_2, su_frac_3_3, su_frac_2_1, su_frac_5_3_2, su_frac_7_2_1};
Volume(5) = {5};

//Volume small block 6
Surface Loop(6) = {su_frac_7_2_3, su_frac_10_2_3_2, su_frac_1_2, su_frac_2_2, su_frac_6_4_4, su_frac_3_4};
Volume(6) = {6};

//Volume small block 7
Surface Loop(7) = {su_frac_2_1, su_frac_1_1, su_frac_9_2_4_2, su_frac_5_3_4, su_frac_6_4_1, su_frac_3_1};
Volume(7) = {7};

//Volume small block 8
Surface Loop(8) = {su_frac_5_3_3, su_frac_8_2_3_1, su_frac_9_2_4_4, su_frac_3_1, su_frac_2_3, su_frac_1_3};
Volume(8) = {8};


//Inlet definition
p79 = newp; Point(p79) = {0.25, 0., 0.,h};
p80 = newp; Point(p80) = {0., 0.25, 0.,h};
p81 = newp; Point(p81) = {0., 0., 0.25,h};
p82 = newp; Point(p82) = {0.25, 0.25, 0.,h};
p83 = newp; Point(p83) = {0.25, 0., 0.25,h};
p84 = newp; Point(p84) = {0., 0.25, 0.25,h};

//Outlet Definition
p85 = newp; Point(p85) = {0.875, 1., 1.,h};
p86 = newp; Point(p86) = {1., 0.875, 1.,h};
p87 = newp; Point(p87) = {1., 1., 0.875,h};
p88 = newp; Point(p88) = {0.875, 0.875, 1.,h};
p89 = newp; Point(p89) = {0.875, 1., 0.875,h};
p90 = newp; Point(p90) = {1., 0.875, 0.875,h};

//Domain definition
p91 = newp; Point(p91) = {1.0, 0.0, 1.0,h};
p92 = newp; Point(p92) = {1.0, 1.0, 1.0,h};
p93 = newp; Point(p93) = {0.0, 0.0, 1.0,h};
p94 = newp; Point(p94) = {0.0, 1.0, 1.0,h};
p95 = newp; Point(p95) = {1.0, 0.0, 0.0,h};
p96 = newp; Point(p96) = {1.0, 1.0, 0.0,h};
p97 = newp; Point(p97) = {0.0, 0.0, 0.0,h};
p98 = newp; Point(p98) = {0.0, 1.0, 0.0,h};

//Inlet lines definition
l249 = newl; Line(l249) = {p92,p87};
l250 = newl; Line(l250) = {p87,p90};
l251 = newl; Line(l251) = {p86,p90};
l252 = newl; Line(l252) = {p86,p92};
l253 = newl; Line(l253) = {p86,p88};
l254 = newl; Line(l254) = {p88,p85};
l255 = newl; Line(l255) = {p85,p92};
l256 = newl; Line(l256) = {p85,p89};
l257 = newl; Line(l257) = {p89,p87};

inlet_iz = newll; Line Loop(inlet_iz) = {l249, l250, l251, l252};
inlet_bo = newll; Line Loop(inlet_bo) = {l249, l255, l256, l257};
inlet_fr = newll; Line Loop(inlet_fr) = {l252, l253, l254, l255};

su_inlet_iz = news; Plane Surface(su_inlet_iz) = {inlet_iz};
su_inlet_bo = news; Plane Surface(su_inlet_bo) = {inlet_bo};
su_inlet_fr = news; Plane Surface(su_inlet_fr) = {inlet_fr};

//Outlet lines definition
l258 = newl; Line(l258) = {p97,p79};
l259 = newl; Line(l259) = {p79,p82};
l260 = newl; Line(l260) = {p82,p80};
l261 = newl; Line(l261) = {p80,p97};
l262 = newl; Line(l262) = {p80,p84};
l263 = newl; Line(l263) = {p84,p81};
l264 = newl; Line(l264) = {p81,p97};
l265 = newl; Line(l265) = {p81,p83};
l266 = newl; Line(l266) = {p83,p79};

outlet_ba = newll; Line Loop(outlet_ba) = {l258, l259, l260, l261};
outlet_ri = newll; Line Loop(outlet_ri) = {l261, l262, l263, l264};
outlet_to = newll; Line Loop(outlet_to) = {l258, l264, l265, l266};

su_outlet_ba = news; Plane Surface(su_outlet_ba) = {outlet_ba};
su_outlet_ri = news; Plane Surface(su_outlet_ri) = {outlet_ri};
su_outlet_to = news; Plane Surface(su_outlet_to) = {outlet_to};


//Domain lines definition

p99 = newp; Point(p99) = {0., 0.25, 0.5,h};
p100 = newp; Point(p100) = {0., 0.5, 0.25,h};
p101 = newp; Point(p101) = {0.25, 0.5, 0.,h};
p102 = newp; Point(p102) = {0.5, 0., 0.25,h};
p103 = newp; Point(p103) = {0.25, 0., 0.5,h};
p104 = newp; Point(p104) = {0.5, 0.25, 0.0,h};

l280 = newl; Line(l280) = {p96,p55};
l281 = newl; Line(l281) = {p55,p98};
l282 = newl; Line(l282) = {p85,p31};
l283 = newl; Line(l283) = {p31,p56};
l284 = newl; Line(l284) = {p56,p94};

l285 = newl; Line(l285) = {p95,p59};
l286 = newl; Line(l286) = {p59,p79};
l287 = newl; Line(l287) = {p91,p60};
l288 = newl; Line(l288) = {p60,p93};

l289 = newl; Line(l289) = {p93,p70};
l290 = newl; Line(l290) = {p70,p94};
l291 = newl; Line(l291) = {p80,p69};
l292 = newl; Line(l292) = {p69,p98};
l293 = newl; Line(l293) = {p91,p68};
l294 = newl; Line(l294) = {p68,p37};
l295 = newl; Line(l295) = {p37,p86};
l296 = newl; Line(l296) = {p95,p67};
l297 = newl; Line(l297) = {p67,p96};

l298 = newl; Line(l298) = {p96,p47};
l299 = newl; Line(l299) = {p47,p22};
l300 = newl; Line(l300) = {p22,p87};
l301 = newl; Line(l301) = {p49,p94};
l302 = newl; Line(l302) = {p49,p98};

l303 = newl; Line(l303) = {p93,p54};
l304 = newl; Line(l304) = {p54,p81};
l305 = newl; Line(l305) = {p95,p52};
l306 = newl; Line(l306) = {p52,p91};

l307 = newl; Line(l307) = {p89,p41};
l308 = newl; Line(l308) = {p89,p46};
l309 = newl; Line(l309) = {p90,p45};
l310 = newl; Line(l310) = {p90,p44};
l311 = newl; Line(l311) = {p88,p42};
l312 = newl; Line(l312) = {p88,p43};
l313 = newl; Line(l313) = {p82,p600};
l314 = newl; Line(l314) = {p82,p601};
l315 = newl; Line(l315) = {p83,p604};
l316 = newl; Line(l316) = {p83,p602};
l317 = newl; Line(l317) = {p84,p606};
l318 = newl; Line(l318) = {p84,p605};

//Medium block 1
Curve Loop(279) = {283, -176, -55, 71, 72};
Plane Surface(280) = {279};
Curve Loop(280) = {90, 91, -96, 170, -283};
Plane Surface(281) = {280};
Surface Loop(9) = {281, 135, 280, 222, 113, 122};
Volume(9) = {9};

//Medium block 2
Curve Loop(281) = {282, 90, -311, 255};
Plane Surface(282) = {281};
Curve Loop(282) = {91, 97, -312, 311};
Plane Surface(283) = {282};
Curve Loop(283) = {295, 254, 312, 98};
Plane Surface(284) = {283};
Curve Loop(284) = {76, 295, 252, 309};
Plane Surface(285) = {284};
Curve Loop(285) = {77, -309, 310, -86};
Plane Surface(286) = {285};
Curve Loop(286) = {282, -72, -307, -257};
Plane Surface(287) = {286};
Curve Loop(287) = {258, -300, -57, -308};
Plane Surface(288) = {287};
Curve Loop(288) = {307, -71, 56, -308};
Plane Surface(289) = {288};
Curve Loop(289) = {310, 87, 300, 251};
Plane Surface(290) = {289};
Surface Loop(10) = {264, 262, 263, 288, 290, 286, 121, 283, 284, 285, 282, 287, 289, 108, 135};
Volume(10) = {10};

//Medium block 3
Curve Loop(290) = {64, 143, 175, -55};
Surface(291) = {290};
Surface Loop(11) = {291, 221, 247, 141, 113, 124, 125, 126, 127};
Volume(11) = {11};

//Medium block 4
Curve Loop(292) = {57, -299, 142, -64, 56};
Plane Surface(292) = {292};
Curve Loop(293) = {78, -152, 299, -87, -86};
Plane Surface(293) = {293};
Surface Loop(12) = {293, 123, 242, 292, 108, 141};
Volume(12) = {12};


//Medium block 5
Curve Loop(294) = {186, 294, -76, -77, -85};
Plane Surface(294) = {294};
Curve Loop(295) = {164, 294, -98, -97, 92};
Plane Surface(295) = {295};
Surface Loop(13) = {295, 202, 294, 107, 140, 121};
Volume(13) = {13};

//Medium block 6
Curve Loop(296) = {163, -92, -96, -169};
Plane Surface(296) = {296};
Surface Loop(14) = {296, 201, 223, 140, 122, 111, 112, 110, 109};
Volume(14) = {14};

//Medium block 7
Curve Loop(297) = {185, 85, 78, 153};
Plane Surface(297) = {297};
Surface Loop(15) = {297, 207, 241, 107, 137, 136, 138, 139, 123};
Volume(15) = {15};


//Big block 1
Curve Loop(298) = {290, -284, -170, -169, -162};
Plane Surface(298) = {298};
Curve Loop(299) = {290, -301, 159, 189};
Plane Surface(299) = {299};
Curve Loop(300) = {301, -284, -176, -175, 144};
Plane Surface(300) = {300};
Surface Loop(16) = {299, 298, 300, 249, 200, 222, 223, 221, 225, 227, 224, 226};
Volume(16) = {16};

//Big block 2
Curve Loop(301) = {162, -168, 288, 289};
Plane Surface(301) = {301};
Curve Loop(302) = {303, -161, -160, 189, -289};
Plane Surface(302) = {302};
Curve Loop(303) = {303, -151, -150, 183, 288};
Plane Surface(303) = {303};
Surface Loop(17) = {248, 303, 302, 301, 229, 200};
Volume(17) = {17};

//Big block 3
Curve Loop(304) = {168, 163, 164, -293, 287};
Plane Surface(304) = {304};
Curve Loop(305) = {287, -183, -149, 306};
Plane Surface(305) = {305};
Curve Loop(306) = {306, 293, -186, -185, 154};
Plane Surface(306) = {306};
Surface Loop(18) = {304, 306, 305, 240, 201, 202, 203, 204, 205, 206, 229, 207};
Volume(18) = {18};

//Big block 4
Curve Loop(307) = {305, -154, -184, -296};
Plane Surface(307) = {307};
Curve Loop(308) = {296, -167, -172, -171, -285};
Plane Surface(308) = {308};
Curve Loop(309) = {305, 149, -182, -181, -285};
Plane Surface(309) = {309};
Surface Loop(19) = {308, 307, 309, 228, 208, 240};
Volume(19) = {19};

//Big block 6
Curve Loop(323) = {159, -188, -187, 292, -302};
Plane Surface(322) = {323};
Curve Loop(324) = {302, -281, 174, 144};
Plane Surface(323) = {324};
Curve Loop(325) = {292, -281, -173, -166, -165};
Plane Surface(324) = {325};
Surface Loop(20) = {322, 324, 323, 220, 249, 209};
Volume(20) = {20};

//Big block 7
Curve Loop(326) = {173, -280, -297, -167};
Plane Surface(325) = {326};
Curve Loop(327) = {174, -143, -142, -298, 280};
Plane Surface(326) = {327};
Curve Loop(328) = {297, 298, 152, 153, -184};
Plane Surface(327) = {328};
Surface Loop(21) = {326, 327, 325, 208, 220, 241, 242, 247, 244, 246, 243, 245};
Volume(21) = {21};

//Big block 8
Curve Loop(329) = {151, 304, 271, 265, -286, 181, 182, 150};
Plane Surface(328) = {329};
Curve Loop(330) = {188, 160, -318, 317};
Plane Surface(329) = {330};
Curve Loop(331) = {187, -317, -269, 291};
Plane Surface(330) = {331};
Curve Loop(332) = {318, 161, 304, -270};
Plane Surface(331) = {332};
Curve Loop(333) = {272, 316, 151, 304};
Plane Surface(332) = {333};
Curve Loop(334) = {316, -150, -182, -315};
Plane Surface(333) = {334};
Curve Loop(335) = {273, -286, 181, -315};
Plane Surface(334) = {335};
Curve Loop(336) = {266, 314, -171, 286};
Plane Surface(335) = {336};
Curve Loop(337) = {313, 166, -172, -314};
Plane Surface(336) = {337};
Curve Loop(338) = {165, -313, 267, 291};
Plane Surface(337) = {338};
Surface Loop(22) = {209, 248, 228, 336, 337, 277, 335, 278, 331, 329, 330, 334, 279, 332, 333};
Volume(22) = {22};

Coherence;

Physical Volume("RockMatrix_1") = {16, 17, 20, 21, 22, 29, 14, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3};
Physical Volume("RockMatrix_2") = {1, 2, 13, 15, 18, 19};
Physical Surface("BCInlet") = {277, 278, 279};
Physical Surface("BCOutlet") = {262, 263, 264};
Physical Surface("BCImpervious") = {322, 299, 298, 300, 323, 326, 330, 329, 331, 302, 301, 304, 295, 296, 281, 280, 291, 292, 288, 289, 287, 282, 284, 283, 290, 286, 285, 293, 297, 294, 306, 307, 327, 325, 324, 308, 336, 335, 337, 334, 333, 332, 309, 305, 303};
Physical Surface("Fractures") = {220, 209, 208, 228, 248, 249, 229, 240, 200, 241, 242, 108, 141, 135, 113, 222, 223, 112, 111, 225, 227, 221, 122, 121, 247, 123, 107, 138, 127, 136, 125, 140, 110, 43, 124, 201, 226, 203, 205, 245, 37, 204, 224, 45, 52, 139, 207, 243, 206, 137, 244, 109, 202, 53, 44, 54, 38, 46, 35, 126, 51, 246, 36};
Physical Curve("FracturesIntersections") = {59,60,61,62,63,65,66,67,68,69,70,73,74,75,79,80,81,82,83,84,88,89,93,94,95,99,145,146,147,148,155,156,157,158,177,178,179,180,58,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};
Physical Point("CrossingIntresections") = {6,71,29,14,7,13,34,8,28,16,4,15,12,21,11,20,5,19,38,1,26,10,2,9,32,3,40};






