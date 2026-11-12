3-RPR:
Macaluay2:
i1 : 
-- ==========================================
-- 3-RPR Groebner Basis Comparison
-- Variables: x > y > c > s
-- Coefficient field: Rational numbers
-- ==========================================


-- 1. LEXICOGRAPHIC ORDER

R1 = QQ[x1,y1,c1,s1, MonomialOrder => Lex];

i2 : 
I1 = ideal(
    x1^2 + y1^2 + 4*x1 - 2*x1*c1 - 4*c1 - 2*y1*s1 + 3,
    x1^2 + y1^2 + 2*x1*c1 - 4*x1 - 4*c1 + 2*y1*s1 + 3,
    x1^2 + y1^2 - 2*x1*s1 + 2*y1*c1 - 6*y1 - 6*c1 + 9,
    c1^2 + s1^2 - 1
);


o2 :  Ideal of R1 
i3 : 
print "========== LEX ==========";

========== LEX ========== 
i4 : 
elapsedTime (G1 = gb I1);

 -- .00390847s elapsed
i5 : 
print "Basis:";

Basis: 
i6 : gens G1;


o6 :  Matrix R1 
1
 ⟵R1 
4
  
i7 : 
print "Number of basis polynomials:";

Number of basis polynomials: 
i8 : numColumns (gens G1);

i9 : 

-- 2. GRADED LEXICOGRAPHIC ORDER

R2 = QQ[x2,y2,c2,s2, MonomialOrder => GLex];

i10 : 
I2 = ideal(
    x2^2 + y2^2 + 4*x2 - 2*x2*c2 - 4*c2 - 2*y2*s2 + 3,
    x2^2 + y2^2 + 2*x2*c2 - 4*x2 - 4*c2 + 2*y2*s2 + 3,
    x2^2 + y2^2 - 2*x2*s2 + 2*y2*c2 - 6*y2 - 6*c2 + 9,
    c2^2 + s2^2 - 1
);


o10 :  Ideal of R2 
i11 : 
print "========== GRLEX ==========";

========== GRLEX ========== 
i12 : 
elapsedTime (G2 = gb I2);

 -- .000582751s elapsed
i13 : 
print "Basis:";

Basis: 
i14 : gens G2;


o14 :  Matrix R2 
1
 ⟵R2 
9
  
i15 : 
print "Number of basis polynomials:";

Number of basis polynomials: 
i16 : numColumns (gens G2);

i17 : 

-- 3. GRADED REVERSE LEXICOGRAPHIC ORDER

R3 = QQ[x3,y3,c3,s3, MonomialOrder => GRevLex];

i18 : 
I3 = ideal(
    x3^2 + y3^2 + 4*x3 - 2*x3*c3 - 4*c3 - 2*y3*s3 + 3,
    x3^2 + y3^2 + 2*x3*c3 - 4*x3 - 4*c3 + 2*y3*s3 + 3,
    x3^2 + y3^2 - 2*x3*s3 + 2*y3*c3 - 6*y3 - 6*c3 + 9,
    c3^2 + s3^2 - 1
);


o18 :  Ideal of R3 
i19 : 
print "========== GREVLEX ==========";

========== GREVLEX ========== 
i20 : 
elapsedTime (G3 = gb I3);

 -- .000638023s elapsed
 
Python:
 [Running] python -u "c:\Users\ayush\Desktop\BTP\Planer_3RPR_Parallel_Robot.py"
=================================================================
1. ALGORITHMIC STEP METRICS (LEX ORDERING)
=================================================================
 � Total S-Polynomial Pairs Evaluated : 120
 � Non-Zero Reductions Added          : 12
 � Pairs Reduced to Zero              : 108
 � Total Division Steps               : 2094
 � Terms Moved to Remainder           : 69
 � New Generators (terms, degree)     : [(3, 2), (5, 2), (10, 2), (6, 2), (5, 2), (4, 3), (11, 5), (5, 5), (9, 5), (4, 5), (4, 4), (3, 6)]
 � Final Intermediate Basis Size      : 16

[Done] exited with code=0 in 8.575 seconds

[Running] python -u "c:\Users\ayush\Desktop\BTP\Planer_3RPR_Parallel_Robot.py"
=================================================================
1. ALGORITHMIC STEP METRICS (GRLEX ORDERING)
=================================================================
 � Total S-Polynomial Pairs Evaluated : 78
 � Non-Zero Reductions Added          : 9
 � Pairs Reduced to Zero              : 69
 � Total Division Steps               : 1338
 � Terms Moved to Remainder           : 44
 � New Generators (terms, degree)     : [(3, 2), (5, 2), (10, 2), (4, 2), (5, 2), (4, 2), (4, 3), (5, 3), (4, 2)]
 � Final Intermediate Basis Size      : 13

[Done] exited with code=0 in 4.37 seconds

[Running] python -u "c:\Users\ayush\Desktop\BTP\Planer_3RPR_Parallel_Robot.py"
=================================================================
1. ALGORITHMIC STEP METRICS (GREVLEX ORDERING)
=================================================================
 � Total S-Polynomial Pairs Evaluated : 66
 � Non-Zero Reductions Added          : 8
 � Pairs Reduced to Zero              : 58
 � Total Division Steps               : 1105
 � Terms Moved to Remainder           : 38
 � New Generators (terms, degree)     : [(3, 2), (5, 2), (9, 2), (4, 2), (4, 2), (5, 3), (4, 2), (4, 2)]
 � Final Intermediate Basis Size      : 12

[Done] exited with code=0 in 3.826 seconds


solution:
Real solution 1:
  x = 0.000000, y = 1.000000
  c = 1.000000, s = 0.000000  ->  theta = 0.00°
  

2R1P:

macaluay2:
========== LEX ========== 
i1 : 
-- ============================================
-- 2R1P: Groebner Basis Comparison
-- Variables: c1 > s1 > c2 > s2 > p
-- Coefficient field: QQ
-- ============================================


-- 1. LEXICOGRAPHIC ORDER

R1 = QQ[c1L,s1L,c2L,s2L,pL, MonomialOrder => Lex];

i2 : 
I1 = ideal(
    c1L^2 + s1L^2 - 1,
    c2L^2 + s2L^2 - 1,
    c1L*c2L - s1L*s2L - 1,
    c1L*s2L + c2L*s1L,
    2*c1L + pL - 1/2,
    2*s1L - 39/20
);


o2 :  Ideal of R1 
i3 : 
print "========== LEX ==========";

========== LEX ========== 
i4 : 
elapsedTime (G1 = gb I1);

 -- .00982479s elapsed
i5 : 
print "Groebner basis:";

Groebner basis: 
i6 : gens G1;


o6 :  Matrix R1 
1
 ⟵R1 
5
  
i7 : 
print "Number of basis polynomials:";

Number of basis polynomials: 
i8 : numColumns (gens G1);

i9 : 

-- 2. GRADED LEXICOGRAPHIC ORDER

R2 = QQ[c1G,s1G,c2G,s2G,pG, MonomialOrder => GLex];

i10 : 
I2 = ideal(
    c1G^2 + s1G^2 - 1,
    c2G^2 + s2G^2 - 1,
    c1G*c2G - s1G*s2G - 1,
    c1G*s2G + c2G*s1G,
    2*c1G + pG - 1/2,
    2*s1G - 39/20
);


o10 :  Ideal of R2 
i11 : 
print "========== GRLEX ==========";

========== GRLEX ========== 
i12 : 
elapsedTime (G2 = gb I2);

 -- .000499041s elapsed
i13 : 
print "Groebner basis:";

Groebner basis: 
i14 : gens G2;


o14 :  Matrix R2 
1
 ⟵R2 
5
  
i15 : 
print "Number of basis polynomials:";

Number of basis polynomials: 
i16 : numColumns (gens G2);

i17 : 

-- 3. GRADED REVERSE LEXICOGRAPHIC ORDER

R3 = QQ[c1R,s1R,c2R,s2R,pR, MonomialOrder => GRevLex];

i18 : 
I3 = ideal(
    c1R^2 + s1R^2 - 1,
    c2R^2 + s2R^2 - 1,
    c1R*c2R - s1R*s2R - 1,
    c1R*s2R + c2R*s1R,
    2*c1R + pR - 1/2,
    2*s1R - 39/20
);


o18 :  Ideal of R3 
i19 : 
print "========== GREVLEX ==========";

========== GREVLEX ========== 
i20 : 
elapsedTime (G3 = gb I3);

 -- .000452246s elapsed


python:
[Running] python -u "c:\Users\ayush\Desktop\BTP\Planar_2R1P_Robot_Diagram.py"
=================================================================
1. ALGORITHMIC STEP METRICS (LEX ORDERING)
=================================================================
 � Total S-Polynomial Pairs Evaluated : 36
 � Non-Zero Reductions Added          : 3
 � Pairs Reduced to Zero              : 33
 � Total Division Steps               : 166
 � Terms Moved to Remainder           : 8
 � New Generators (terms, degree)     : [(3, 1), (2, 1), (3, 2)]

[Done] exited with code=0 in 1.025 seconds

[Running] python -u "c:\Users\ayush\Desktop\BTP\Planar_2R1P_Robot_Diagram.py"
=================================================================
1. ALGORITHMIC STEP METRICS (GRLEX ORDERING)
=================================================================
 � Total S-Polynomial Pairs Evaluated : 36
 � Non-Zero Reductions Added          : 3
 � Pairs Reduced to Zero              : 33
 � Total Division Steps               : 166
 � Terms Moved to Remainder           : 8
 � New Generators (terms, degree)     : [(3, 1), (2, 1), (3, 2)]

[Done] exited with code=0 in 0.931 seconds

[Running] python -u "c:\Users\ayush\Desktop\BTP\Planar_2R1P_Robot_Diagram.py"
=================================================================
1. ALGORITHMIC STEP METRICS (GREVLEX ORDERING)
=================================================================
 � Total S-Polynomial Pairs Evaluated : 36
 � Non-Zero Reductions Added          : 3
 � Pairs Reduced to Zero              : 33
 � Total Division Steps               : 157
 � Terms Moved to Remainder           : 8
 � New Generators (terms, degree)     : [(3, 1), (3, 2), (2, 1)]

[Done] exited with code=0 in 1.015 seconds


solution:
Solution 1: p = 0.055590
 -> VALID: Satisfies physical limits (0.0 <= p <= 1.0)
 c1 = 0.2222, s1 = 0.9750, c2 = 0.2222, s2 = -0.9750
 • Revolute Joint 1 (theta1) = 77.16°
 • Revolute Joint 2 (theta2) = -77.16°
 • Prismatic Extension (p) = 0.0556 units

Solution 2: p = 0.944410
 -> VALID: Satisfies physical limits (0.0 <= p <= 1.0)
 c1 = -0.2222, s1 = 0.9750, c2 = -0.2222, s2 = -0.9750
 • Revolute Joint 1 (theta1) = 102.84°
 • Revolute Joint 2 (theta2) = -102.84°
 • Prismatic Extension (p) = 0.9444 units
 
 
 octahedral_stewart_platform:
 macaluay2:
 
 i1 : 
-- =====================================================
-- Groebner Basis Comparison
-- 3D Polynomial System
-- Variables: x1 > y1 > z1 > x2 > y2 > z2 > x3 > y3 > z3
-- Coefficient field: QQ
-- =====================================================


-- =====================================================
-- 1. LEXICOGRAPHIC ORDER
-- =====================================================

R1 = QQ[x1L,y1L,z1L,x2L,y2L,z2L,x3L,y3L,z3L,
        MonomialOrder => Lex];

i2 : 
I1 = ideal(

    x1L^2 + y1L^2 + z1L^2 - 21,

    x1L^2 - 12*x1L + y1L^2 + z1L^2 - 9,

    x2L^2 - 12*x2L + y2L^2 + z2L^2 - 69/5,

    x2L^2 - 4*x2L + y2L^2 - 10*y2L + z2L^2 + 11,

    x3L^2 - 4*x3L + y3L^2 - 10*y3L + z3L^2 - 271/13,

    x3L^2 + y3L^2 + z3L^2 - 3801/65,

    (x1L-x2L)^2 + (y1L-y2L)^2 + (z1L-z2L)^2 - 16,

    (x2L-x3L)^2 + (y2L-y3L)^2 + (z2L-z3L)^2 - 18,

    (x1L-x3L)^2 + (y1L-y3L)^2 + (z1L-z3L)^2 - 10
);


o2 :  Ideal of R1 
i3 : 
print "========== LEX ==========";

========== LEX ========== 
i4 : 
elapsedTime (G1 = gb I1);

 -- .0346577s elapsed
i5 : 
print "Groebner basis:";

Groebner basis: 
i6 : gens G1;


o6 :  Matrix R1 
1
 ⟵R1 
9
  
i7 : 
print "Number of basis polynomials:";

Number of basis polynomials: 
i8 : numColumns (gens G1);

i9 : 

-- =====================================================
-- 2. GRADED LEXICOGRAPHIC ORDER
-- =====================================================

R2 = QQ[x1G,y1G,z1G,x2G,y2G,z2G,x3G,y3G,z3G,
        MonomialOrder => GLex];

i10 : 
I2 = ideal(

    x1G^2 + y1G^2 + z1G^2 - 21,

    x1G^2 - 12*x1G + y1G^2 + z1G^2 - 9,

    x2G^2 - 12*x2G + y2G^2 + z2G^2 - 69/5,

    x2G^2 - 4*x2G + y2G^2 - 10*y2G + z2G^2 + 11,

    x3G^2 - 4*x3G + y3G^2 - 10*y3G + z3G^2 - 271/13,

    x3G^2 + y3G^2 + z3G^2 - 3801/65,

    (x1G-x2G)^2 + (y1G-y2G)^2 + (z1G-z2G)^2 - 16,

    (x2G-x3G)^2 + (y2G-y3G)^2 + (z2G-z3G)^2 - 18,

    (x1G-x3G)^2 + (y1G-y3G)^2 + (z1G-z3G)^2 - 10
);


o10 :  Ideal of R2 
i11 : 
print "========== GRLEX ==========";

========== GRLEX ========== 
i12 : 
elapsedTime (G2 = gb I2);

 -- .0166764s elapsed
i13 : 
print "Groebner basis:";

Groebner basis: 
i14 : gens G2;


o14 :  Matrix R2 
1
 ⟵R2 
27
  
i15 : 
print "Number of basis polynomials:";

Number of basis polynomials: 
i16 : numColumns (gens G2);

i17 : 

-- =====================================================
-- 3. GRADED REVERSE LEXICOGRAPHIC ORDER
-- =====================================================

R3 = QQ[x1R,y1R,z1R,x2R,y2R,z2R,x3R,y3R,z3R,
        MonomialOrder => GRevLex];

i18 : 
I3 = ideal(

    x1R^2 + y1R^2 + z1R^2 - 21,

    x1R^2 - 12*x1R + y1R^2 + z1R^2 - 9,

    x2R^2 - 12*x2R + y2R^2 + z2R^2 - 69/5,

    x2R^2 - 4*x2R + y2R^2 - 10*y2R + z2R^2 + 11,

    x3R^2 - 4*x3R + y3R^2 - 10*y3R + z3R^2 - 271/13,

    x3R^2 + y3R^2 + z3R^2 - 3801/65,

    (x1R-x2R)^2 + (y1R-y2R)^2 + (z1R-z2R)^2 - 16,

    (x2R-x3R)^2 + (y2R-y3R)^2 + (z2R-z3R)^2 - 18,

    (x1R-x3R)^2 + (y1R-y3R)^2 + (z1R-z3R)^2 - 10
);


o18 :  Ideal of R3 
i19 : 
print "========== GREVLEX ==========";

========== GREVLEX ========== 
i20 : 
elapsedTime (G3 = gb I3);

 -- .00961632s elapsed
i21 : 
print "Groebner basis:";

Groebner basis: 
i22 : gens G3;


o22 :  Matrix R3 
1
 ⟵R3 
26
  
i23 : 
print "Number of basis polynomials:";

Number of basis polynomials: 
i24 : numColumns (gens G3);

python:
lex - more than 10 minutes
[Running] python -u "c:\Users\ayush\Desktop\BTP\temp.py"
=================================================================
1. ALGORITHMIC STEP METRICS (GRLEX ORDERING)
=================================================================
 � S-Polynomial Pairs Evaluated       : 891
 � Coprime Pairs Skipped (1st Crit.)  : 1062
 � Non-Zero Reductions Added          : 54
 � Pairs Reduced to Zero (Total)      : 1899
 � Total Division Steps               : 22970
 � Terms Moved to Remainder           : 455
 � New Generators (terms, degree)     : [(2, 1), (4, 2), (5, 2), (3, 1), (7, 2), (3, 1), (3, 2), (4, 2), (4, 2), (5, 2), (6, 3), (7, 3), (8, 3), (6, 3), (6, 3), (7, 3), (7, 3), (7, 3), (10, 4), (11, 4), (8, 4), (10, 4), (11, 4), (8, 4), (11, 4), (8, 4), (10, 3), (13, 4), (9, 4), (9, 2), (11, 3), (11, 3), (11, 3), (11, 3), (13, 3), (13, 3), (9, 4), (8, 2), (9, 3), (9, 3), (9, 3), (12, 3), (11, 3), (11, 3), (9, 2), (9, 2), (9, 3), (9, 3), (11, 3), (11, 3), (9, 3), (10, 2), (9, 2), (9, 3)]

[Done] exited with code=0 in 95.099 seconds

[Running] python -u "c:\Users\ayush\Desktop\BTP\temp.py"
=================================================================
1. ALGORITHMIC STEP METRICS (GREVLEX ORDERING)
=================================================================
 � S-Polynomial Pairs Evaluated       : 859
 � Coprime Pairs Skipped (1st Crit.)  : 1032
 � Non-Zero Reductions Added          : 53
 � Pairs Reduced to Zero (Total)      : 1838
 � Total Division Steps               : 23614
 � Terms Moved to Remainder           : 444
 � New Generators (terms, degree)     : [(2, 1), (4, 2), (5, 2), (3, 1), (7, 2), (3, 1), (3, 2), (4, 2), (4, 2), (5, 2), (6, 3), (7, 3), (8, 3), (6, 3), (6, 3), (7, 3), (7, 3), (7, 3), (11, 4), (11, 4), (8, 4), (11, 4), (10, 4), (8, 4), (10, 4), (8, 4), (10, 3), (13, 4), (9, 4), (9, 2), (11, 3), (11, 3), (11, 3), (11, 3), (13, 3), (13, 3), (9, 4), (8, 2), (9, 3), (9, 3), (9, 3), (10, 3), (11, 3), (11, 3), (9, 2), (9, 2), (9, 3), (9, 3), (11, 3), (11, 3), (9, 3), (10, 2), (9, 2)]

[Done] exited with code=0 in 92.248 seconds

solution:
Solution 1:
x1 = 1.0000000000
y1 = 2.0000000000
z1 = 4.0000000000
x2 = 3.4000000000
y2 = 5.2000000000
z2 = 4.0000000000
x3 = 0.6769230769
y3 = 3.4923076923
z3 = 6.7692307692
Maximum equation residual = 7.105e-15

Solution 2:
x1 = 1.0000000000
y1 = 2.0000000000
z1 = -4.0000000000
x2 = 3.4000000000
y2 = 5.2000000000
z2 = -4.0000000000
x3 = 0.6769230769
y3 = 3.4923076923
z3 = -6.7692307692
Maximum equation residual = 1.066e-14

Solution 3:
x1 = 1.0000000000
y1 = 2.4043175490
z1 = 3.7708430256
x2 = 3.6489603312
y2 = 5.3991682650
z2 = 3.8886494471
x3 = 1.5321103224
y3 = 3.1502327941
z3 = 6.7974697042
Maximum equation residual = 5.329e-15

Solution 4:
x1 = 1.0000000000
y1 = 2.4043175490
z1 = -3.7708430256
x2 = 3.6489603312
y2 = 5.3991682650
z2 = -3.8886494471
x3 = 1.5321103224
y3 = 3.1502327941
z3 = -6.7974697042
Maximum equation residual = 7.105e-15