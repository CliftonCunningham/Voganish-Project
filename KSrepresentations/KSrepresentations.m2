------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- File: KSrepresentations.m2
-- Author: Nicole Kitt
--
-- The functions found in this file are used in the computations of Ev in the KS project.
--
-- The following functions can be applied to arbitrary Vogan varieties:
--
-- getCompare(List,RankConditions)
-- getSupstrataNumb(List, HashTable)
-- getSupstrata(List, HashTable)
-- getDualStrata(List)
-- getSupstrata2(List, HashTable)
--
-- The following functions are specific to the Vogan variety V in the KS project:
-- x0Matrix()
-- x1Matrix()
-- x2Matrix()
-- x3Matrix()
-- coverCR()
-- coverCr()
-- coverCm()
-- coverCpsi()
-- smallCR()
-- smallCr()
-- smallCm()
-- smallCpsi()
-- JacCR()
-- JacCr()
-- JacCm()
-- JacCpsi()
-- subJacCR()
-- subJacCr()
-- subJacCm()
-- subJacCpsi()
-- InfinityJacCr2()
-- InfinityJacCm2()
-- InfinityJacCm3()
-- InfinityJacCm4()
-- InfinityJacCpsi2()
-- InfinityJacCpsi3()
-- InfinityJacCpsi4()
-- InfinityJacCpsi5()
-- InfinityJacCpsi6()
-- InfinityJacCpsi7()
-- InfinityJacCpsi8()
-- InfinityJacCpsi9()
-- InfinityJacCpsi10()
-- InfinityJacCpsi11()
-- InfinityJacCpsi12()
-- InfinityJacCpsi13()
-- InfinityJacCpsi14()
-- InfinityJacCpsi15()
-- InfinityJacCpsi16()
-- subInfinity(Matrix)
------------------------------------------------------------------------------------------------------------------------------------------------------------------------

-- getCompare is a function that compares a list of orbits  of a given Vogan variety to a particualar orbit C. This comparison 
-- tells us which elements in the list have closure which contain/don't contain C, are equal to C, or incomparable. 
--
-- Example:
-- 
-- i1 : getAllStrataR({1,3,1})
--
-- o1 = {1  3  1, 1  3  1, 1  3  1, 1  3  1, 1  3  1}
--        0  0     0   1    1   0    1   1    1   1
--         0         0        0        0        1
--
-- o1 : List
--
-- i2 : C = new RankConditions from o1#3
--
-- o2 = 1  3  1
--       1   1
--         0
--
-- o2 : RankConditions
--
-- i3 : getCompare(o1,C)
--
-- o3 = {>, >, >, ==, <}
--
-- o3 : List
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- getSupstrataNumb is a function that inputs a list of orbits of a Vogan variety and a paricular orbit C. This 
-- function outputs a sequence of numbers which represent (the location of) orbits in the list
-- which have closures containing C. 
--
-- Example: Using previous inputs/outputs in prev eg.
--
-- i4 : getSubstrataNumb(o1,C)
--
-- o4 : 1 : (4)
--
-- o4 : Sequence
--
-- In this example there is only one orbit in the list o1 whose closure contains C. That being, the 5th 
-- element in o1 (we start at zero in lists/sequences). 
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- getSupstrata takes in a list of orbits from a Vogan variety and a particular orbit C. The function outputs orbits
-- (in triangle format) from the list whose closures contain C.
--
-- Example: Contining outputs/inputs from above eg.
--
-- i5 : getSupstrata(o1,C)
--
-- o5 = {1  3  1}
--        1   1
--          1
--
-- o5 : List
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- getDualstrata is a function that inputs a list L of orbits. It outputs a list of the corresponding dual orbits in L.
--
-- Example:
--
-- i1 : getAllStrataR({1,3,1})
--
-- o1 = {1  3  1, 1  3  1, 1  3  1, 1  3  1, 1  3  1}
--        0  0     0   1     1  0    1   1    1   1
--         0         0        0        0        1
--
-- o1 : List
--
-- i2 : getDualstrata(o1)
--
-- o2 = {1  3  1, 1  3  1, 1  3  1, 1  3  1, 1  3  1}
--        1   1    1  0     0   1    1   1     0  0
--          1       0         0        0         0
--
-- o2 : List
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- getSupstrata2 is a function that takes on a list of orbits L and a particular orbit C. The function finds all orbits C' in 
-- the list L such that C < C' AND C < \hat{C'}. Here \hat{C} denotes the transpose of the dual of C.
--
-- Example:
--
-- i1: getAllStrataR({1,3,1})
--
-- o1: {1  3  1, 1  3  1, 1  3  1, 1  3  1, 1  3  1}
--       0  0     0  1     1   0    1   1    1  1
--        0         0        0        0       1
--
-- o1: List
--
-- i2: C = new RankConditions from o1#3
--
-- o2 = 1  3  1
--       0  1
--        0
--
-- o2 : RankConditions
--
-- i3 : getSupstrata2(o1,C)
--
-- o3 = {1  3  1}
--        1   1
--          0
--
-- o3 : List
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The functions x0Matrix(), x1Matrix(), x2Matrix(), x3Matrix() output a fixed matrix. 
-- The tuple (x0Matrix(), x1Matrix(), x2Matrix(), x3Matrix()) corresponds to an arbitrary element (x_4,x_3,x_2,x_1) of V in the KS project.
-- These functions are used in the following: coverCr(), coverCR(), coverCm(), coverCpsi().
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The functions coverCR(), coverCr(), coverCm(), coverCpsi() output a cover for the closures of the orbits C_R,C_r,C_m, and C_\psi, respectively. 
-- These are the main covers used in the KS project. The subspace Eij in the code is equal to the subspace E_{\lambda_i}^j in the paper. 
-- Note that (x0,x1,x2,x3) in the code corresponds to (x_4,x_3,x_2,x_1) in the paper.
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The functions smallCR(), smallCr(), smallCm(), smallCpsi() output information that tells us if the cover for the closures of CR, Cr, Cm, and Cpsi,  respectively, 
-- is small or semi-small or neither. 
-- The output consists of two lists. Note that the first list will always contain the orbit itself. 
-- We describe these lists below: Let C be one of CR,Cr,Cm Cpsi. And \rho:\widetilde{C}\rightarrow\overline{C} be a cover of \overline{C}.
-- The first list of smallC() consists of all suborbits C' \subseteq \overline{C} such that 2\dim(\rho^{-1}(x)) + dim C' = dim\widetilde{C}, where x\in C' is arbitrary.
-- The second list of smallC() consists of all suborbirs C'\subseteq \overline{C} such that 2\dim(\rho^{-1}(x)) + \dim C' > \dim\widetilde{C}, where x\in C' is arbirary.
-- If the first list only contains the orbit itself and the second list is empty, then the cover is small.
-- If the first list contains orbits other than itself and the second list is empty, then the cover is semi-small.
-- If the second list is non-empty, then the cover is neither small nor semi-small.
-- Note that \dim\widetilde{C} = \dim\overline{C} = \dim C for the H-orbits C = CR,Cr,Cm,Cpsi.
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The function JacCR() outputs the Jacobian for the equations that describe the main cover for CR ( i.e., coverCR() ), C_{KS}*, and \tilde{f}_R^{-1}(0). 
-- The function JacCr() outputs the Jacobian for the equations that describe the main cover for Cr ( i.e., coverCr() ), C_{KS}*, and \tilde{f}_r^{-1}(0). 
-- The function JacCm() outputs the Jacobian for the equations that describe the main cover for Cm ( i.e., coverCm() ), C_{KS}*, and \tilde{f}_m^{-1}(0). 
-- The function JacCpsi() outputs the Jacobian for the equations that describe the main cover for Cpsi ( i.e., coverCpsi() ), C_{KS}*, and \tilde{f}_{\psi}^{-1}(0). 
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The function subJacCR() plugs the regular point (x_{\KS},y_{\KS}(t_1,t_2), with (t_1t_2,t_1+t_2)\in Z_{reg}, and fibre of \rho_R above C_{\KS} into JacCR().
-- The function subJacCr() plugs the regular point (x_{\KS},y_{\KS}(t_1,t_2), with (t_1t_2,t_1+t_2)\in Z_{reg}, and fibre of \rho_r above C_{\KS} into JacCr().
-- The function subJacCm() plugs the regular point (x_{\KS},y_{\KS}(t_1,t_2), with (t_1t_2,t_1+t_2)\in Z_{reg}, and fibre of \rho_m above C_{\KS} into JacCm().
-- The function subJacCpsi() plugs the regular point (x_{\KS},y_{\KS}(t_1,t_2), with (t_1t_2,t_1+t_2)\in Z_{reg}, and fibre of \rho_{psi} above C_{\KS} into JacCpsi().
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The function InfinityJacCr2() outputs a hard-coded Jacobian needed to check smoothness at a particular point at infinity. The point being checked is commented 
-- directly above the function in the code.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The functions InfinityJacCm2(), InfinityJacCm3(), InfinityJacCm4(), InfinityJacCpsi2() output hard-coded Jacobians needed to check smoothness at particular
-- points at inifinity. The point being checked is commented directly above the function in the code.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The functions InfinityJacCpsi3(), InfinityJacCpsi4(), InfinityJacCpsi5(), InfinityJacCpsi6(), InfinityJacCpsi7(), InfinityJacCpsi8(), InfinityJacCpsi9(),
-- InfinityJacCpsi10(), InfinityJacCpsi11(), InfinityJacCpsi12(), InfinityJacCpsi13(), InfinityJacCpsi14(), InfinityJacCpsi15(), InfinityJacCpsi16() output hard-coded
-- Jacobians needed to check smoothness at particular points at inifinity. The point being checked is commented directly above the function in the code.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The function subInfinity(Matrix) only takes in matrices over a polynomial ring with variables x_{i,j,k}, y_{i,j,k}. 
-- So, this function should only be applied to the hard-coded Jacobians in this file.
-- This function plugs in our reg pt (x_KS,y_KS(t_1,t_2)), with (t_1t_2,t_1+t_2)\in Z_{\reg}. It then sets all other available variables to 0 (i.e., plugs in point at infinity).
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- Note: That (y3,y2,y1,y0) in the code corresponds to (y_1,y_2,y_3,y_4) in the paper.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

needs "ComputeDuals.m2" 
needs "NetworkFlow.m2"
needs "ComputeCover.m2" 
needs "VanishingCycles.m2"
needs "VoganV.m2"
needs "PSNF.m2"

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

getCompare = method()
getCompare (List,RankConditions) := (L,C) -> (
    for i when i < #L list C ? new RankConditions from L#i
)

getSupstrataNumb = method()
getSupstrataNumb (List, HashTable) := (L,C) -> (
  delete(null,(0..#L-1) / ( i -> if (getCompare(L,C))#i == symbol < then i))
)

getSupstrata = method()
getSupstrata (List, HashTable) := (L,C) -> (
    N := #getSupstrataNumb(L,C);
    ls := getSupstrataNumb(L,C);
    for i when i < N list L#(ls#i)
 )
 
getDualstrata = method()
getDualstrata (List) := (L) -> (
    for i when i < #L list computeDual L#i
)

getSupstrata2 = method()
getSupstrata2 (List, HashTable) := (L,C) -> (
    sup := getSupstrata (L,C);
    dualsup := getDualstrata (sup);
    supdualsup := getSupstrata (dualsup, computeDual C);
    getDualstrata (supdualsup)
)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

x0Matrix = () -> (
  QQ[x_{0,0,0},x_{0,0,1},x_{0,0,2},x_{0,0,3},x_{0,1,0},x_{0,1,1},x_{0,1,2},x_{0,1,3}];
  matrix{{x_{0,0,0},x_{0,0,1},x_{0,0,2},x_{0,0,3}},{x_{0,1,0},x_{0,1,1},x_{0,1,2},x_{0,1,3}}}
)
  
x1Matrix = () -> (
  QQ[x_{1,0,0},x_{1,0,1},x_{1,0,2},x_{1,0,3},x_{1,1,0},x_{1,1,1},x_{1,1,2},x_{1,1,3},x_{1,2,0},x_{1,2,1},x_{1,2,2},x_{1,2,3},x_{1,3,0},x_{1,3,1},x_{1,3,2},x_{1,3,3}];
  matrix{{x_{1,0,0},x_{1,0,1},x_{1,0,2},x_{1,0,3}},{x_{1,1,0},x_{1,1,1},x_{1,1,2},x_{1,1,3}},{x_{1,2,0},x_{1,2,1},x_{1,2,2},x_{1,2,3}},{x_{1,3,0},x_{1,3,1},x_{1,3,2},x_{1,3,3}}}
)

x2Matrix = () -> (
  use QQ[x_{2,0,0},x_{2,0,1},x_{2,0,2},x_{2,0,3},x_{2,1,0},x_{2,1,1},x_{2,1,2},x_{2,1,3},x_{2,2,0},x_{2,2,1},x_{2,2,2},x_{2,2,3},x_{2,3,0},x_{2,3,1},x_{2,3,2},x_{2,3,3}];
  matrix{{x_{2,0,0},x_{2,0,1},x_{2,0,2},x_{2,0,3}},{x_{2,1,0},x_{2,1,1},x_{2,1,2},x_{2,1,3}},{x_{2,2,0},x_{2,2,1},x_{2,2,2},x_{2,2,3}},{x_{2,3,0},x_{2,3,1},x_{2,3,2},x_{2,3,3}}}
)

x3Matrix = () -> (
  QQ[x_{3,0,0},x_{3,0,1},x_{3,1,0},x_{3,1,1},x_{3,2,0},x_{3,2,1},x_{3,3,0},x_{3,3,1}];
  matrix{{x_{3,0,0},x_{3,0,1}},{x_{3,1,0},x_{3,1,1}},{x_{3,2,0},x_{3,2,1}},{x_{3,3,0},x_{3,3,1}}}
)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

coverCR = () -> (
	CR = new RankConditions from ({2,4,4,4,2},{{2,2,4,2},{0,2,2},{0,0},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4];
	bigRing = ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E22 := matrix{{b1,b2},{b3,b4},{1,0},{0,1}};
	E22 = sub(E22, bigRing);
	E32 := matrix{{1,0},{0,1},{c1,c2},{c3,c4}};
	E32 = sub(E32, bigRing);
	x0 := sub(x0Matrix(), bigRing);
	x1 := sub(x1Matrix(), bigRing);
	x2 := sub(x2Matrix(), bigRing);
	x3 := sub(x3Matrix(), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2(E12) \subset E22
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1,b2},{M56_(1,0),b3,b4},{M56_(2,0),1,0},{M56_(3,0),0,1}};
	eq26 := matrix{{M56_(0,1),b1,b2},{M56_(1,1),b3,b4},{M56_(2,1),1,0},{M56_(3,1),0,1}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq25))|(flatten entries gens minors(3,eq26));
-- Flag eqn x0(E32)=0
	M3 := x0*E32;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x1 \subset E32
	eq31 := matrix{{x_{1,0,0},1,0},{x_{1,1,0},0,1},{x_{1,2,0},c1,c2},{x_{1,3,0},c3,c4}};
	eq32 := matrix{{x_{1,0,1},1,0},{x_{1,1,1},0,1},{x_{1,2,1},c1,c2},{x_{1,3,1},c3,c4}};
	eq33 := matrix{{x_{1,0,2},1,0},{x_{1,1,2},0,1},{x_{1,2,2},c1,c2},{x_{1,3,2},c3,c4}};
	eq34 := matrix{{x_{1,0,3},1,0},{x_{1,1,3},0,1},{x_{1,2,3},c1,c2},{x_{1,3,3},c3,c4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq31))|(flatten entries gens minors(3,eq32));
	equationsCover = equationsCover|(flatten entries gens minors(3,eq33))|(flatten entries gens minors(3,eq34));
-- Flag eqn x1(E22)=0
	M5 := x1*E22;
	equationsCover = equationsCover|(flatten entries M5);
-- NO flag containtment equations
	equationsCover = ideal equationsCover
)

coverCr = () -> (
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4];
	bigRing = ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{1},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,1},{c1,c2,c3}};
	E23 = sub(E23, bigRing);
	E32 := matrix{{1,0},{0,1},{d1,d2},{d3,d4}};
	E32 = sub(E32, bigRing);
	x0 := sub(x0Matrix(), bigRing);
	x1 := sub(x1Matrix(), bigRing);
	x2 := sub(x2Matrix(), bigRing);
	x3 := sub(x3Matrix(), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),b2},{M56_(2,0),1},{M56_(3,0),b3}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),b2},{M56_(2,1),1},{M56_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x2 \subset E23
	eq1  := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,1},{x_{2,3,0},c1,c2,c3}};
	eq2  := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,1},{x_{2,3,1},c1,c2,c3}};
	eq3  := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,1},{x_{2,3,2},c1,c2,c3}};
	eq4  := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,1},{x_{2,3,3},c1,c2,c3}};
	equationsCover = append(equationsCover, det eq1);
	equationsCover = append(equationsCover, det eq2);
	equationsCover = append(equationsCover, det eq3);
	equationsCover = append(equationsCover, det eq4);
-- Flag eqn x0(E32)=0
	M3 := x0*E32;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x1(E21)=0
	M5 := x1*E21;
	equationsCover = equationsCover|(flatten entries M5);
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq66 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),d1,d2},{M6_(3,0),d3,d4}};
	eq67 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),d1,d2},{M6_(3,1),d3,d4}};
	eq68 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),d1,d2},{M6_(3,2),d3,d4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq66))|(flatten entries gens minors(3,eq67))|(flatten entries gens minors(3,eq68));
-- Flag eqn x1 \subset E32
	eq76 := matrix{{x_{1,0,0},1,0},{x_{1,1,0},0,1},{x_{1,2,0},d1,d2},{x_{1,3,0},d3,d4}};
	eq77 := matrix{{x_{1,0,1},1,0},{x_{1,1,1},0,1},{x_{1,2,1},d1,d2},{x_{1,3,1},d3,d4}};
	eq78 := matrix{{x_{1,0,2},1,0},{x_{1,1,2},0,1},{x_{1,2,2},d1,d2},{x_{1,3,2},d3,d4}};
	eq79 := matrix{{x_{1,0,3},1,0},{x_{1,1,3},0,1},{x_{1,2,3},d1,d2},{x_{1,3,3},d3,d4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq76))|(flatten entries gens minors(3,eq77));
	equationsCover = equationsCover|(flatten entries gens minors(3,eq78))|(flatten entries gens minors(3,eq79));
-- flag containtment equations
	con := matrix{{b1,1,0,0},{b2,0,1,0},{1,0,0,1},{b3,c1,c2,c3}};
	equationsCover = append(equationsCover, det con);
	ideal equationsCover	
)

coverCm = () -> (
	Cm := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{0,0},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4,e1,e2,e3,g];
	bigRing = ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{1},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,1},{c1,c2,c3}};
	E23 = sub(E23, bigRing);
	E32 := matrix{{1,0},{0,1},{d1,d2},{d3,d4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,1},{e1,e2,e3}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{1},{g}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(), bigRing);
	x1 := sub(x1Matrix(), bigRing);
	x2 := sub(x2Matrix(), bigRing);
	x3 := sub(x3Matrix(), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2(E12) \subset E21
	M1 := x2*E12;
	eq13 := matrix{{M1_(0,0),b1},{M1_(1,0),b2},{M1_(2,0),1},{M1_(3,0),b3}};
	eq14 := matrix{{M1_(0,1),b1},{M1_(1,1),b2},{M1_(2,1),1},{M1_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq13))|(flatten entries gens minors(2,eq14));
-- Flag eqn x2 \subset E23
	eq31 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,1},{x_{2,3,0},c1,c2,c3}};
	eq32 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,1},{x_{2,3,1},c1,c2,c3}};
	eq33 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,1},{x_{2,3,2},c1,c2,c3}};
	eq34 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,1},{x_{2,3,3},c1,c2,c3}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) = 0
	M2 := x1*E21;
	equationsCover = equationsCover|(flatten entries M2);
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq66 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),d1,d2},{M6_(3,0),d3,d4}};
	eq67 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),d1,d2},{M6_(3,1),d3,d4}};
	eq68 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),d1,d2},{M6_(3,2),d3,d4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq66))|(flatten entries gens minors(3,eq67))|(flatten entries gens minors(3,eq68));
-- Flag eqn x1 \subset E33
	eq1  := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,1},{x_{1,3,0},e1,e2,e3}};
	eq2  := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,1},{x_{1,3,1},e1,e2,e3}};
	eq3  := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,1},{x_{1,3,2},e1,e2,e3}};
	eq4  := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,1},{x_{1,3,3},e1,e2,e3}};
	equationsCover = append(equationsCover, det eq1);
	equationsCover = append(equationsCover, det eq2);
	equationsCover = append(equationsCover, det eq3);
	equationsCover = append(equationsCover, det eq4);
-- Flag eqn x0(E32) = 0
	M4 := x0*E32;
	equationsCover = equationsCover|(flatten entries M4);
-- Flag eqn x0(E33) \subset E41
	M5 := x0*E33;
	eq5  := matrix{{M5_(0,0),1},{M5_(1,0),g}};
	eq6  := matrix{{M5_(0,1),1},{M5_(1,1),g}};
	eq7  := matrix{{M5_(0,2),1},{M5_(1,2),g}};
	equationsCover = append(equationsCover, det eq5);
	equationsCover = append(equationsCover, det eq6);
	equationsCover = append(equationsCover, det eq7);
-- Flag containtment equations
	M9 := matrix{{b1,1,0,0},{b2,0,1,0},{1,0,0,1},{b3,c1,c2,c3}};
	equationsCover = append(equationsCover, det M9);
	M7 := matrix{{1,1,0,0},{0,0,1,0},{d1,0,0,1},{d3,e1,e2,e3}};
	M8 := matrix{{0,1,0,0},{1,0,1,0},{d2,0,0,1},{d4,e1,e2,e3}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = append(equationsCover, det M8);
	ideal equationsCover
	
)

coverCpsi = () -> (
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{1},{b1},{b2},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,1},{c1,c2,c3}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{1},{d1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,1},{f1,f2,f3}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{1},{g}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(), bigRing);
	x1 := sub(x1Matrix(), bigRing);
	x2 := sub(x2Matrix(), bigRing);
	x3 := sub(x3Matrix(), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,1},{x_{2,3,0},c1,c2,c3}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,1},{x_{2,3,1},c1,c2,c3}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,1},{x_{2,3,2},c1,c2,c3}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,1},{x_{2,3,3},c1,c2,c3}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),1},{M56_(1,0),b1},{M56_(2,0),b2},{M56_(3,0),b3}};
	eq26 := matrix{{M56_(0,1),1},{M56_(1,1),b1},{M56_(2,1),b2},{M56_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),1},{M12_(1,0),g}};
	eq42 := matrix{{M12_(0,1),1},{M12_(1,1),g}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),1},{M456_(1,0),g}};
	eq45 := matrix{{M456_(0,1),1},{M456_(1,1),g}};
	eq46 := matrix{{M456_(0,2),1},{M456_(1,2),g}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,1},{x_{1,3,0},f1,f2,f3}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,1},{x_{1,3,1},f1,f2,f3}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,1},{x_{1,3,2},f1,f2,f3}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,1},{x_{1,3,3},f1,f2,f3}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),1},{M5_(1,0),d1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{1,1,0,0},{b1,0,1,0},{b2,0,0,1},{b3,c1,c2,c3}};
	M8 := matrix{{1,1,0},{d1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,1},{e3,f1,f2,f3}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,1},{e4,f1,f2,f3}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	ideal equationsCover

)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

smallCR = () -> (
CR := new RankConditions from ({2,4,4,4,2},{{2,2,4,2},{0,2,2},{0,0},{0}});
semiList := {};
notList := {};
orbitDim := dim(getEquations CR);
CoverCR := coverCR();
for S in getAllSubstrata(CR) do (
	rep := getOtherMatrixRep S;
	A := QQ[a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4];
	re := gens A;
	for M in rep do re = join(re,flatten(entries(flatten(transpose M))));
	fibre := sub(gens CoverCR,matrix{re});
	fibre = sub(fibre,A);
	subOrbitDim := dim(getEquations S);
	fibreDim := dim(ideal fibre);
	if 2*fibreDim + subOrbitDim == orbitDim then semiList = append(semiList,S);
	if 2*fibreDim + subOrbitDim > orbitDim then notList = append(notList,S);
	);
	(semiList,notList)
)

smallCr = () -> (
Cr := new RankConditions from ({2,4,4,4,2},{{2,2,3,2},{0,2,1},{0,0},{0}});
semiList := {};
notList := {};
orbitDim := dim(getEquations Cr);
CoverCr := coverCr();
for S in getAllSubstrata(Cr) do (
	rep := getOtherMatrixRep S;
	A := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4];
	re := gens A;
	for M in rep do re = join(re,flatten(entries(flatten(transpose M))));
	fibre := sub(gens CoverCr,matrix{re});
	fibre = sub(fibre,A);
	subOrbitDim := dim(getEquations S);
	fibreDim := dim(ideal fibre);
	if 2*fibreDim + subOrbitDim == orbitDim then semiList = append(semiList,S);
	if 2*fibreDim + subOrbitDim > orbitDim then notList = append(notList,S);
	);
	(semiList,notList)
)

smallCm = () -> (
Cm := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{0,0},{0}});
semiList := {};
notList := {};
orbitDim := dim(getEquations Cm);
CoverCm := coverCm();
for S in getAllSubstrata(Cm) do (
	rep := getOtherMatrixRep S;
	A := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4,e1,e2,e3,g];
	re := gens A;
	for M in rep do re = join(re,flatten(entries(flatten(transpose M))));
	fibre := sub(gens CoverCm,matrix{re});
	fibre = sub(fibre,A);
	subOrbitDim := dim(getEquations S);
	fibreDim := dim(ideal fibre);
	if 2*fibreDim + subOrbitDim == orbitDim then semiList = append(semiList,S);
	if 2*fibreDim + subOrbitDim > orbitDim then notList = append(notList,S);
	);
	(semiList,notList)
)
 
smallCpsi = () -> (
Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
semiList := {};
notList := {};
orbitDim := dim(getEquations Cpsi);
CoverCpsi := coverCpsi();
for S in getAllSubstrata(Cpsi) do (
	rep := getOtherMatrixRep S;
	A := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	re := gens A;
	for M in rep do re = join(re,flatten(entries(flatten(transpose M))));
	fibre := sub(gens CoverCpsi,matrix{re});
	fibre = sub(fibre,A);
	subOrbitDim := dim(getEquations S);
	fibreDim := dim(ideal fibre);
	if 2*fibreDim + subOrbitDim == orbitDim then semiList = append(semiList,S);
	if 2*fibreDim + subOrbitDim > orbitDim then notList = append(notList,S);
	);
	(semiList,notList)
)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- Ev computation for CR
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

JacCR = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	CR = new RankConditions from ({2,4,4,4,2},{{2,2,4,2},{0,2,2},{0,0},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E22 := matrix{{b1,b2},{b3,b4},{1,0},{0,1}};
	E22 = sub(E22, bigRing);
	E32 := matrix{{1,0},{0,1},{c1,c2},{c3,c4}};
	E32 = sub(E32, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2(E12) \subset E22
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1,b2},{M56_(1,0),b3,b4},{M56_(2,0),1,0},{M56_(3,0),0,1}};
	eq26 := matrix{{M56_(0,1),b1,b2},{M56_(1,1),b3,b4},{M56_(2,1),1,0},{M56_(3,1),0,1}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq25))|(flatten entries gens minors(3,eq26));
-- Flag eqn x0(E32)=0
	M3 := x0*E32;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x1 \subset E32
	eq31 := matrix{{x_{1,0,0},1,0},{x_{1,1,0},0,1},{x_{1,2,0},c1,c2},{x_{1,3,0},c3,c4}};
	eq32 := matrix{{x_{1,0,1},1,0},{x_{1,1,1},0,1},{x_{1,2,1},c1,c2},{x_{1,3,1},c3,c4}};
	eq33 := matrix{{x_{1,0,2},1,0},{x_{1,1,2},0,1},{x_{1,2,2},c1,c2},{x_{1,3,2},c3,c4}};
	eq34 := matrix{{x_{1,0,3},1,0},{x_{1,1,3},0,1},{x_{1,2,3},c1,c2},{x_{1,3,3},c3,c4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq31))|(flatten entries gens minors(3,eq32));
	equationsCover = equationsCover|(flatten entries gens minors(3,eq33))|(flatten entries gens minors(3,eq34));
-- Flag eqn x1(E22)=0
	M5 := x1*E22;
	equationsCover = equationsCover|(flatten entries M5);
-- NO flag containtment equations
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF CR;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCRCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCRCk
)

subJacCR = () -> (
	M := JacCR();
	use ring M;
	M = substitute(M,{ x_{0,0,0} => 0, x_{0,0,1} => 0, x_{0,0,2} => 1, x_{0,0,3} => 0 });
	M = substitute(M,{ x_{0,1,0} => 0, x_{0,1,1} => 0, x_{0,1,2} => 0, x_{0,1,3} => 1 });
	M = substitute(M,{ x_{1,0,0} => 1, x_{1,0,1} => 0, x_{1,0,2} => 0, x_{1,0,3} => 0 });
	M = substitute(M,{ x_{1,1,0} => 0, x_{1,1,1} => 1, x_{1,1,2} => 0, x_{1,1,3} => 0 });
	M = substitute(M,{ x_{1,2,0} => 0, x_{1,2,1} => 0, x_{1,2,2} => 0, x_{1,2,3} => 0 });
	M = substitute(M,{ x_{1,3,0} => 0, x_{1,3,1} => 0, x_{1,3,2} => 0, x_{1,3,3} => 0 });
	M = substitute(M,{ x_{2,0,0} => 0, x_{2,0,1} => 0, x_{2,0,2} => 1, x_{2,0,3} => 0 });
	M = substitute(M,{ x_{2,1,0} => 0, x_{2,1,1} => 0, x_{2,1,2} => 0, x_{2,1,3} => 1 });
	M = substitute(M,{ x_{2,2,0} => 0, x_{2,2,1} => 0, x_{2,2,2} => 0, x_{2,2,3} => 0 });
	M = substitute(M,{ x_{2,3,0} => 0, x_{2,3,1} => 0, x_{2,3,2} => 0, x_{2,3,3} => 0 });
	M = substitute(M,{ x_{3,0,0} => 1, x_{3,0,1} => 0 });
	M = substitute(M,{ x_{3,1,0} => 0, x_{3,1,1} => 1 });
	M = substitute(M,{ x_{3,2,0} => 0, x_{3,2,1} => 0 });
	M = substitute(M,{ x_{3,3,0} => 0, x_{3,3,1} => 0 });
-- subbing in y0
	M = substitute(M,{ y_{0,0,0} => 1, y_{0,0,1} => 0} );
	M = substitute(M,{ y_{0,1,0} => 0, y_{0,1,1} => 1 });
	M = substitute(M,{ y_{0,2,0} => 0, y_{0,2,1} => 0 });
	M = substitute(M,{ y_{0,3,0} => 0, y_{0,3,1} => 0 });
-- subbing in y1
	M = substitute(M,{ y_{1,0,0} => 0, y_{1,0,1} => 0, y_{1,0,2} => 1, y_{1,0,3} => 0 });
	M = substitute(M,{ y_{1,1,0} => 0, y_{1,1,1} => 0, y_{1,1,2} => 0, y_{1,1,3} => 1 });
	M = substitute(M,{ y_{1,2,0} => 0, y_{1,2,1} => 0, y_{1,2,2} => 1, y_{1,2,3} => 0 });
	M = substitute(M,{ y_{1,3,0} => 0, y_{1,3,1} => 0, y_{1,3,2} => 0, y_{1,3,3} => 1 });
-- subbing in y2
	M = substitute(M,{ y_{2,0,0} => 1, y_{2,0,1} => 0, y_{2,0,3} => 1 });
	M = substitute(M,{ y_{2,1,0} => 0, y_{2,1,1} => 1, y_{2,1,2} => 0});
	M = substitute(M,{ y_{2,2,0} => 0, y_{2,2,1} => 0, y_{2,2,2} => 0, y_{2,2,3} => 0 });
	M = substitute(M,{ y_{2,3,0} => 0, y_{2,3,1} => 0, y_{2,3,2} => 0, y_{2,3,3} => 0 });
-- subbing in y3
	M = substitute(M,{ y_{3,0,0} => 0, y_{3,0,1} => 0, y_{3,0,2} => 1, y_{3,0,3} => 0 });
	M = substitute(M,{ y_{3,1,0} => 0, y_{3,1,1} => 0, y_{3,1,2} => 0, y_{3,1,3} =>1 });
-- subbing in fibre flag variables
	M = substitute(M,{ a1=>0, a2=>0, a3=>0, a4=>0, b1=>0, b2=>0, b3=>0, b4=>0, c1=>0, c2=>0, c3=>0, c4=>0 });
	M
)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- Ev computation for Cr
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

JacCr = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cr = new RankConditions from ({2,4,4,4,2},{{2,2,3,2},{0,2,1},{0,0},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
 	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{1},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,1},{c1,c2,c3}};
	E23 = sub(E23, bigRing);
	E32 := matrix{{1,0},{0,1},{d1,d2},{d3,d4}};
	E32 = sub(E32, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),b2},{M56_(2,0),1},{M56_(3,0),b3}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),b2},{M56_(2,1),1},{M56_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x2 \subset E23
	eq1  := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,1},{x_{2,3,0},c1,c2,c3}};
	eq2  := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,1},{x_{2,3,1},c1,c2,c3}};
	eq3  := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,1},{x_{2,3,2},c1,c2,c3}};
	eq4  := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,1},{x_{2,3,3},c1,c2,c3}};
	equationsCover = append(equationsCover, det eq1);
	equationsCover = append(equationsCover, det eq2);
	equationsCover = append(equationsCover, det eq3);
	equationsCover = append(equationsCover, det eq4);
-- Flag eqn x0(E32)=0
	M3 := x0*E32;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x1(E21)=0
	M5 := x1*E21;
	equationsCover = equationsCover|(flatten entries M5);
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq66 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),d1,d2},{M6_(3,0),d3,d4}};
	eq67 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),d1,d2},{M6_(3,1),d3,d4}};
	eq68 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),d1,d2},{M6_(3,2),d3,d4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq66))|(flatten entries gens minors(3,eq67))|(flatten entries gens minors(3,eq68));
-- Flag eqn x1 \subset E32
	eq76 := matrix{{x_{1,0,0},1,0},{x_{1,1,0},0,1},{x_{1,2,0},d1,d2},{x_{1,3,0},d3,d4}};
	eq77 := matrix{{x_{1,0,1},1,0},{x_{1,1,1},0,1},{x_{1,2,1},d1,d2},{x_{1,3,1},d3,d4}};
	eq78 := matrix{{x_{1,0,2},1,0},{x_{1,1,2},0,1},{x_{1,2,2},d1,d2},{x_{1,3,2},d3,d4}};
	eq79 := matrix{{x_{1,0,3},1,0},{x_{1,1,3},0,1},{x_{1,2,3},d1,d2},{x_{1,3,3},d3,d4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq76))|(flatten entries gens minors(3,eq77));
	equationsCover = equationsCover|(flatten entries gens minors(3,eq78))|(flatten entries gens minors(3,eq79));
-- flag containtment equations
	con := matrix{{b1,1,0,0},{b2,0,1,0},{1,0,0,1},{b3,c1,c2,c3}};
	equationsCover = append(equationsCover, det con);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cr;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCrCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCrCk
)	

subJacCr = () -> (
	M := JacCr();
	use ring M;
	M = substitute(M,{ x_{0,0,0} => 0, x_{0,0,1} => 0, x_{0,0,2} => 1, x_{0,0,3} => 0 });
	M = substitute(M,{ x_{0,1,0} => 0, x_{0,1,1} => 0, x_{0,1,2} => 0, x_{0,1,3} => 1 });
	M = substitute(M,{ x_{1,0,0} => 1, x_{1,0,1} => 0, x_{1,0,2} => 0, x_{1,0,3} => 0 });
	M = substitute(M,{ x_{1,1,0} => 0, x_{1,1,1} => 1, x_{1,1,2} => 0, x_{1,1,3} => 0 });
	M = substitute(M,{ x_{1,2,0} => 0, x_{1,2,1} => 0, x_{1,2,2} => 0, x_{1,2,3} => 0 });
	M = substitute(M,{ x_{1,3,0} => 0, x_{1,3,1} => 0, x_{1,3,2} => 0, x_{1,3,3} => 0 });
	M = substitute(M,{ x_{2,0,0} => 0, x_{2,0,1} => 0, x_{2,0,2} => 1, x_{2,0,3} => 0 });
	M = substitute(M,{ x_{2,1,0} => 0, x_{2,1,1} => 0, x_{2,1,2} => 0, x_{2,1,3} => 1 });
	M = substitute(M,{ x_{2,2,0} => 0, x_{2,2,1} => 0, x_{2,2,2} => 0, x_{2,2,3} => 0 });
	M = substitute(M,{ x_{2,3,0} => 0, x_{2,3,1} => 0, x_{2,3,2} => 0, x_{2,3,3} => 0 });
	M = substitute(M,{ x_{3,0,0} => 1, x_{3,0,1} => 0 });
	M = substitute(M,{ x_{3,1,0} => 0, x_{3,1,1} => 1 });
	M = substitute(M,{ x_{3,2,0} => 0, x_{3,2,1} => 0 });
	M = substitute(M,{ x_{3,3,0} => 0, x_{3,3,1} => 0 });
-- subbing in y0
	M = substitute(M,{ y_{0,0,0} => 1, y_{0,0,1} => 0} );
	M = substitute(M,{ y_{0,1,0} => 0, y_{0,1,1} => 1 });
	M = substitute(M,{ y_{0,2,0} => 0, y_{0,2,1} => 0 });
	M = substitute(M,{ y_{0,3,0} => 0, y_{0,3,1} => 0 });
-- subbing in y1
	M = substitute(M,{ y_{1,0,0} => 0, y_{1,0,1} => 0, y_{1,0,2} => 1, y_{1,0,3} => 0 });
	M = substitute(M,{ y_{1,1,0} => 0, y_{1,1,1} => 0, y_{1,1,2} => 0, y_{1,1,3} => 1 });
	M = substitute(M,{ y_{1,2,0} => 0, y_{1,2,1} => 0, y_{1,2,2} => 1, y_{1,2,3} => 0 });
	M = substitute(M,{ y_{1,3,0} => 0, y_{1,3,1} => 0, y_{1,3,2} => 0, y_{1,3,3} => 1 });
-- subbing in y2
	M = substitute(M,{ y_{2,0,0} => 1, y_{2,0,1} => 0, y_{2,0,3} => 1 });
	M = substitute(M,{ y_{2,1,0} => 0, y_{2,1,1} => 1, y_{2,1,2} => 0});
	M = substitute(M,{ y_{2,2,0} => 0, y_{2,2,1} => 0, y_{2,2,2} => 0, y_{2,2,3} => 0 });
	M = substitute(M,{ y_{2,3,0} => 0, y_{2,3,1} => 0, y_{2,3,2} => 0, y_{2,3,3} => 0 });
-- subbing in y3
	M = substitute(M,{ y_{3,0,0} => 0, y_{3,0,1} => 0, y_{3,0,2} => 1, y_{3,0,3} => 0 });
	M = substitute(M,{ y_{3,1,0} => 0, y_{3,1,1} => 0, y_{3,1,2} => 0, y_{3,1,3} => 1 });
-- subbing in fibre flag variables
	M = substitute(M,{ a1=>0, a2=>0, a3=>0, a4=>0, b1=>0, b2=>0, b3=>c3, c1=>0, c2=>0, d1=>0, d2=>0, d3=>0, d4=>0 });
	M
)

-- for checking smoothness at pt at infinity [a:b] = [0:1]. 
			
InfinityJacCr2 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cr = new RankConditions from ({2,4,4,4,2},{{2,2,3,2},{0,2,1},{0,0},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
 	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{b3},{1}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,c3},{c1,c2,1}};
	E23 = sub(E23, bigRing);
	E32 := matrix{{1,0},{0,1},{d1,d2},{d3,d4}};
	E32 = sub(E32, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),b2},{M56_(2,0),b3},{M56_(3,0),1}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),b2},{M56_(2,1),b3},{M56_(3,1),1}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x2 \subset E23
	eq1  := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,c3},{x_{2,3,0},c1,c2,1}};
	eq2  := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,c3},{x_{2,3,1},c1,c2,1}};
	eq3  := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,c3},{x_{2,3,2},c1,c2,1}};
	eq4  := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,c3},{x_{2,3,3},c1,c2,1}};
	equationsCover = append(equationsCover, det eq1);
	equationsCover = append(equationsCover, det eq2);
	equationsCover = append(equationsCover, det eq3);
	equationsCover = append(equationsCover, det eq4);
-- Flag eqn x0(E32)=0
	M3 := x0*E32;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x1(E21)=0
	M5 := x1*E21;
	equationsCover = equationsCover|(flatten entries M5);
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq66 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),d1,d2},{M6_(3,0),d3,d4}};
	eq67 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),d1,d2},{M6_(3,1),d3,d4}};
	eq68 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),d1,d2},{M6_(3,2),d3,d4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq66))|(flatten entries gens minors(3,eq67))|(flatten entries gens minors(3,eq68));
-- Flag eqn x1 \subset E32
	eq76 := matrix{{x_{1,0,0},1,0},{x_{1,1,0},0,1},{x_{1,2,0},d1,d2},{x_{1,3,0},d3,d4}};
	eq77 := matrix{{x_{1,0,1},1,0},{x_{1,1,1},0,1},{x_{1,2,1},d1,d2},{x_{1,3,1},d3,d4}};
	eq78 := matrix{{x_{1,0,2},1,0},{x_{1,1,2},0,1},{x_{1,2,2},d1,d2},{x_{1,3,2},d3,d4}};
	eq79 := matrix{{x_{1,0,3},1,0},{x_{1,1,3},0,1},{x_{1,2,3},d1,d2},{x_{1,3,3},d3,d4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq76))|(flatten entries gens minors(3,eq77));
	equationsCover = equationsCover|(flatten entries gens minors(3,eq78))|(flatten entries gens minors(3,eq79));
-- flag containtment equations
	con := matrix{{b1,1,0,0},{b2,0,1,0},{b3,0,0,c3},{1,c1,c2,1}};
	equationsCover = append(equationsCover, det con);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cr;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCrCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCrCk
)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- Ev computation for Cm
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

JacCm = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cm = new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{0,0},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4,e1,e2,e3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{1},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,1},{c1,c2,c3}};
	E23 = sub(E23, bigRing);
	E32 := matrix{{1,0},{0,1},{d1,d2},{d3,d4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,1},{e1,e2,e3}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{1},{g}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2(E12) \subset E21
	M1 := x2*E12;
	eq13 := matrix{{M1_(0,0),b1},{M1_(1,0),b2},{M1_(2,0),1},{M1_(3,0),b3}};
	eq14 := matrix{{M1_(0,1),b1},{M1_(1,1),b2},{M1_(2,1),1},{M1_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq13))|(flatten entries gens minors(2,eq14));
-- Flag eqn x2 \subset E23
	eq31 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,1},{x_{2,3,0},c1,c2,c3}};
	eq32 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,1},{x_{2,3,1},c1,c2,c3}};
	eq33 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,1},{x_{2,3,2},c1,c2,c3}};
	eq34 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,1},{x_{2,3,3},c1,c2,c3}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) = 0
	M2 := x1*E21;
	equationsCover = equationsCover|(flatten entries M2);
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq66 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),d1,d2},{M6_(3,0),d3,d4}};
	eq67 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),d1,d2},{M6_(3,1),d3,d4}};
	eq68 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),d1,d2},{M6_(3,2),d3,d4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq66))|(flatten entries gens minors(3,eq67))|(flatten entries gens minors(3,eq68));
-- Flag eqn x1 \subset E33
	eq1  := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,1},{x_{1,3,0},e1,e2,e3}};
	eq2  := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,1},{x_{1,3,1},e1,e2,e3}};
	eq3  := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,1},{x_{1,3,2},e1,e2,e3}};
	eq4  := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,1},{x_{1,3,3},e1,e2,e3}};
	equationsCover = append(equationsCover, det eq1);
	equationsCover = append(equationsCover, det eq2);
	equationsCover = append(equationsCover, det eq3);
	equationsCover = append(equationsCover, det eq4);
-- Flag eqn x0(E32) = 0
	M4 := x0*E32;
	equationsCover = equationsCover|(flatten entries M4);
-- Flag eqn x0(E33) \subset E41
	M5 := x0*E33;
	eq5  := matrix{{M5_(0,0),1},{M5_(1,0),g}};
	eq6  := matrix{{M5_(0,1),1},{M5_(1,1),g}};
	eq7  := matrix{{M5_(0,2),1},{M5_(1,2),g}};
	equationsCover = append(equationsCover, det eq5);
	equationsCover = append(equationsCover, det eq6);
	equationsCover = append(equationsCover, det eq7);
-- Flag containtment equations
	M9 := matrix{{b1,1,0,0},{b2,0,1,0},{1,0,0,1},{b3,c1,c2,c3}};
	equationsCover = append(equationsCover, det M9);
	M7 := matrix{{1,1,0,0},{0,0,1,0},{d1,0,0,1},{d3,e1,e2,e3}};
	M8 := matrix{{0,1,0,0},{1,0,1,0},{d2,0,0,1},{d4,e1,e2,e3}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = append(equationsCover, det M8);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cm;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCmCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCmCk
)

subJacCm = () -> (
	M := JacCm();
	use ring M;
	M = substitute(M,{ x_{0,0,0} => 0, x_{0,0,1} => 0, x_{0,0,2} => 1, x_{0,0,3} => 0 });
	M = substitute(M,{ x_{0,1,0} => 0, x_{0,1,1} => 0, x_{0,1,2} => 0, x_{0,1,3} => 1 });
	M = substitute(M,{ x_{1,0,0} => 1, x_{1,0,1} => 0, x_{1,0,2} => 0, x_{1,0,3} => 0 });
	M = substitute(M,{ x_{1,1,0} => 0, x_{1,1,1} => 1, x_{1,1,2} => 0, x_{1,1,3} => 0 });
	M = substitute(M,{ x_{1,2,0} => 0, x_{1,2,1} => 0, x_{1,2,2} => 0, x_{1,2,3} => 0 });
	M = substitute(M,{ x_{1,3,0} => 0, x_{1,3,1} => 0, x_{1,3,2} => 0, x_{1,3,3} => 0 });
	M = substitute(M,{ x_{2,0,0} => 0, x_{2,0,1} => 0, x_{2,0,2} => 1, x_{2,0,3} => 0 });
	M = substitute(M,{ x_{2,1,0} => 0, x_{2,1,1} => 0, x_{2,1,2} => 0, x_{2,1,3} => 1 });
	M = substitute(M,{ x_{2,2,0} => 0, x_{2,2,1} => 0, x_{2,2,2} => 0, x_{2,2,3} => 0 });
	M = substitute(M,{ x_{2,3,0} => 0, x_{2,3,1} => 0, x_{2,3,2} => 0, x_{2,3,3} => 0 });
	M = substitute(M,{ x_{3,0,0} => 1, x_{3,0,1} => 0 });
	M = substitute(M,{ x_{3,1,0} => 0, x_{3,1,1} => 1 });
	M = substitute(M,{ x_{3,2,0} => 0, x_{3,2,1} => 0 });
	M = substitute(M,{ x_{3,3,0} => 0, x_{3,3,1} => 0 });
-- subbing in y0
	M = substitute(M,{ y_{0,0,0} => 1, y_{0,0,1} => 0} );
	M = substitute(M,{ y_{0,1,0} => 0, y_{0,1,1} => 1 });
	M = substitute(M,{ y_{0,2,0} => 0, y_{0,2,1} => 0 });
	M = substitute(M,{ y_{0,3,0} => 0, y_{0,3,1} => 0 });
-- subbing in y1
	M = substitute(M,{ y_{1,0,0} => 0, y_{1,0,1} => 0, y_{1,0,2} => 1, y_{1,0,3} => 0 });
	M = substitute(M,{ y_{1,1,0} => 0, y_{1,1,1} => 0, y_{1,1,2} => 0, y_{1,1,3} => 1 });
	M = substitute(M,{ y_{1,2,0} => 0, y_{1,2,1} => 0, y_{1,2,2} => 1, y_{1,2,3} => 0 });
	M = substitute(M,{ y_{1,3,0} => 0, y_{1,3,1} => 0, y_{1,3,2} => 0, y_{1,3,3} => 1 });
-- subbing in y2
	M = substitute(M,{ y_{2,0,0} => 1, y_{2,0,1} => 0, y_{2,0,3} => 1 });
	M = substitute(M,{ y_{2,1,0} => 0, y_{2,1,1} => 1, y_{2,1,2} => 0});
	M = substitute(M,{ y_{2,2,0} => 0, y_{2,2,1} => 0, y_{2,2,2} => 0, y_{2,2,3} => 0 });
	M = substitute(M,{ y_{2,3,0} => 0, y_{2,3,1} => 0, y_{2,3,2} => 0, y_{2,3,3} => 0 });
-- subbing in y3
	M = substitute(M,{ y_{3,0,0} => 0, y_{3,0,1} => 0, y_{3,0,2} => 1, y_{3,0,3} => 0 });
	M = substitute(M,{ y_{3,1,0} => 0, y_{3,1,1} => 0, y_{3,1,2} => 0, y_{3,1,3} =>1 });
-- subbing in fibre flag variables
	M = substitute(M,{ a1=>0, a2=>0, a3=>0, a4=>0, b1=>0, b2=>0, b3=>c3, c1=>0, c2=>0, d1=>0, d2=>0, d3=>0, d4=>0, e1=>0, e2=>0, e3=>g });
	M
)

-- (2) [a:b] = [1:0] and [c:d]=[0:1]

 InfinityJacCm2 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cm = new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{0,0},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4,e1,e2,e3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{1},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,1},{c1,c2,c3}};
	E23 = sub(E23, bigRing);
	E32 := matrix{{1,0},{0,1},{d1,d2},{d3,d4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,e3},{e1,e2,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2(E12) \subset E21
	M1 := x2*E12;
	eq13 := matrix{{M1_(0,0),b1},{M1_(1,0),b2},{M1_(2,0),1},{M1_(3,0),b3}};
	eq14 := matrix{{M1_(0,1),b1},{M1_(1,1),b2},{M1_(2,1),1},{M1_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq13))|(flatten entries gens minors(2,eq14));
-- Flag eqn x2 \subset E23
	eq31 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,1},{x_{2,3,0},c1,c2,c3}};
	eq32 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,1},{x_{2,3,1},c1,c2,c3}};
	eq33 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,1},{x_{2,3,2},c1,c2,c3}};
	eq34 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,1},{x_{2,3,3},c1,c2,c3}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) = 0
	M2 := x1*E21;
	equationsCover = equationsCover|(flatten entries M2);
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq66 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),d1,d2},{M6_(3,0),d3,d4}};
	eq67 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),d1,d2},{M6_(3,1),d3,d4}};
	eq68 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),d1,d2},{M6_(3,2),d3,d4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq66))|(flatten entries gens minors(3,eq67))|(flatten entries gens minors(3,eq68));
-- Flag eqn x1 \subset E33
	eq1  := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,e3},{x_{1,3,0},e1,e2,1}};
	eq2  := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,e3},{x_{1,3,1},e1,e2,1}};
	eq3  := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,e3},{x_{1,3,2},e1,e2,1}};
	eq4  := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,e3},{x_{1,3,3},e1,e2,1}};
	equationsCover = append(equationsCover, det eq1);
	equationsCover = append(equationsCover, det eq2);
	equationsCover = append(equationsCover, det eq3);
	equationsCover = append(equationsCover, det eq4);
-- Flag eqn x0(E32) = 0
	M4 := x0*E32;
	equationsCover = equationsCover|(flatten entries M4);
-- Flag eqn x0(E33) \subset E41
	M5 := x0*E33;
	eq5  := matrix{{M5_(0,0),g},{M5_(1,0),1}};
	eq6  := matrix{{M5_(0,1),g},{M5_(1,1),1}};
	eq7  := matrix{{M5_(0,2),g},{M5_(1,2),1}};
	equationsCover = append(equationsCover, det eq5);
	equationsCover = append(equationsCover, det eq6);
	equationsCover = append(equationsCover, det eq7);
-- Flag containtment equations
	M9 := matrix{{b1,1,0,0},{b2,0,1,0},{1,0,0,1},{b3,c1,c2,c3}};
	equationsCover = append(equationsCover, det M9);
	M7 := matrix{{1,1,0,0},{0,0,1,0},{d1,0,0,e3},{d3,e1,e2,1}};
	M8 := matrix{{0,1,0,0},{1,0,1,0},{d2,0,0,e3},{d4,e1,e2,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = append(equationsCover, det M8);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cm;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCmCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCmCk
	
)

 -- (3) for checking smoothness at pt at infinity [a:b]=[0:1] and [c:d] = [1:0]

InfinityJacCm3 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cm = new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{0,0},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4,e1,e2,e3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{b3},{1}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,c3},{c1,c2,1}};
	E23 = sub(E23, bigRing);
	E32 := matrix{{1,0},{0,1},{d1,d2},{d3,d4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,1},{e1,e2,e3}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{1},{g}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2(E12) \subset E21
	M1 := x2*E12;
	eq13 := matrix{{M1_(0,0),b1},{M1_(1,0),b2},{M1_(2,0),b3},{M1_(3,0),1}};
	eq14 := matrix{{M1_(0,1),b1},{M1_(1,1),b2},{M1_(2,1),b3},{M1_(3,1),1}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq13))|(flatten entries gens minors(2,eq14));
-- Flag eqn x2 \subset E23
	eq31 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,c3},{x_{2,3,0},c1,c2,1}};
	eq32 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,c3},{x_{2,3,1},c1,c2,1}};
	eq33 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,c3},{x_{2,3,2},c1,c2,1}};
	eq34 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,c3},{x_{2,3,3},c1,c2,1}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) = 0
	M2 := x1*E21;
	equationsCover = equationsCover|(flatten entries M2);
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq66 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),d1,d2},{M6_(3,0),d3,d4}};
	eq67 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),d1,d2},{M6_(3,1),d3,d4}};
	eq68 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),d1,d2},{M6_(3,2),d3,d4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq66))|(flatten entries gens minors(3,eq67))|(flatten entries gens minors(3,eq68));
-- Flag eqn x1 \subset E33
	eq1  := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,1},{x_{1,3,0},e1,e2,e3}};
	eq2  := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,1},{x_{1,3,1},e1,e2,e3}};
	eq3  := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,1},{x_{1,3,2},e1,e2,e3}};
	eq4  := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,1},{x_{1,3,3},e1,e2,e3}};
	equationsCover = append(equationsCover, det eq1);
	equationsCover = append(equationsCover, det eq2);
	equationsCover = append(equationsCover, det eq3);
	equationsCover = append(equationsCover, det eq4);
-- Flag eqn x0(E32) = 0
	M4 := x0*E32;
	equationsCover = equationsCover|(flatten entries M4);
-- Flag eqn x0(E33) \subset E41
	M5 := x0*E33;
	eq5  := matrix{{M5_(0,0),1},{M5_(1,0),g}};
	eq6  := matrix{{M5_(0,1),1},{M5_(1,1),g}};
	eq7  := matrix{{M5_(0,2),1},{M5_(1,2),g}};
	equationsCover = append(equationsCover, det eq5);
	equationsCover = append(equationsCover, det eq6);
	equationsCover = append(equationsCover, det eq7);
-- Flag containtment equations
	M9 := matrix{{b1,1,0,0},{b2,0,1,0},{b3,0,0,c3},{1,c1,c2,1}};
	equationsCover = append(equationsCover, det M9);
	M7 := matrix{{1,1,0,0},{0,0,1,0},{d1,0,0,1},{d3,e1,e2,e3}};
	M8 := matrix{{0,1,0,0},{1,0,1,0},{d2,0,0,1},{d4,e1,e2,e3}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = append(equationsCover, det M8);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cm;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCmCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCmCk
	
)

-- (4) for checking smoothness at pt at infinity [a:b]=[0:1] and [c:d]=[0:1]

InfinityJacCm4 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cm = new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{0,0},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4,e1,e2,e3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{b3},{1}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,c3},{c1,c2,1}};
	E23 = sub(E23, bigRing);
	E32 := matrix{{1,0},{0,1},{d1,d2},{d3,d4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,e3},{e1,e2,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2(E12) \subset E21
	M1 := x2*E12;
	eq13 := matrix{{M1_(0,0),b1},{M1_(1,0),b2},{M1_(2,0),b3},{M1_(3,0),1}};
	eq14 := matrix{{M1_(0,1),b1},{M1_(1,1),b2},{M1_(2,1),b3},{M1_(3,1),1}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq13))|(flatten entries gens minors(2,eq14));
-- Flag eqn x2 \subset E23
	eq31 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,c3},{x_{2,3,0},c1,c2,1}};
	eq32 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,c3},{x_{2,3,1},c1,c2,1}};
	eq33 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,c3},{x_{2,3,2},c1,c2,1}};
	eq34 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,c3},{x_{2,3,3},c1,c2,1}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) = 0
	M2 := x1*E21;
	equationsCover = equationsCover|(flatten entries M2);
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq66 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),d1,d2},{M6_(3,0),d3,d4}};
	eq67 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),d1,d2},{M6_(3,1),d3,d4}};
	eq68 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),d1,d2},{M6_(3,2),d3,d4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq66))|(flatten entries gens minors(3,eq67))|(flatten entries gens minors(3,eq68));
-- Flag eqn x1 \subset E33
	eq1  := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,e3},{x_{1,3,0},e1,e2,1}};
	eq2  := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,e3},{x_{1,3,1},e1,e2,1}};
	eq3  := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,e3},{x_{1,3,2},e1,e2,1}};
	eq4  := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,e3},{x_{1,3,3},e1,e2,1}};
	equationsCover = append(equationsCover, det eq1);
	equationsCover = append(equationsCover, det eq2);
	equationsCover = append(equationsCover, det eq3);
	equationsCover = append(equationsCover, det eq4);
-- Flag eqn x0(E32) = 0
	M4 := x0*E32;
	equationsCover = equationsCover|(flatten entries M4);
-- Flag eqn x0(E33) \subset E41
	M5 := x0*E33;
	eq5  := matrix{{M5_(0,0),g},{M5_(1,0),1}};
	eq6  := matrix{{M5_(0,1),g},{M5_(1,1),1}};
	eq7  := matrix{{M5_(0,2),g},{M5_(1,2),1}};
	equationsCover = append(equationsCover, det eq5);
	equationsCover = append(equationsCover, det eq6);
	equationsCover = append(equationsCover, det eq7);
-- Flag containtment equations
	M9 := matrix{{b1,1,0,0},{b2,0,1,0},{b3,0,0,c3},{1,c1,c2,1}};
	equationsCover = append(equationsCover, det M9);
	M7 := matrix{{1,1,0,0},{0,0,1,0},{d1,0,0,e3},{d3,e1,e2,1}};
	M8 := matrix{{0,1,0,0},{1,0,1,0},{d2,0,0,e3},{d4,e1,e2,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = append(equationsCover, det M8);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cm;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCmCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCmCk
)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- Ev computation for Cpsi
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

JacCpsi = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{1},{b1},{b2},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,1},{c1,c2,c3}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{1},{d1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,1},{f1,f2,f3}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{1},{g}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,1},{x_{2,3,0},c1,c2,c3}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,1},{x_{2,3,1},c1,c2,c3}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,1},{x_{2,3,2},c1,c2,c3}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,1},{x_{2,3,3},c1,c2,c3}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),1},{M56_(1,0),b1},{M56_(2,0),b2},{M56_(3,0),b3}};
	eq26 := matrix{{M56_(0,1),1},{M56_(1,1),b1},{M56_(2,1),b2},{M56_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),1},{M12_(1,0),g}};
	eq42 := matrix{{M12_(0,1),1},{M12_(1,1),g}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),1},{M456_(1,0),g}};
	eq45 := matrix{{M456_(0,1),1},{M456_(1,1),g}};
	eq46 := matrix{{M456_(0,2),1},{M456_(1,2),g}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,1},{x_{1,3,0},f1,f2,f3}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,1},{x_{1,3,1},f1,f2,f3}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,1},{x_{1,3,2},f1,f2,f3}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,1},{x_{1,3,3},f1,f2,f3}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),1},{M5_(1,0),d1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{1,1,0,0},{b1,0,1,0},{b2,0,0,1},{b3,c1,c2,c3}};
	M8 := matrix{{1,1,0},{d1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,1},{e3,f1,f2,f3}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,1},{e4,f1,f2,f3}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi = () -> (
	M := JacCpsi();
	use ring M;
-- subbing in X_KS
	M = substitute(M,{ x_{0,0,0} => 0, x_{0,0,1} => 0, x_{0,0,2} => 1, x_{0,0,3} => 0 });
	M = substitute(M,{ x_{0,1,0} => 0, x_{0,1,1} => 0, x_{0,1,2} => 0, x_{0,1,3} => 1 });
	M = substitute(M,{ x_{1,0,0} => 1, x_{1,0,1} => 0, x_{1,0,2} => 0, x_{1,0,3} => 0 });
	M = substitute(M,{ x_{1,1,0} => 0, x_{1,1,1} => 1, x_{1,1,2} => 0, x_{1,1,3} => 0 });
	M = substitute(M,{ x_{1,2,0} => 0, x_{1,2,1} => 0, x_{1,2,2} => 0, x_{1,2,3} => 0 });
	M = substitute(M,{ x_{1,3,0} => 0, x_{1,3,1} => 0, x_{1,3,2} => 0, x_{1,3,3} => 0 });
	M = substitute(M,{ x_{2,0,0} => 0, x_{2,0,1} => 0, x_{2,0,2} => 1, x_{2,0,3} => 0 });
	M = substitute(M,{ x_{2,1,0} => 0, x_{2,1,1} => 0, x_{2,1,2} => 0, x_{2,1,3} => 1 });
	M = substitute(M,{ x_{2,2,0} => 0, x_{2,2,1} => 0, x_{2,2,2} => 0, x_{2,2,3} => 0 });
	M = substitute(M,{ x_{2,3,0} => 0, x_{2,3,1} => 0, x_{2,3,2} => 0, x_{2,3,3} => 0 });
	M = substitute(M,{ x_{3,0,0} => 1, x_{3,0,1} => 0 });
	M = substitute(M,{ x_{3,1,0} => 0, x_{3,1,1} => 1 });
	M = substitute(M,{ x_{3,2,0} => 0, x_{3,2,1} => 0 });
	M = substitute(M,{ x_{3,3,0} => 0, x_{3,3,1} => 0 });
-- subbing in y0
	M = substitute(M,{ y_{0,0,0} => 1, y_{0,0,1} => 0} );
	M = substitute(M,{ y_{0,1,0} => 0, y_{0,1,1} => 1 });
	M = substitute(M,{ y_{0,2,0} => 0, y_{0,2,1} => 0 });
	M = substitute(M,{ y_{0,3,0} => 0, y_{0,3,1} => 0 });
-- subbing in y1	
	M = substitute(M,{ y_{1,0,0} => 0, y_{1,0,1} => 0, y_{1,0,2} => 1, y_{1,0,3} => 0 });
	M = substitute(M,{ y_{1,1,0} => 0, y_{1,1,1} => 0, y_{1,1,2} => 0, y_{1,1,3} => 1 });
	M = substitute(M,{ y_{1,2,0} => 0, y_{1,2,1} => 0, y_{1,2,2} => 1, y_{1,2,3} => 0 });
	M = substitute(M,{ y_{1,3,0} => 0, y_{1,3,1} => 0, y_{1,3,2} => 0, y_{1,3,3} => 1 });
-- subbing in y2 (y_{2,0,2} = t_1, y_{2,1,3} = t_2)
        M = substitute(M,{ y_{2,0,0} => 1, y_{2,0,1} => 0, y_{2,0,3} => 1 });
	M = substitute(M,{ y_{2,1,0} => 0, y_{2,1,1} => 1, y_{2,1,2} => 0});
	M = substitute(M,{ y_{2,2,0} => 0, y_{2,2,1} => 0, y_{2,2,2} => 0, y_{2,2,3} => 0 });
	M = substitute(M,{ y_{2,3,0} => 0, y_{2,3,1} => 0, y_{2,3,2} => 0, y_{2,3,3} => 0 });
-- subbing in y3
	M = substitute(M,{ y_{3,0,0} => 0, y_{3,0,1} => 0, y_{3,0,2} => 1, y_{3,0,3} => 0 });
	M = substitute(M,{ y_{3,1,0} => 0, y_{3,1,1} => 0, y_{3,1,2} => 0, y_{3,1,3} =>1 });
-- subbing in fibre flag variables
	M = substitute(M,{ a1=>0, a2=>0, a3=>0, a4=>0});
	M = substitute(M,{ d1=>b1, b3=>c3*b2});
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0, g=>f3});
	M
)


-- (2) for checking smoothness at pt at infinity [a:b]=[1:0], [c:d]=[1:0], [e:h]=[1:0], [beta,alpha]=[0:1]

InfinityJacCpsi2 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{1},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,1},{c1,c2,c3}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{1},{d1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,1},{f1,f2,f3}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{1},{g}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,1},{x_{2,3,0},c1,c2,c3}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,1},{x_{2,3,1},c1,c2,c3}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,1},{x_{2,3,2},c1,c2,c3}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,1},{x_{2,3,3},c1,c2,c3}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),b2},{M56_(2,0),1},{M56_(3,0),b3}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),b2},{M56_(2,1),1},{M56_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),1},{M12_(1,0),g}};
	eq42 := matrix{{M12_(0,1),1},{M12_(1,1),g}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),1},{M456_(1,0),g}};
	eq45 := matrix{{M456_(0,1),1},{M456_(1,1),g}};
	eq46 := matrix{{M456_(0,2),1},{M456_(1,2),g}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,1},{x_{1,3,0},f1,f2,f3}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,1},{x_{1,3,1},f1,f2,f3}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,1},{x_{1,3,2},f1,f2,f3}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,1},{x_{1,3,3},f1,f2,f3}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),1},{M5_(1,0),d1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{b1,1,0,0},{b2,0,1,0},{1,0,0,1},{b3,c1,c2,c3}};
	M8 := matrix{{1,1,0},{d1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,1},{e3,f1,f2,f3}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,1},{e4,f1,f2,f3}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

-- (3) for checking smoothness at pt at infinity [a:b]=[1:0], [c:d]=[1:0], [e:h]=[0:1], [beta,alpha]=[1:0]

InfinityJacCpsi3 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{1},{b1},{b2},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,1},{c1,c2,c3}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{1},{d1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,f3},{f1,f2,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,1},{x_{2,3,0},c1,c2,c3}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,1},{x_{2,3,1},c1,c2,c3}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,1},{x_{2,3,2},c1,c2,c3}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,1},{x_{2,3,3},c1,c2,c3}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),1},{M56_(1,0),b1},{M56_(2,0),b2},{M56_(3,0),b3}};
	eq26 := matrix{{M56_(0,1),1},{M56_(1,1),b1},{M56_(2,1),b2},{M56_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),g},{M12_(1,0),1}};
	eq42 := matrix{{M12_(0,1),g},{M12_(1,1),1}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),g},{M456_(1,0),1}};
	eq45 := matrix{{M456_(0,1),g},{M456_(1,1),1}};
	eq46 := matrix{{M456_(0,2),g},{M456_(1,2),1}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,f3},{x_{1,3,0},f1,f2,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,f3},{x_{1,3,1},f1,f2,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,f3},{x_{1,3,2},f1,f2,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,f3},{x_{1,3,3},f1,f2,1}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),1},{M5_(1,0),d1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{1,1,0,0},{b1,0,1,0},{b2,0,0,1},{b3,c1,c2,c3}};
	M8 := matrix{{1,1,0},{d1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,f3},{e3,f1,f2,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,f3},{e4,f1,f2,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

-- (4) for checking smoothness at pt at infinity [a:b]=[1:0], [c:d]=[1:0], [e:h]=[0:1], [beta,alpha]=[0:1]

InfinityJacCpsi4 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{1},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,1},{c1,c2,c3}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{1},{d1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,f3},{f1,f2,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,1},{x_{2,3,0},c1,c2,c3}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,1},{x_{2,3,1},c1,c2,c3}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,1},{x_{2,3,2},c1,c2,c3}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,1},{x_{2,3,3},c1,c2,c3}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),b2},{M56_(2,0),1},{M56_(3,0),b3}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),b2},{M56_(2,1),1},{M56_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),g},{M12_(1,0),1}};
	eq42 := matrix{{M12_(0,1),g},{M12_(1,1),1}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),g},{M456_(1,0),1}};
	eq45 := matrix{{M456_(0,1),g},{M456_(1,1),1}};
	eq46 := matrix{{M456_(0,2),g},{M456_(1,2),1}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,f3},{x_{1,3,0},f1,f2,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,f3},{x_{1,3,1},f1,f2,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,f3},{x_{1,3,2},f1,f2,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,f3},{x_{1,3,3},f1,f2,1}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),1},{M5_(1,0),d1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{b1,1,0,0},{b1,0,1,0},{1,0,0,1},{b3,c1,c2,c3}};
	M8 := matrix{{1,1,0},{d1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,f3},{e3,f1,f2,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,f3},{e4,f1,f2,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

-- (5) for checking smoothness at pt at infinity [a:b]=[1:0], [c:d]=[0:1], [e:h]=[1:0], [beta,alpha]=[1:0]

InfinityJacCpsi5 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{1},{b1},{b2},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,c3},{c1,c2,1}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{1},{d1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,1},{f1,f2,f3}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{1},{g}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,c3},{x_{2,3,0},c1,c2,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,c3},{x_{2,3,1},c1,c2,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,c3},{x_{2,3,2},c1,c2,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,c3},{x_{2,3,3},c1,c2,1}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),1},{M56_(1,0),b1},{M56_(2,0),b2},{M56_(3,0),b3}};
	eq26 := matrix{{M56_(0,1),1},{M56_(1,1),b1},{M56_(2,1),b2},{M56_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),1},{M12_(1,0),g}};
	eq42 := matrix{{M12_(0,1),1},{M12_(1,1),g}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),1},{M456_(1,0),g}};
	eq45 := matrix{{M456_(0,1),1},{M456_(1,1),g}};
	eq46 := matrix{{M456_(0,2),1},{M456_(1,2),g}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,1},{x_{1,3,0},f1,f2,f3}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,1},{x_{1,3,1},f1,f2,f3}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,1},{x_{1,3,2},f1,f2,f3}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,1},{x_{1,3,3},f1,f2,f3}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),1},{M5_(1,0),d1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{1,1,0,0},{b1,0,1,0},{b2,0,0,c3},{b3,c1,c2,1}};
	M8 := matrix{{1,1,0},{d1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,1},{e3,f1,f2,f3}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,1},{e4,f1,f2,f3}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

-- (6) for checking smoothness at pt at infinity [a:b]=[1:0], [c:d]=[0:1], [e:h]=[1:0], [beta,alpha]=[0:1]

InfinityJacCpsi6 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{b3},{1}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,c3},{c1,c2,1}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{1},{d1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,1},{f1,f2,f3}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{1},{g}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,c3},{x_{2,3,0},c1,c2,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,c3},{x_{2,3,1},c1,c2,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,c3},{x_{2,3,2},c1,c2,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,c3},{x_{2,3,3},c1,c2,1}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),b2},{M56_(2,0),b3},{M56_(3,0),1}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),b2},{M56_(2,1),b3},{M56_(3,1),1}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),1},{M12_(1,0),g}};
	eq42 := matrix{{M12_(0,1),1},{M12_(1,1),g}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),1},{M456_(1,0),g}};
	eq45 := matrix{{M456_(0,1),1},{M456_(1,1),g}};
	eq46 := matrix{{M456_(0,2),1},{M456_(1,2),g}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,1},{x_{1,3,0},f1,f2,f3}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,1},{x_{1,3,1},f1,f2,f3}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,1},{x_{1,3,2},f1,f2,f3}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,1},{x_{1,3,3},f1,f2,f3}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),1},{M5_(1,0),d1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{b1,1,0,0},{b2,0,1,0},{b3,0,0,c3},{1,c1,c2,1}};
	M8 := matrix{{1,1,0},{d1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,1},{e3,f1,f2,f3}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,1},{e4,f1,f2,f3}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

-- (7) for checking smoothness at pt at infinity [a:b]=[1:0], [c:d]=[0:1], [e:h]=[0:1], [beta,alpha]=[1:0]

InfinityJacCpsi7 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{1},{b1},{b2},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,c3},{c1,c2,1}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{1},{d1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,f3},{f1,f2,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,c3},{x_{2,3,0},c1,c2,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,c3},{x_{2,3,1},c1,c2,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,c3},{x_{2,3,2},c1,c2,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,c3},{x_{2,3,3},c1,c2,1}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),1},{M56_(1,0),b1},{M56_(2,0),b2},{M56_(3,0),b3}};
	eq26 := matrix{{M56_(0,1),1},{M56_(1,1),b1},{M56_(2,1),b2},{M56_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),g},{M12_(1,0),1}};
	eq42 := matrix{{M12_(0,1),g},{M12_(1,1),1}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),g},{M456_(1,0),1}};
	eq45 := matrix{{M456_(0,1),g},{M456_(1,1),1}};
	eq46 := matrix{{M456_(0,2),g},{M456_(1,2),1}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,f3},{x_{1,3,0},f1,f2,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,f3},{x_{1,3,1},f1,f2,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,f3},{x_{1,3,2},f1,f2,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,f3},{x_{1,3,3},f1,f2,1}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),1},{M5_(1,0),d1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{1,1,0,0},{b1,0,1,0},{b2,0,0,c3},{b3,c1,c2,1}};
	M8 := matrix{{1,1,0},{d1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,f3},{e3,f1,f2,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,f3},{e4,f1,f2,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

-- (8) for checking smoothness at pt at infinity [a:b]=[1:0], [c:d]=[0:1], [e:h]=[0:1], [beta,alpha]=[0:1]

InfinityJacCpsi8 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{b3},{1}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,c3},{c1,c2,1}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{1},{d1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,f3},{f1,f2,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,c3},{x_{2,3,0},c1,c2,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,c3},{x_{2,3,1},c1,c2,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,c3},{x_{2,3,2},c1,c2,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,c3},{x_{2,3,3},c1,c2,1}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),b2},{M56_(2,0),b3},{M56_(3,0),1}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),b2},{M56_(2,1),b3},{M56_(3,1),1}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),g},{M12_(1,0),1}};
	eq42 := matrix{{M12_(0,1),g},{M12_(1,1),1}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),g},{M456_(1,0),1}};
	eq45 := matrix{{M456_(0,1),g},{M456_(1,1),1}};
	eq46 := matrix{{M456_(0,2),g},{M456_(1,2),1}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,f3},{x_{1,3,0},f1,f2,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,f3},{x_{1,3,1},f1,f2,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,f3},{x_{1,3,2},f1,f2,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,f3},{x_{1,3,3},f1,f2,1}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),1},{M5_(1,0),d1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{b1,1,0,0},{b2,0,1,0},{b3,0,0,c3},{1,c1,c2,1}};
	M8 := matrix{{1,1,0},{d1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,f3},{e3,f1,f2,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,f3},{e4,f1,f2,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

-- (9) for checking smoothness at pt at infinity [a:b]=[0:1], [c:d]=[1:0], [e:h]=[1:0], [beta,alpha]=[1:0]

InfinityJacCpsi9 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{1},{b2},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,1},{c1,c2,c3}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{d1},{1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,1},{f1,f2,f3}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{1},{g}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,1},{x_{2,3,0},c1,c2,c3}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,1},{x_{2,3,1},c1,c2,c3}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,1},{x_{2,3,2},c1,c2,c3}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,1},{x_{2,3,3},c1,c2,c3}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),1},{M56_(2,0),b2},{M56_(3,0),b3}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),1},{M56_(2,1),b2},{M56_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),1},{M12_(1,0),g}};
	eq42 := matrix{{M12_(0,1),1},{M12_(1,1),g}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),1},{M456_(1,0),g}};
	eq45 := matrix{{M456_(0,1),1},{M456_(1,1),g}};
	eq46 := matrix{{M456_(0,2),1},{M456_(1,2),g}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,1},{x_{1,3,0},f1,f2,f3}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,1},{x_{1,3,1},f1,f2,f3}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,1},{x_{1,3,2},f1,f2,f3}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,1},{x_{1,3,3},f1,f2,f3}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),d1},{M5_(1,0),1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{b1,1,0,0},{1,0,1,0},{b2,0,0,1},{b3,c1,c2,c3}};
	M8 := matrix{{d1,1,0},{1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,1},{e3,f1,f2,f3}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,1},{e4,f1,f2,f3}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

-- (10) for checking smoothness at pt at infinity [a:b]=[0:1], [c:d]=[1:0], [e:h]=[1:0], [beta,alpha]=[0:1]

InfinityJacCpsi10 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{1},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,1},{c1,c2,c3}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{d1},{1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,1},{f1,f2,f3}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{1},{g}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,1},{x_{2,3,0},c1,c2,c3}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,1},{x_{2,3,1},c1,c2,c3}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,1},{x_{2,3,2},c1,c2,c3}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,1},{x_{2,3,3},c1,c2,c3}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),b2},{M56_(2,0),1},{M56_(3,0),b3}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),b2},{M56_(2,1),1},{M56_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),1},{M12_(1,0),g}};
	eq42 := matrix{{M12_(0,1),1},{M12_(1,1),g}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),1},{M456_(1,0),g}};
	eq45 := matrix{{M456_(0,1),1},{M456_(1,1),g}};
	eq46 := matrix{{M456_(0,2),1},{M456_(1,2),g}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,1},{x_{1,3,0},f1,f2,f3}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,1},{x_{1,3,1},f1,f2,f3}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,1},{x_{1,3,2},f1,f2,f3}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,1},{x_{1,3,3},f1,f2,f3}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),d1},{M5_(1,0),1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{b1,1,0,0},{b2,0,1,0},{1,0,0,1},{b3,c1,c2,c3}};
	M8 := matrix{{d1,1,0},{1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,1},{e3,f1,f2,f3}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,1},{e4,f1,f2,f3}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

-- (11) for checking smoothness at pt at infinity [a:b]=[0:1], [c:d]=[1:0], [e:h]=[0:1], [beta,alpha]=[1:0]

InfinityJacCpsi11 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{1},{b2},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,1},{c1,c2,c3}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{d1},{1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,f3},{f1,f2,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,1},{x_{2,3,0},c1,c2,c3}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,1},{x_{2,3,1},c1,c2,c3}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,1},{x_{2,3,2},c1,c2,c3}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,1},{x_{2,3,3},c1,c2,c3}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),1},{M56_(2,0),b2},{M56_(3,0),b3}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),1},{M56_(2,1),b2},{M56_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),g},{M12_(1,0),1}};
	eq42 := matrix{{M12_(0,1),g},{M12_(1,1),1}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),g},{M456_(1,0),1}};
	eq45 := matrix{{M456_(0,1),g},{M456_(1,1),1}};
	eq46 := matrix{{M456_(0,2),g},{M456_(1,2),1}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,f3},{x_{1,3,0},f1,f2,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,f3},{x_{1,3,1},f1,f2,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,f3},{x_{1,3,2},f1,f2,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,f3},{x_{1,3,3},f1,f2,1}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),d1},{M5_(1,0),1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{b1,1,0,0},{1,0,1,0},{b2,0,0,1},{b3,c1,c2,c3}};
	M8 := matrix{{d1,1,0},{1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,f3},{e3,f1,f2,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,f3},{e4,f1,f2,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

-- (12) for checking smoothness at pt at infinity [a:b]=[0:1], [c:d]=[1:0], [e:h]=[0:1], [beta,alpha]=[0:1]

InfinityJacCpsi12 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{1},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,1},{c1,c2,c3}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{d1},{1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,f3},{f1,f2,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,1},{x_{2,3,0},c1,c2,c3}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,1},{x_{2,3,1},c1,c2,c3}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,1},{x_{2,3,2},c1,c2,c3}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,1},{x_{2,3,3},c1,c2,c3}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),b2},{M56_(2,0),1},{M56_(3,0),b3}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),b2},{M56_(2,1),1},{M56_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),g},{M12_(1,0),1}};
	eq42 := matrix{{M12_(0,1),g},{M12_(1,1),1}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),g},{M456_(1,0),1}};
	eq45 := matrix{{M456_(0,1),g},{M456_(1,1),1}};
	eq46 := matrix{{M456_(0,2),g},{M456_(1,2),1}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,f3},{x_{1,3,0},f1,f2,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,f3},{x_{1,3,1},f1,f2,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,f3},{x_{1,3,2},f1,f2,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,f3},{x_{1,3,3},f1,f2,1}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),d1},{M5_(1,0),1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{b1,1,0,0},{b2,0,1,0},{1,0,0,1},{b3,c1,c2,c3}};
	M8 := matrix{{d1,1,0},{1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,f3},{e3,f1,f2,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,f3},{e4,f1,f2,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

-- (13) for checking smoothness at pt at infinity [a:b]=[0:1], [c:d]=[0:1], [e:h]=[1:0], [beta,alpha]=[1:0]

InfinityJacCpsi13 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{1},{b2},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,c3},{c1,c2,1}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{d1},{1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,1},{f1,f2,f3}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{1},{g}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,c3},{x_{2,3,0},c1,c2,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,c3},{x_{2,3,1},c1,c2,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,c3},{x_{2,3,2},c1,c2,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,c3},{x_{2,3,3},c1,c2,1}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),1},{M56_(2,0),b2},{M56_(3,0),b3}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),1},{M56_(2,1),b2},{M56_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),1},{M12_(1,0),g}};
	eq42 := matrix{{M12_(0,1),1},{M12_(1,1),g}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),1},{M456_(1,0),g}};
	eq45 := matrix{{M456_(0,1),1},{M456_(1,1),g}};
	eq46 := matrix{{M456_(0,2),1},{M456_(1,2),g}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,1},{x_{1,3,0},f1,f2,f3}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,1},{x_{1,3,1},f1,f2,f3}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,1},{x_{1,3,2},f1,f2,f3}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,1},{x_{1,3,3},f1,f2,f3}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),d1},{M5_(1,0),1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{b1,1,0,0},{1,0,1,0},{b2,0,0,c3},{b3,c1,c2,1}};
	M8 := matrix{{d1,1,0},{1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,1},{e3,f1,f2,f3}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,1},{e4,f1,f2,f3}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

-- (14) for checking smoothness at pt at infinity [a:b]=[0:1], [c:d]=[0:1], [e:h]=[1:0], [beta,alpha]=[1:0]

InfinityJacCpsi14 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{b3},{1}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,c3},{c1,c2,1}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{d1},{1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,1},{f1,f2,f3}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{1},{g}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,c3},{x_{2,3,0},c1,c2,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,c3},{x_{2,3,1},c1,c2,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,c3},{x_{2,3,2},c1,c2,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,c3},{x_{2,3,3},c1,c2,1}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),b2},{M56_(2,0),b3},{M56_(3,0),1}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),b2},{M56_(2,1),b3},{M56_(3,1),1}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),1},{M12_(1,0),g}};
	eq42 := matrix{{M12_(0,1),1},{M12_(1,1),g}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),1},{M456_(1,0),g}};
	eq45 := matrix{{M456_(0,1),1},{M456_(1,1),g}};
	eq46 := matrix{{M456_(0,2),1},{M456_(1,2),g}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,1},{x_{1,3,0},f1,f2,f3}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,1},{x_{1,3,1},f1,f2,f3}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,1},{x_{1,3,2},f1,f2,f3}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,1},{x_{1,3,3},f1,f2,f3}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),d1},{M5_(1,0),1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{b1,1,0,0},{b2,0,1,0},{b3,0,0,c3},{1,c1,c2,1}};
	M8 := matrix{{d1,1,0},{1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,1},{e3,f1,f2,f3}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,1},{e4,f1,f2,f3}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

-- (15) for checking smoothness at pt at infinity [a:b]=[0:1], [c:d]=[0:1], [e:h]=[0:1], [beta,alpha]=[1:0]

InfinityJacCpsi15 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{1},{b2},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,c3},{c1,c2,1}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{d1},{1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,f3},{f1,f2,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,c3},{x_{2,3,0},c1,c2,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,c3},{x_{2,3,1},c1,c2,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,c3},{x_{2,3,2},c1,c2,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,c3},{x_{2,3,3},c1,c2,1}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),1},{M56_(2,0),b2},{M56_(3,0),b3}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),1},{M56_(2,1),b2},{M56_(3,1),b3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),g},{M12_(1,0),1}};
	eq42 := matrix{{M12_(0,1),g},{M12_(1,1),1}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),g},{M456_(1,0),1}};
	eq45 := matrix{{M456_(0,1),g},{M456_(1,1),1}};
	eq46 := matrix{{M456_(0,2),g},{M456_(1,2),1}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,f3},{x_{1,3,0},f1,f2,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,f3},{x_{1,3,1},f1,f2,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,f3},{x_{1,3,2},f1,f2,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,f3},{x_{1,3,3},f1,f2,1}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),d1},{M5_(1,0),1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{b1,1,0,0},{1,0,1,0},{b2,0,0,c3},{b3,c1,c2,1}};
	M8 := matrix{{d1,1,0},{1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,f3},{e3,f1,f2,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,f3},{e4,f1,f2,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

-- (16) for checking smoothness at pt at infinity [a:b]=[0:1], [c:d]=[0:1], [e:h]=[0:1], [beta,alpha]=[0:1]

InfinityJacCpsi16 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{b3},{1}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,c3},{c1,c2,1}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{d1},{1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{0,0,f3},{f1,f2,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
	E41 = sub(E41, bigRing);
	x0 := sub(x0Matrix(0), bigRing);
	x1 := sub(x1Matrix(0), bigRing);
	x2 := sub(x2Matrix(0), bigRing);
	x3 := sub(x3Matrix(0), bigRing);
	equationsCover :={};
	use bigRing;
-- Flag eqn x3 \subset E12
	eq11 := matrix{{x_{3,0,0},1,0},{x_{3,1,0},0,1},{x_{3,2,0},a1,a2},{x_{3,3,0},a3,a4}};
	eq12 := matrix{{x_{3,0,1},1,0},{x_{3,1,1},0,1},{x_{3,2,1},a1,a2},{x_{3,3,1},a3,a4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq11))|(flatten entries gens minors(3,eq12));
-- Flag eqn x2 \subset E23
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},0,0,c3},{x_{2,3,0},c1,c2,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},0,0,c3},{x_{2,3,1},c1,c2,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},0,0,c3},{x_{2,3,2},c1,c2,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},0,0,c3},{x_{2,3,3},c1,c2,1}};
	equationsCover = append(equationsCover, det eq21);
	equationsCover = append(equationsCover, det eq22);
	equationsCover = append(equationsCover, det eq23);
	equationsCover = append(equationsCover, det eq24);
-- Flag eqn x2(E12) \subset E21
	M56 := x2*E12;
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),b2},{M56_(2,0),b3},{M56_(3,0),1}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),b2},{M56_(2,1),b3},{M56_(3,1),1}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x0(E32) \subset E41
	M12 := x0*E32;
	eq41 := matrix{{M12_(0,0),g},{M12_(1,0),1}};
	eq42 := matrix{{M12_(0,1),g},{M12_(1,1),1}};
	equationsCover = append(equationsCover, det eq41);
	equationsCover = append(equationsCover, det eq42);
-- Flag eqn x0(E31)=0
	M3 := x0*E31;
	equationsCover = equationsCover|(flatten entries M3);
-- Flag eqn x0(E33) \subset E41
	M456 := x0*E33;
	eq44 := matrix{{M456_(0,0),g},{M456_(1,0),1}};
	eq45 := matrix{{M456_(0,1),g},{M456_(1,1),1}};
	eq46 := matrix{{M456_(0,2),g},{M456_(1,2),1}};
	equationsCover = append(equationsCover, det eq44);
	equationsCover = append(equationsCover, det eq45);
	equationsCover = append(equationsCover, det eq46);
-- Flag eqn x1 \subset E33
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},0,0,f3},{x_{1,3,0},f1,f2,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},0,0,f3},{x_{1,3,1},f1,f2,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},0,0,f3},{x_{1,3,2},f1,f2,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},0,0,f3},{x_{1,3,3},f1,f2,1}};
	equationsCover = append(equationsCover, det eq31);
	equationsCover = append(equationsCover, det eq32);
	equationsCover = append(equationsCover, det eq33);
	equationsCover = append(equationsCover, det eq34);
-- Flag eqn x1(E21) \subset E31
	M5 := x1*E21;
	eq35 := matrix{{M5_(0,0),d1},{M5_(1,0),1},{M5_(2,0),d2},{M5_(3,0),d3}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq35));
-- Flag eqn x1(E23) \subset E32
	M6 := x1*E23;
	eq36 := matrix{{M6_(0,0),1,0},{M6_(1,0),0,1},{M6_(2,0),e1,e2},{M6_(3,0),e3,e4}};
	eq37 := matrix{{M6_(0,1),1,0},{M6_(1,1),0,1},{M6_(2,1),e1,e2},{M6_(3,1),e3,e4}};
	eq38 := matrix{{M6_(0,2),1,0},{M6_(1,2),0,1},{M6_(2,2),e1,e2},{M6_(3,2),e3,e4}};
	equationsCover = equationsCover|(flatten entries gens minors(3,eq36))|(flatten entries gens minors(3,eq37))|(flatten entries gens minors(3,eq38));
-- flag containtment equations
	M7 := matrix{{b1,1,0,0},{b2,0,1,0},{b3,0,0,c3},{1,c1,c2,1}};
	M8 := matrix{{d1,1,0},{1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,0,0,f3},{e3,f1,f2,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,0,0,f3},{e4,f1,f2,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together!
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

subInfinity = method()
subInfinity(Matrix) := (M) -> (
	use ring M;
-- subbing in x_{KS}
	M = substitute(M,{ x_{0,0,0} => 0, x_{0,0,1} => 0, x_{0,0,2} => 1, x_{0,0,3} => 0 });
	M = substitute(M,{ x_{0,1,0} => 0, x_{0,1,1} => 0, x_{0,1,2} => 0, x_{0,1,3} => 1 });
	M = substitute(M,{ x_{1,0,0} => 1, x_{1,0,1} => 0, x_{1,0,2} => 0, x_{1,0,3} => 0 });
	M = substitute(M,{ x_{1,1,0} => 0, x_{1,1,1} => 1, x_{1,1,2} => 0, x_{1,1,3} => 0 });
	M = substitute(M,{ x_{1,2,0} => 0, x_{1,2,1} => 0, x_{1,2,2} => 0, x_{1,2,3} => 0 });
	M = substitute(M,{ x_{1,3,0} => 0, x_{1,3,1} => 0, x_{1,3,2} => 0, x_{1,3,3} => 0 });
	M = substitute(M,{ x_{2,0,0} => 0, x_{2,0,1} => 0, x_{2,0,2} => 1, x_{2,0,3} => 0 });
	M = substitute(M,{ x_{2,1,0} => 0, x_{2,1,1} => 0, x_{2,1,2} => 0, x_{2,1,3} => 1 });
	M = substitute(M,{ x_{2,2,0} => 0, x_{2,2,1} => 0, x_{2,2,2} => 0, x_{2,2,3} => 0 });
	M = substitute(M,{ x_{2,3,0} => 0, x_{2,3,1} => 0, x_{2,3,2} => 0, x_{2,3,3} => 0 });
	M = substitute(M,{ x_{3,0,0} => 1, x_{3,0,1} => 0 });
	M = substitute(M,{ x_{3,1,0} => 0, x_{3,1,1} => 1 });
	M = substitute(M,{ x_{3,2,0} => 0, x_{3,2,1} => 0 });
	M = substitute(M,{ x_{3,3,0} => 0, x_{3,3,1} => 0 });
-- subbing in y0
	M = substitute(M,{ y_{0,0,0} => 1, y_{0,0,1} => 0} );
	M = substitute(M,{ y_{0,1,0} => 0, y_{0,1,1} => 1 });
	M = substitute(M,{ y_{0,2,0} => 0, y_{0,2,1} => 0 });
	M = substitute(M,{ y_{0,3,0} => 0, y_{0,3,1} => 0 });
-- subbing in y1
	M = substitute(M,{ y_{1,0,0} => 0, y_{1,0,1} => 0, y_{1,0,2} => 1, y_{1,0,3} => 0 });
	M = substitute(M,{ y_{1,1,0} => 0, y_{1,1,1} => 0, y_{1,1,2} => 0, y_{1,1,3} => 1 });
	M = substitute(M,{ y_{1,2,0} => 0, y_{1,2,1} => 0, y_{1,2,2} => 1, y_{1,2,3} => 0 });
	M = substitute(M,{ y_{1,3,0} => 0, y_{1,3,1} => 0, y_{1,3,2} => 0, y_{1,3,3} => 1 });
-- subbing in y2
	M = substitute(M,{ y_{2,0,0} => 1, y_{2,0,1} => 0, y_{2,0,3} => 1 });
	M = substitute(M,{ y_{2,1,0} => 0, y_{2,1,1} => 1, y_{2,1,2} => 0});
	M = substitute(M,{ y_{2,2,0} => 0, y_{2,2,1} => 0, y_{2,2,2} => 0, y_{2,2,3} => 0 });
	M = substitute(M,{ y_{2,3,0} => 0, y_{2,3,1} => 0, y_{2,3,2} => 0, y_{2,3,3} => 0 });
-- subbing in y3
	M = substitute(M,{ y_{3,0,0} => 0, y_{3,0,1} => 0, y_{3,0,2} => 1, y_{3,0,3} => 0 });
	M = substitute(M,{ y_{3,1,0} => 0, y_{3,1,1} => 0, y_{3,1,2} => 0, y_{3,1,3} =>1 });
-- list of variables in M
	V = support M;
-- deleting varaibles y_{2,0,2},y_{2,1,3} from V. The variabes in V are still in "ring M"
	V = delete(y_{2,0,2},V);
	V = delete(y_{2,1,3},V);
-- setting every variable in V to 0
 	for i in 0..(#V - 1) do (
     		 M = substitute(M,{ V_i => 0});
      	);
	M
)





