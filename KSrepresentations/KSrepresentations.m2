------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- File: KSrepresentations.m2
-- Author: Nicole Kitt
--
-- The functions found in this file are used in the computations of Ev in the KS project.
-- The functions are specific to the Vogan variety V in the KS project.
--
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
-- JacCr2()
-- JacCm2()
-- JacCm3()
-- JacCm4()
-- JacCpsi2()
-- JacCpsi3()
-- JacCpsi4()
-- JacCpsi5()
-- JacCpsi6()
-- JacCpsi7()
-- JacCpsi8()
-- JacCpsi9()
-- JacCpsi10()
-- JacCpsi11()
-- JacCpsi12()
-- JacCpsi13()
-- JacCpsi14()
-- JacCpsi15()
-- JacCpsi16()
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The functions x0Matrix(), x1Matrix(), x2Matrix(), x3Matrix() output a fixed matrix. 
-- The tuple (x0Matrix(), x1Matrix(), x2Matrix(), x3Matrix()) corresponds to an arbitrary element (x_4,x_3,x_2,x_1) of V.
-- These functions are used in the following: coverCr(), coverCR(), coverCm(), coverCpsi().
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The functions coverCR(), coverCr(), coverCm(), coverCpsi() output equations that describe an affine chart of the cover space \widetilde{C} for the closures of the
-- orbits C_R,C_r,C_m, and C_\psi, respectively. 
-- These are the main affine charts used in the KS project. The subspace Eij in the code is equal to the subspace E_{\lambda_i}^j in the paper. 
-- Note that (x0,x1,x2,x3) in the code corresponds to (x_4,x_3,x_2,x_1) in the paper.
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The functions smallCR(), smallCr(), smallCm(), smallCpsi() output information that tells us if the cover for the closures of CR, Cr, Cm, and Cpsi,  respectively, 
-- is small or semi-small or neither. 
-- The output consists of three lists. Note that the first list will always contain the orbit itself. 
--  NOTE: The third lists contains orbits that should be manually checked to verify if they should be in list 1 or 2.
--        In practice it is empty.
-- We describe what the first two lists contain below: Let C be one of CR,Cr,Cm Cpsi. And \rho:\widetilde{C}\rightarrow\overline{C} be a cover of \overline{C}.
-- The first list of smallC() consists of all suborbits C' \subseteq \overline{C} such that 2\dim(\rho^{-1}(x)) + dim C' = dim\widetilde{C}, where x\in C' is arbitrary.
-- The second list of smallC() consists of all suborbirs C'\subseteq \overline{C} such that 2\dim(\rho^{-1}(x)) + \dim C' > \dim\widetilde{C}, where x\in C' is arbirary.
-- If the first list only contains the orbit itself and the second list is empty, then the cover is small.
-- If the first list contains orbits other than itself and the second list is empty, then the cover is semi-small.
-- If the second list is non-empty, then the cover is neither small nor semi-small.
-- Note that \dim\widetilde{C} = \dim\overline{C} = \dim C for the H-orbits C = CR,Cr,Cm,Cpsi.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The function JacCR() outputs the Jacobian for the equations that describe the main affine choice of the cover for CR ( i.e., coverCR() ), C_{KS}*, and \tilde{f}_R^{-1}(0). 
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The function JacCr() outputs the Jacobian for the equations that describe the main affine choice of the cover for Cr ( i.e., coverCr() ), C_{KS}*, and \tilde{f}_r^{-1}(0). 
-- The function JacCr2() does the same, but with a different affine choice for the cover.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The function JacCm() outputs the Jacobian for the equations that describe the main affine choice of the cover for Cm ( i.e., coverCm() ), C_{KS}*, and \tilde{f}_m^{-1}(0). 
-- The functions JacCm2(), JacCm3(), JacCm4() do the same, but with a different affine choice for the cover.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The function JacCpsi() outputs the Jacobian for the equations that describe the main cover for Cpsi ( i.e., coverCpsi() ), C_{KS}*, and \tilde{f}_{\psi}^{-1}(0). 
-- The functions JacCpsi2(), JacCpsi3(), JacCpsi4(), JacCpsi5(), JacCpsi6(), JacCpsi7(),JacCpsi8(), JacCpsi9(),
-- JacCpsi10(), JacCpsi11(), JacCpsi12(), JacCpsi13(), JacCpsi14(), JacCpsi15(), JacCpsi16() do the same, but with a different affine choice for the cover.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The function subJacCR() evaluates the Jacobian JacCR() at the generic point (x_{\KS},y_{\KS}(t_1,t_2), with diag(t_1,t_2)\in T'_{reg}, and adds the conditions that
-- describe U\cap \rho_R^{-1}(x_{\KS}. Here U=coverCR(). The only input subJacCR() takes is JacCR().
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The function subJacCr() evaluates the Jacobian JacCr() at the generic point (x_{\KS},y_{\KS}(t_1,t_2), with diag(t_1,t_2)\in T'_{reg}, and adds the conditions that
-- describe U\cap \rho_r^{-1}(x_{\KS}. Here U=coverCr(). 
-- The function subJacCr2() does the same, but describes U\cap\rho_r^{-1}(x_{\KS}) where U is the affine choice of the cover found in JacCr2(). 
-- The only input subJacCr() and subJacCr2() take is JacCr().
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- The function subJacCm() evaluates the Jacobian JacCm() at the generic point (x_{\KS},y_{\KS}(t_1,t_2), with diag(t_1,t_2)\in T'_{reg}, and adds the conditions that
-- describe U\cap \rho_m^{-1}(x_{\KS}. Here U=coverCm(). The only input subJacCm() takes is JacCm().
-- The function subJacCm#() does the same, but describes U\cap\rho_m^{-1}(x_{\KS}) where U is the affine choice of the cover found in JacCm#(). Here #=2-4.
-- The only input subJacCm#() takes is JacCm#().
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
--The function subJacCpsi() evaluates the Jacobian JacCpsi() at the generic point (x_{\KS},y_{\KS}(t_1,t_2), with diag(t_1,t_2)\in T'_{reg}, and adds the conditions that
-- describe U\cap \rho_{psi}^{-1}(x_{\KS}. Here U=coverCpsi(). The only input subJacCpsi() takes is JacCpsi().
-- The function subJacCpsi#() does the same, but describes U\cap\rho_{\psi}^{-1}(x_{\KS}) where U is the affine choice of the cover found in JacCpsi#(). Here #=2-16.
-- The only input subJacCpsi#() takes is JacCpsi#().
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- Note: That (y3,y2,y1,y0) in the code corresponds to (y_1,y_2,y_3,y_4) in the paper.
-- Note: That (x0,x1,x2,x3) in the code corresponds to (x_4,x_3,x_2,z_1) in the paper.
-- Note: That y_{2,0,2} corresponds to t_1 in the paer, and y_{2,1,3} corresponds to t_2 in the paper.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

needs "ComputeDuals.m2" 
needs "NetworkFlow.m2"
needs "ComputeCover.m2" 
needs "VanishingCycles.m2"
needs "VoganV.m2"
needs "PSNF.m2"
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
-- making affine choice for flag subspaces
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E22 := matrix{{b1,b2},{b3,b4},{1,0},{0,1}};
	E22 = sub(E22, bigRing);
	E32 := matrix{{1,0},{0,1},{c1,c2},{c3,c4}};
	E32 = sub(E32, bigRing);
-- pushing matrices into bigRing
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
-- making affine choice for flag subspaces
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{1},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{0,0,1},{c1,c2,c3}};
	E23 = sub(E23, bigRing);
	E32 := matrix{{1,0},{0,1},{d1,d2},{d3,d4}};
	E32 = sub(E32, bigRing);
-- pushing matrices into bigRing
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
-- making affine choice for flag subspaces
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
-- pushing matrices into bigRing
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
-- making affine choice for flag subspaces
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
-- pushing matrices into bigRing
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

-- Idea is that:
--   fibreDim  is a lower bound on the dimension of fiber (because we only look at one affine chart)
--   (pullbackDim - subOrbitDim) is an upper bound on the dimension of the fibre.
--  We compute pullbackDim by treating the pullback as a vector bundle over the flag variety and 
--    finding the rank at a chosen point (this is very fast)
--  We use both an upper and lower bound precisely because pullbackDim gives us
--    dim( \rho^{-1}(\overline{S}) ) when we want dim( \rho_r^{-1}(S) ) and there are
--    orbits where a boundary strata is higher dimensional.
--  NOTE: If the map is semi-small then this higher dimensional sub-strata cant actually trigger
--    us to think a strata is relevant (because the innequality would need to be violated at the substrata).

smallCR = () -> (
CR := new RankConditions from ({2,4,4,4,2},{{2,2,4,2},{0,2,2},{0,0},{0}});
semiList := {};
notList := {};
fList := {};
orbitDim := dim(getEquations CR);
CoverCR := coverCR();
for S in getAllSubstrata(CR) do (
	rep := getOtherMatrixRep S;
	A := QQ[a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4];
        eqS := getEquations S;
	subOrbitDim := dim(eqS);
        bigRing := A ** (ring eqS);

        -- We compute the dimension of the fiber (in our affine chart) over a "random" representative 
	re := gens A;
	for M in rep do re = join(re,flatten(entries(flatten(transpose M))));
	fibre := sub(gens CoverCR,matrix{re});
	fibre = sub(fibre,A);
	fibreDim := dim(ideal fibre);

        -- This computes \rho^{-1}(\overline{S})
        pullback := sub( CoverCR, bigRing) + sub(eqS, bigRing);

        -- This looks at a fibre of the pullback under projection to the partial flag variety
        use ring pullback;
        pullbackfibre := sub( pullback,  {a1=>0,a2=>0,a3=>0,a4=>0,b1=>0,b2=>0,b3=>0,b4=>0,c1=>0,c2=>0,c3=>0,c4=>0} );
        pullbackfibre = sub( pullbackfibre, ring eqS );

        -- The dimension of \rho_r^{-1}(S) is the dimension of the flag variety (12) plus dimension of 
        --     the vector bundle over it.
        pullbackDim := 12 + dim( pullbackfibre );
        --   You can use the following to verify, but note... it is slow.
        -- pullbackDim := dim( pullback );

        -- Sanity check
        if pullbackDim - subOrbitDim < fibreDim then (
           print "Failure:";
           print S;
        );

        if pullbackDim != fibreDim + subOrbitDim then (
        --       print "upper != lower";
               if 2* pullbackDim - subOrbitDim >= orbitDim then (
                    fList = append(fList, S);
               )
        );

	if 2*fibreDim + subOrbitDim == orbitDim then semiList = append(semiList,S);
	if 2*fibreDim + subOrbitDim > orbitDim then notList = append(notList,S);
	);
	(semiList,notList,fList)
)

-- Idea is that:
--   fibreDim  is a lower bound on the dimension of fiber (because we only look at one affine chart)
--   (pullbackDim - subOrbitDim) is an upper bound on the dimension of the fibre.
--  We compute pullbackDim by treating the pullback as a vector bundle over the flag variety and 
--    finding the rank at a chosen point (this is very fast)
--  We use both an upper and lower bound precisely because pullbackDim gives us
--    dim( \rho^{-1}(\overline{S}) ) when we want dim( \rho_r^{-1}(S) ) and there are
--    orbits where a boundary strata is higher dimensional.
--  NOTE: If the map is semi-small then this higher dimensional sub-strata cant actually trigger
--    us to think a strata is relevant (because the innequality would need to be violated at the substrata).

smallCr = () -> (
Cr := new RankConditions from ({2,4,4,4,2},{{2,2,3,2},{0,2,1},{0,0},{0}});
semiList := {};
notList := {};
fList := {};
orbitDim := dim(getEquations Cr);
CoverCr := coverCr();
for S in getAllSubstrata(Cr) do (
	rep := getOtherMatrixRep S;
	A := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4];
        eqS := getEquations S;
	subOrbitDim := dim(eqS);
        bigRing := A ** (ring eqS);

        -- We compute the dimension of the fiber (in our affine chart) over a "random" representative 
	re := gens A;
	for M in rep do re = join(re,flatten(entries(flatten(transpose M))));
	fibre := sub(gens CoverCr,matrix{re});
	fibre = sub(fibre,A);
	fibreDim := dim(ideal fibre);

        -- This computes \rho^{-1}(\overline{S})
        pullback := sub( CoverCr, bigRing) + sub(eqS, bigRing);

        -- This looks at a fibre of the pullback under projection to the partial flag variety
        use ring pullback;
        pullbackfibre := sub( pullback,  {a1=>0,a2=>0,a3=>0,a4=>0,b1=>0,b2=>0,b3=>0,c1=>0,c2=>0,c3=>0,d1=>0,d2=>0,d3=>0,d4=>0} );
        pullbackfibre = sub( pullbackfibre, ring eqS );

        -- The dimension of \rho_r^{-1}(S) is the dimension of the flag variety (13) plus dimension of 
        --     the vector bundle over it.
        pullbackDim := 13 + dim(  pullbackfibre );
        --   You can use the following to verify, but note... it is slow.
        -- pullbackDim := dim( pullback );


        -- Sanity check
        if pullbackDim - subOrbitDim < fibreDim then (
           print "Failure:";
           print S;
        );

        if pullbackDim != fibreDim + subOrbitDim then (
        --       print "upper != lower";
               if 2* pullbackDim - subOrbitDim >= orbitDim then (
                    fList = append(fList, S);
               )
        );

	if 2*fibreDim + subOrbitDim == orbitDim then semiList = append(semiList,S);
	if 2*fibreDim + subOrbitDim > orbitDim then notList = append(notList,S);
	);
	(semiList,notList,fList)
)

-- Idea is that:
--   fibreDim  is a lower bound on the dimension of fiber (because we only look at one affine chart)
--   (pullbackDim - subOrbitDim) is an upper bound on the dimension of the fibre.
--  We compute pullbackDim by treating the pullback as a vector bundle over the flag variety and 
--    finding the rank at a chosen point (this is very fast)
--  We use both an upper and lower bound precisely because pullbackDim gives us
--    dim( \rho^{-1}(\overline{S}) ) when we want dim( \rho_r^{-1}(S) ) and there are
--    orbits where a boundary strata is higher dimensional.
--  NOTE: If the map is semi-small then this higher dimensional sub-strata cant actually trigger
--    us to think a strata is relevant (because the innequality would need to be violated at the substrata).

smallCm = () -> (
Cm := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{0,0},{0}});
semiList := {};
notList := {};
fList := {};
orbitDim := dim(getEquations Cm);
CoverCm := coverCm();
for S in getAllSubstrata(Cm) do (
	rep := getOtherMatrixRep S;
	A := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4,e1,e2,e3,g];
        eqS := getEquations S;
	subOrbitDim := dim(eqS);

        -- We compute the dimension of the fiber (in our affine chart) over a "random" representative 
        bigRing := A ** (ring eqS);
	re := gens A;
	for M in rep do re = join(re,flatten(entries(flatten(transpose M))));
	fibre := sub(gens CoverCm,matrix{re});
	fibre = sub(fibre,A);
	fibreDim := dim(ideal fibre);

        -- This computes \rho^{-1}(\overline{S})
        pullback := sub( CoverCm, bigRing) + sub(eqS, bigRing);

        -- This looks at a fibre of the pullback under projection to the partial flag variety
        use ring pullback;
        pullbackfibre := sub( pullback,  {a1=>0,a2=>0,a3=>0,a4=>0,b1=>0,b2=>0,b3=>0,c1=>0,c2=>0,c3=>0,d1=>0,d2=>0,d3=>0,d4=>0,e1=>0,e2=>0,e3=>0,g=>0} );
        pullbackfibre = sub( pullbackfibre, ring eqS );

        -- The dimension of \rho_r^{-1}(S) is the dimension of the flag variety (15) plus dimension of 
        --     the vector bundle over it.
        pullbackDim := 15 + dim( pullbackfibre );
        --   You can use the following to verify, but note... it is slow.
        -- pullbackDim := dim( pullback );

        -- Sanity check
        if pullbackDim - subOrbitDim < fibreDim then (
           print "Failure:";
           print S;
        );

        if pullbackDim != fibreDim + subOrbitDim then (
        --       print "upper != lower";
               if 2* pullbackDim - subOrbitDim >= orbitDim then (
                    fList = append(fList, S);
               )
        );


	if 2*fibreDim + subOrbitDim == orbitDim then semiList = append(semiList,S);
	if 2*fibreDim + subOrbitDim > orbitDim then notList = append(notList,S);
	);
	(semiList,notList,fList)
)
 
-- Idea is that:
--   fibreDim  is a lower bound on the dimension of fiber (because we only look at one affine chart)
--   (pullbackDim - subOrbitDim) is an upper bound on the dimension of the fibre.
--  We compute pullbackDim by treating the pullback as a vector bundle over the flag variety and 
--    finding the rank at a chosen point (this is very fast)
--  We use both an upper and lower bound precisely because pullbackDim gives us
--    dim( \rho^{-1}(\overline{S}) ) when we want dim( \rho_r^{-1}(S) ) and there are
--    orbits where a boundary strata is higher dimensional.
--  NOTE: If the map is semi-small then this higher dimensional sub-strata cant actually trigger
--    us to think a strata is relevant (because the innequality would need to be violated at the substrata).
smallCpsi = () -> (
Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
semiList := {};
notList := {};
fList := {};
orbitDim := dim(getEquations Cpsi);
CoverCpsi := coverCpsi();
for S in getAllSubstrata(Cpsi) do (
	rep := getOtherMatrixRep S;
	A := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
        eqS := getEquations S;
	subOrbitDim := dim(eqS);

        -- We compute the dimension of the fiber (in our affine chart) over a "random" representative 
        bigRing := A ** (ring eqS);
	re := gens A;
	for M in rep do re = join(re,flatten(entries(flatten(transpose M))));
	fibre := sub(gens CoverCpsi,matrix{re});
	fibre = sub(fibre,A);
	fibreDim := dim(ideal fibre);

        -- This computes \rho^{-1}(\overline{S})
        pullback := sub( CoverCpsi, bigRing) + sub(eqS, bigRing);

        -- This looks at a fibre of the pullback under projection to the partial flag variety
        use ring pullback;
        pullbackfibre := sub( pullback, {a1=>0,a2=>0,a3=>0,a4=>0,b1=>0,b2=>0,b3=>0,c1=>0,c2=>0,c3=>0,d1=>0,d2=>0,d3=>0,e1=>0,e2=>0,e3=>0,e4=>0,f1=>0,f2=>0,f3=>0,g=>0} );
        pullbackfibre = sub( pullbackfibre, ring eqS );

        -- The dimension of \rho_r^{-1}(S) is the dimension of the flag variety (16) plus dimension of 
        --     the vector bundle over it.
        pullbackDim := 16 + dim( pullbackfibre );
        --   You can use the following to verify, but note... it is slow.
        -- pullbackDim := dim( pullback );

        -- Sanity check
        if pullbackDim - subOrbitDim < fibreDim then (
           print "Failure:";
           print S;
        );

        if pullbackDim != fibreDim + subOrbitDim then (
        --       print "upper != lower";
               if 2* pullbackDim - subOrbitDim >= orbitDim then (
                    fList = append(fList, S);
               )
        );

	if 2*fibreDim + subOrbitDim == orbitDim then semiList = append(semiList,S);
	if 2*fibreDim + subOrbitDim > orbitDim then notList = append(notList,S);
	);
	(semiList,notList,fList)
)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- Ev computation for CR
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

-- JacCR(): uses affine chart containing fibre \rho_R^{-1}(x_{\KS}).
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
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF CR;
	F = sub(F, bigRing);
-- Putting all of the equations together
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

-- (1) JacCr(): uses affine chart containing {[1:b]}.

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
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cr;
	F = sub(F, bigRing);
-- Putting all of the equations together
	idealCrCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCrCk
)	

-- (2) JacCr2(): uses affine chart containing {[a:1]}.
			
JacCr2 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cr = new RankConditions from ({2,4,4,4,2},{{2,2,3,2},{0,2,1},{0,0},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
 	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{b3},{1}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{c1,c2,c3},{0,0,1}};
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
	eq25 := matrix{{M56_(0,0),b1},{M56_(1,0),b2},{M56_(2,0),b3},{M56_(3,0),1}};
	eq26 := matrix{{M56_(0,1),b1},{M56_(1,1),b2},{M56_(2,1),b3},{M56_(3,1),1}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq25))|(flatten entries gens minors(2,eq26));
-- Flag eqn x2 \subset E23
	eq1  := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},c1,c2,c3},{x_{2,3,0},0,0,1}};
	eq2  := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},c1,c2,c3},{x_{2,3,1},0,0,1}};
	eq3  := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},c1,c2,c3},{x_{2,3,2},0,0,1}};
	eq4  := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},c1,c2,c3},{x_{2,3,3},0,0,1}};
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
	con := matrix{{b1,1,0,0},{b2,0,1,0},{b3,c1,c2,c3},{1,0,0,1}};
	equationsCover = append(equationsCover, det con);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cr;
	F = sub(F, bigRing);
-- Putting all of the equations together
	idealCrCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCrCk
)

subJacCr = method()
subJacCr(Matrix) := (M) -> (
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

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- Ev computation for Cm
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

-- (1) JacCm():  uses affine chart containing ({[1:b],[1:d]}).

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
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cm;
	F = sub(F, bigRing);
-- Putting all of the equations together
	idealCmCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCmCk
)

-- (2) JacCm2():  uses affine chart containing ({[1:b]},{[c:1]}).

 JacCm2 = () -> (
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
	E33 := matrix{{1,0,0},{0,1,0},{e1,e2,e3},{0,0,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
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
	eq1  := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},e1,e2,e3},{x_{1,3,0},0,0,1}};
	eq2  := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},e1,e2,e3},{x_{1,3,1},0,0,1}};
	eq3  := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},e1,e2,e3},{x_{1,3,2},0,0,1}};
	eq4  := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},e1,e2,e3},{x_{1,3,3},0,0,1}};
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
	M7 := matrix{{1,1,0,0},{0,0,1,0},{d1,e1,e2,e3},{d3,0,0,1}};
	M8 := matrix{{0,1,0,0},{1,0,1,0},{d2,e1,e2,e3},{d4,0,0,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = append(equationsCover, det M8);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cm;
	F = sub(F, bigRing);
-- Putting all of the equations together
	idealCmCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCmCk
	
)

 -- (3) JacCms(): uses affine chart containing ({[a:1]},{[1:d]}).

JacCm3 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cm = new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{0,0},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4,e1,e2,e3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{b3},{1}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{c1,c2,c3},{0,0,1}};
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
	eq13 := matrix{{M1_(0,0),b1},{M1_(1,0),b2},{M1_(2,0),b3},{M1_(3,0),1}};
	eq14 := matrix{{M1_(0,1),b1},{M1_(1,1),b2},{M1_(2,1),b3},{M1_(3,1),1}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq13))|(flatten entries gens minors(2,eq14));
-- Flag eqn x2 \subset E23
	eq31 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},c1,c2,c3},{x_{2,3,0},0,0,1}};
	eq32 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},c1,c2,c3},{x_{2,3,1},0,0,1}};
	eq33 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},c1,c2,c3},{x_{2,3,2},0,0,1}};
	eq34 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},c1,c2,c3},{x_{2,3,3},0,0,1}};
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
	M9 := matrix{{b1,1,0,0},{b2,0,1,0},{b3,c1,c2,c3},{1,0,0,1}};
	equationsCover = append(equationsCover, det M9);
	M7 := matrix{{1,1,0,0},{0,0,1,0},{d1,0,0,1},{d3,e1,e2,e3}};
	M8 := matrix{{0,1,0,0},{1,0,1,0},{d2,0,0,1},{d4,e1,e2,e3}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = append(equationsCover, det M8);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cm;
	F = sub(F, bigRing);
-- Putting all of the equations together
	idealCmCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCmCk
	
)

-- (4) JacCm4(): uses affine chart containing ({[a:1]},{[c:1]}).

JacCm4 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cm = new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{0,0},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,d4,e1,e2,e3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{b3},{1}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{c1,c2,c3},{0,0,1}};
	E23 = sub(E23, bigRing);
	E32 := matrix{{1,0},{0,1},{d1,d2},{d3,d4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{e1,e2,e3},{0,0,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
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
	eq13 := matrix{{M1_(0,0),b1},{M1_(1,0),b2},{M1_(2,0),b3},{M1_(3,0),1}};
	eq14 := matrix{{M1_(0,1),b1},{M1_(1,1),b2},{M1_(2,1),b3},{M1_(3,1),1}};
	equationsCover = equationsCover|(flatten entries gens minors(2,eq13))|(flatten entries gens minors(2,eq14));
-- Flag eqn x2 \subset E23
	eq31 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},c1,c2,c3},{x_{2,3,0},0,0,1}};
	eq32 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},c1,c2,c3},{x_{2,3,1},0,0,1}};
	eq33 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},c1,c2,c3},{x_{2,3,2},0,0,1}};
	eq34 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},c1,c2,c3},{x_{2,3,3},0,0,1}};
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
	eq1  := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},e1,e2,e3},{x_{1,3,0},0,0,1}};
	eq2  := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},e1,e2,e3},{x_{1,3,1},0,0,1}};
	eq3  := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},e1,e2,e3},{x_{1,3,2},0,0,1}};
	eq4  := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},e1,e2,e3},{x_{1,3,3},0,0,1}};
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
	M9 := matrix{{b1,1,0,0},{b2,0,1,0},{b3,c1,c2,c3},{1,0,0,1}};
	equationsCover = append(equationsCover, det M9);
	M7 := matrix{{1,1,0,0},{0,0,1,0},{d1,e1,e2,e3},{d3,0,0,1}};
	M8 := matrix{{0,1,0,0},{1,0,1,0},{d2,e1,e2,e3},{d4,0,0,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = append(equationsCover, det M8);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cm;
	F = sub(F, bigRing);
-- Putting all of the equations together
	idealCmCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCmCk
)

subJacCm = method()
subJacCm(Matrix) := (M) -> (
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

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-- Ev computation for Cpsi
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

-- (1) JacCpsi():  uses affine chart containing  ({[1:b]}, {[1:d]} , {[1:h]}, {1:alpha]}).

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
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi = method()
subJacCpsi(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ d1=>b1, b3=>c3*b2, g=>f3});
	M
)

-- (2) JacCpsi2(): uses affine chart containing ({[1:b]}, {[1:d]}, {[1:h]}, {[beta:1]}).

JacCpsi2 = () -> (
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
-- Putting all of the equations together
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi2 = method()
subJacCpsi2(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ b2=>b1*d1, b3=>c3, g=>f3});
	M
)


-- (3) JacCpsi3(): uses affine chart containing ({[1:b]}, {[1:d]}, {[e:1]}, {[1:alpha]}).

JacCpsi3 = () -> (
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
	E33 := matrix{{1,0,0},{0,1,0},{f1,f2,f3},{0,0,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
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
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},f1,f2,f3},{x_{1,3,0},0,0,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},f1,f2,f3},{x_{1,3,1},0,0,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},f1,f2,f3},{x_{1,3,2},0,0,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},f1,f2,f3},{x_{1,3,3},0,0,1}};
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
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,f1,f2,f3},{e3,0,0,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,f1,f2,f3},{e4,0,0,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi3 = method()
subJacCpsi3(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ b3=>b2*c3, d1=>b1, g=>f3});
	M
)

-- (4) JacCpsi4(): uses affine chart containing ({[1:b]}, {[1:d]}, {[e:1]}, {[beta:1]}).

JacCpsi4 = () -> (
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
	E33 := matrix{{1,0,0},{0,1,0},{f1,f2,f3},{0,0,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
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
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},f1,f2,f3},{x_{1,3,0},0,0,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},f1,f2,f3},{x_{1,3,1},0,0,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},f1,f2,f3},{x_{1,3,2},0,0,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},f1,f2,f3},{x_{1,3,3},0,0,1}};
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
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,f1,f2,f3},{e3,0,0,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,f1,f2,f3},{e4,0,0,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi4 = method()
subJacCpsi4(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ b2=>b1*d1, b3=>c3, g=>f3});
	M
)

-- (5) JacCpsi5(): uses affine chart containing  ({[1:b]}, {[c:1]}, {[1:h]}, {[1:alpha]}).

JacCpsi5 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{1},{b1},{b2},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{c1,c2,c3},{0,0,1}};
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
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},c1,c2,c3},{x_{2,3,0},0,0,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},c1,c2,c3},{x_{2,3,1},0,0,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},c1,c2,c3},{x_{2,3,2},0,0,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},c1,c2,c3},{x_{2,3,3},0,0,1}};
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
	M7 := matrix{{1,1,0,0},{b1,0,1,0},{b2,c1,c2,c3},{b3,0,0,1}};
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
-- Putting all of the equations together
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi5 = method()
subJacCpsi5(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ b2=>b3*c3, d1=>b1, g=>f3});
	M
)

-- (6) JacCpsi6(): uses affine chart containing ({[1:b]}, {[c:1]}, {[1:h]}, {[beta:1]}).

JacCpsi6 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{b3},{1}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{c1,c2,c3},{0,0,1}};
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
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},c1,c2,c3},{x_{2,3,0},0,0,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},c1,c2,c3},{x_{2,3,1},0,0,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},c1,c2,c3},{x_{2,3,2},0,0,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},c1,c2,c3},{x_{2,3,3},0,0,1}};
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
	M7 := matrix{{b1,1,0,0},{b2,0,1,0},{b3,c1,c2,c3},{1,0,0,1}};
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
-- Putting all of the equations together
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi6 = method()
subJacCpsi6(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ b2=>b1*d1, b3=>c3, g=>f3});
	M
)

-- (7) JacCpsi7(): uses affine chart containing ({[1:b]}, {[c:1]}, {[e:1]}, {[1:alpha]}).

JacCpsi7 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{1},{b1},{b2},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{c1,c2,c3},{0,0,1}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{1},{d1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{f1,f2,f3},{0,0,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
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
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},c1,c2,c3},{x_{2,3,0},0,0,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},c1,c2,c3},{x_{2,3,1},0,0,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},c1,c2,c3},{x_{2,3,2},0,0,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},c1,c2,c3},{x_{2,3,3},0,0,1}};
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
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},f1,f2,f3},{x_{1,3,0},0,0,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},f1,f2,f3},{x_{1,3,1},0,0,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},f1,f2,f3},{x_{1,3,2},0,0,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},f1,f2,f3},{x_{1,3,3},0,0,1}};
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
	M7 := matrix{{1,1,0,0},{b1,0,1,0},{b2,c1,c2,c3},{b3,0,0,1}};
	M8 := matrix{{1,1,0},{d1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,f1,f2,f3},{e3,0,0,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,f1,f2,f3},{e4,0,0,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi7 = method()
subJacCpsi7(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ d1=>b1, b2=>b3*c3, g=>f3});
	M
)

-- (8) JacCpsi8(): uses affine chart containing ({[1:b]}, {[c:1]}, {[e:1]}, {[beta:1]}).

JacCpsi8 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{b3},{1}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{c1,c2,c3},{0,0,1}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{1},{d1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{f1,f2,f3},{0,0,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
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
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},c1,c2,c3},{x_{2,3,0},0,0,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},c1,c2,c3},{x_{2,3,1},0,0,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},c1,c2,c3},{x_{2,3,2},0,0,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},c1,c2,c3},{x_{2,3,3},0,0,1}};
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
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},f1,f2,f3},{x_{1,3,0},0,0,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},f1,f2,f3},{x_{1,3,1},0,0,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},f1,f2,f3},{x_{1,3,2},0,0,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},f1,f2,f3},{x_{1,3,3},0,0,1}};
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
	M7 := matrix{{b1,1,0,0},{b2,0,1,0},{b3,c1,c2,c3},{1,0,0,1}};
	M8 := matrix{{1,1,0},{d1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,f1,f2,f3},{e3,0,0,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,f1,f2,f3},{e4,0,0,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi8 = method()
subJacCpsi8(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ b2=>b1*d1, b3=>c3, g=>f3});
	M
)

-- (9) JacCpsi9(): uses affine chart containing ({[a:1]}, {[1:d]}, {[1:h]}, {[1,alpha]}).

JacCpsi9 = () -> (
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
-- Putting all of the equations together
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi9 = method()
subJacCpsi9(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ b3=>b2*c3, d1=>b1, g=>f3});
	M
)


-- (10) JacCpsi10(): uses affine chart containing ({[a:1]}, {[1:d]}, {[1:h]}, {[beta:1]}).

JacCpsi10 = () -> (
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
-- Putting all of the equations together
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi10 = method()
subJacCpsi10(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ b1=>b2*d1, b3=>c3, g=>f3});
	M
)

-- (11) JacCpsi11(): uses affine chart containing ({[a:1]}, {[1:d]}, {[e:1]}, {[1:alpha]}).

JacCpsi11 = () -> (
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
	E33 := matrix{{1,0,0},{0,1,0},{f1,f2,f3},{0,0,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
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
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},f1,f2,f3},{x_{1,3,0},0,0,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},f1,f2,f3},{x_{1,3,1},0,0,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},f1,f2,f3},{x_{1,3,2},0,0,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},f1,f2,f3},{x_{1,3,3},0,0,1}};
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
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,f1,f2,f3},{e3,0,0,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,f1,f2,f3},{e4,0,0,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi11 = method()
subJacCpsi11(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ d1=>b1, b3=>c3*b2, g=>f3});
	M
)

-- (12) JacCpsi12(): uses affine chart containing ({[a:1]}, {[1:d]}, {[e:1]}, {[beta:1]}).

JacCpsi12 = () -> (
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
	E33 := matrix{{1,0,0},{0,1,0},{f1,f2,f3},{0,0,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
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
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},f1,f2,f3},{x_{1,3,0},0,0,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},f1,f2,f3},{x_{1,3,1},0,0,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},f1,f2,f3},{x_{1,3,2},0,0,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},f1,f2,f3},{x_{1,3,3},0,0,1}};
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
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,f1,f2,f3},{e3,0,0,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,f1,f2,f3},{e4,0,0,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi12 = method()
subJacCpsi12(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ b1=>b2*d1, b3=>c3, g=>f3});
	M
)

-- (13) JacCpsi13(): uses affine chart containing ({[a:1]}, {[c:1]}, {[1:h]}, {[1:alpha]}).

JacCpsi13 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{1},{b2},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{c1,c2,c3},{0,0,1}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{d1},{1},{d2},{d3}};
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
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},c1,c2,c3},{x_{2,3,0},0,0,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},c1,c2,c3},{x_{2,3,1},0,0,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},c1,c2,c3},{x_{2,3,2},0,0,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},c1,c2,c3},{x_{2,3,3},0,0,1}};
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
	M7 := matrix{{b1,1,0,0},{1,0,1,0},{b2,c1,c2,c3},{b3,0,0,1}};
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
-- Putting all of the equations together
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi13 = method()
subJacCpsi13(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ d1=>b1, b2=>b3*c3, g=>f3});
	M
)

-- (14) JacCpsi14(): uses affine chart containing ({[a:1]}, {[c:1]}, {[1:h]}, {[1:alpha]}).

JacCpsi14 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{b3},{1}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{c1,c2,c3},{0,0,1}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{d1},{1},{d2},{d3}};
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
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},c1,c2,c3},{x_{2,3,0},0,0,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},c1,c2,c3},{x_{2,3,1},0,0,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},c1,c2,c3},{x_{2,3,2},0,0,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},c1,c2,c3},{x_{2,3,3},0,0,1}};
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
	M7 := matrix{{b1,1,0,0},{b2,0,1,0},{b3,c1,c2,c3},{1,0,0,1}};
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

subJacCpsi14 = method()
subJacCpsi14(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ b1=>b2*d1, b3=>c3, g=>f3});
	M
)

-- (15) JacCpsi15(): uses affine chart containing ({[a:1]}, {[c:1]}, {[e:1]}, {[1:alpha]}).

JacCpsi15 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{1},{b2},{b3}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{c1,c2,c3},{0,0,1}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{d1},{1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{f1,f2,f3},{0,0,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
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
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},c1,c2,c3},{x_{2,3,0},0,0,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},c1,c2,c3},{x_{2,3,1},0,0,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},c1,c2,c3},{x_{2,3,2},0,0,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},c1,c2,c3},{x_{2,3,3},0,0,1}};
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
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},f1,f2,f3},{x_{1,3,0},0,0,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},f1,f2,f3},{x_{1,3,1},0,0,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},f1,f2,f3},{x_{1,3,2},0,0,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},f1,f2,f3},{x_{1,3,3},0,0,1}};
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
	M7 := matrix{{b1,1,0,0},{1,0,1,0},{b2,c1,c2,c3},{b3,0,0,1}};
	M8 := matrix{{d1,1,0},{1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,f1,f2,f3},{e3,0,0,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,f1,f2,f3},{e4,0,0,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi15 = method()
subJacCpsi15(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ d1=>b1, b2=>b3*c3, g=>f3});
	M
)

-- (16) JacCpsi16(): uses affine chart containing ({[a:1]}, {[c:1]}, {[e:1]}, {[beta:1]}).

JacCpsi16 = () -> (
	Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
	Cpsi := new RankConditions from ({2,4,4,4,2},{{2,3,3,2},{1,2,1},{1,1},{0}});
	ringE := QQ[a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3,e4,f1,f2,f3,g];
	bigRing = (transposeVoganV Ck) ** ringE ** (ring x0Matrix()) ** (ring x1Matrix()) ** (ring x2Matrix()) ** (ring x3Matrix());
	equationsCk := sub(getTransposeEquations(computeDual Ck), bigRing);
 	E12 := matrix{{1,0},{0,1},{a1,a2},{a3,a4}};
	E12 = sub(E12, bigRing);
	E21 := matrix{{b1},{b2},{b3},{1}};
	E21 = sub(E21, bigRing);
	E23 := matrix{{1,0,0},{0,1,0},{c1,c2,c3},{0,0,1}};
	E23 = sub(E23, bigRing);
	E31 := matrix{{d1},{1},{d2},{d3}};
	E31 = sub( E31, bigRing);
	E32 := matrix{{1,0},{0,1},{e1,e2},{e3,e4}};
	E32 = sub(E32, bigRing);
	E33 := matrix{{1,0,0},{0,1,0},{f1,f2,f3},{0,0,1}};
	E33 = sub(E33, bigRing);
	E41 := matrix{{g},{1}};
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
	eq21 := matrix{{x_{2,0,0},1,0,0},{x_{2,1,0},0,1,0},{x_{2,2,0},c1,c2,c3},{x_{2,3,0},0,0,1}};
	eq22 := matrix{{x_{2,0,1},1,0,0},{x_{2,1,1},0,1,0},{x_{2,2,1},c1,c2,c3},{x_{2,3,1},0,0,1}};
	eq23 := matrix{{x_{2,0,2},1,0,0},{x_{2,1,2},0,1,0},{x_{2,2,2},c1,c2,c3},{x_{2,3,2},0,0,1}};
	eq24 := matrix{{x_{2,0,3},1,0,0},{x_{2,1,3},0,1,0},{x_{2,2,3},c1,c2,c3},{x_{2,3,3},0,0,1}};
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
	eq31 := matrix{{x_{1,0,0},1,0,0},{x_{1,1,0},0,1,0},{x_{1,2,0},f1,f2,f3},{x_{1,3,0},0,0,1}};
	eq32 := matrix{{x_{1,0,1},1,0,0},{x_{1,1,1},0,1,0},{x_{1,2,1},f1,f2,f3},{x_{1,3,1},0,0,1}};
	eq33 := matrix{{x_{1,0,2},1,0,0},{x_{1,1,2},0,1,0},{x_{1,2,2},f1,f2,f3},{x_{1,3,2},0,0,1}};
	eq34 := matrix{{x_{1,0,3},1,0,0},{x_{1,1,3},0,1,0},{x_{1,2,3},f1,f2,f3},{x_{1,3,3},0,0,1}};
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
	M7 := matrix{{b1,1,0,0},{b2,0,1,0},{b3,c1,c2,c3},{1,0,0,1}};
	M8 := matrix{{d1,1,0},{1,0,1},{d2,e1,e2},{d3,e3,e4}};
	M9 := matrix{{1,1,0,0},{0,0,1,0},{e1,f1,f2,f3},{e3,0,0,1}};
	M10:= matrix{{0,1,0,0},{1,0,1,0},{e2,f1,f2,f3},{e4,0,0,1}};
	equationsCover = append(equationsCover, det M7);
	equationsCover = equationsCover|(flatten entries gens minors(3,M8));
	equationsCover = append(equationsCover, det M9);
	equationsCover = append(equationsCover, det M10);
	equationsCover = ideal equationsCover;
-- Killing form for vogan V
	F := buildF Cpsi;
	F = sub(F, bigRing);
-- Putting all of the equations together
	idealCpsiCk := (ideal F) + equationsCk + equationsCover;
	jacobian idealCpsiCk
)

subJacCpsi16 = method()
subJacCpsi16(Matrix) := (M) -> (
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
	M = substitute(M,{ c1=>0,c2=>0,d2=>0,d3=>0});
	M = substitute(M,{ e1=>0, e2=>0, e3=>0, e4=>0});
	M = substitute(M,{ f1=>0, f2=>0});
	M = substitute(M,{ b1=>b2*d1, b3=>c3, g=>f3});
	M
)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------







