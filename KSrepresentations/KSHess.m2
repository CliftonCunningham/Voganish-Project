-- This code is designed to verify the Morse condition at "a smooth point"
-- Basically run the following in order:
--   load "KSHess.m2"
--   eqconf = EqCpsi();
--   locvarinf = getLocalCoordList( eqconf#1, subXCpsi );
--   impvarpart = getImpVarPartials ( eqconf#1 , locvarinf, subXCpsi );
--   hess = getHessian(eqconf#0, locvarinf, impvarpart);
--   rnk = rank( subXCpsi( hess ) );
-- You want the rank to be the codimension of the singular locus in order to satisfy the Hessian condition.
--        (Dim C_{KS} + Dim C_{Psi} - Dim V)
-- The following attempts to verify if the Hessian determinant is a square
--   subhess = getSubHessian( hess, subXCpsi );
--   iso=findIso(subhess#0, subXCpsi);
-- You want the numcolumns iso#0 to be half of the rank of the Hessian
--    NOTE: if this does not happen you can basically conclude nothing. This checks a sufficient but not necissary condition.
--

-- TODO - add a function EqCks which sets up equations for Cks
--        add a function SubXCks which substitutes a vaguely generic point

needs "ComputeDuals.m2" 
needs "NetworkFlow.m2"
needs "ComputeCover.m2" 
needs "VanishingCycles.m2"
needs "VoganV.m2"
needs "PSNF.m2"
needs "KSrepresentations.m2"

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------


-- Constructions the Equations relevant for computing Ev_{C_KS}(IC(1_{C_\psi}))
EqCpsi = () -> (
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
-- Putting all of the equations together!
        {F,equationsCk + equationsCover}
)

-- evaluates M at a point in the singular locus of f relevant to computing Ev_{C_KS}(IC(1_{C_\psi}))
--   The choice of point is a bit arbitrary, though this feels sketchy if everything "runs to completion" its fine
--   as the Hessian computations are valid on an open nhd of the point so are valid at a regular point
--   (Note: If this point wasnt regular we could possibly get a rank that is too low to get desired conclusion)
-- TODO: This function is hardcoded in a way which is strange
--       It might be better to create a RingMap which encodes this and pass it as a parameter to later functions
subXCpsi = (M) -> (
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
-- THIS IS THE FIRST POINT
        M = substitute(M,{ b1=>0, f3=>0, c3=>0 } );
-- THis is a choice of t1...
        M = substitute(M,{ b2=>1, y_{2,0,2}=>-1,  y_{2,1,3}=>7 } );

-- NOW DONE
        M
);


--- Gets equations for Ck times Ck*
EqCk = () -> (
        Ck = new RankConditions from ({2,4,4,4,2},{{2,2,2,2},{0,2,0},{0,0},{0}});
        bigRing = (voganV Ck) ** (transposeVoganV Ck);
        equationsCk :=  sub(getEquations(Ck), bigRing);
        equationsCkDual := sub(getTransposeEquations(computeDual Ck), bigRing);

        F := buildF Ck;
        F = sub(F, bigRing);

        {F, equationsCk + equationsCkDual}
);

subXCks = (M) -> (
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
-- This is a choice of t1...
        M = substitute(M,{ y_{2,0,2}=>-1,  y_{2,1,3}=>7 } );

-- NOW DONE
        M
);


-- Computes a set of local coordinates using "jacobian" condition
--   Finds a minimal list "LC" such that Jac( LC + eqns ) has maximal rank, which will be equal to dim(R) assuming the point described by subXCpsi
--          Is smooth. The list LC describes functions whose image in the completed local ring in a neighbourhood of subXCpsi point are local coordinates.
-- NOTE: This is not optimally implemented
--        We recompute the jacobian and redo all the substitutions for every loop iteration
--        This is unneccissary, we could compute one jacobian and focus on selecting equations by 
--        "adding/deleting rows"
getLocalCoordList = ( eqns, ptfunc ) -> (
   R = ring eqns;
   use R;
   jac = jacobian ( eqns );
   evjac = ptfunc( jac );
   oldrank = rank( evjac );
   varlist = gens R;
   lclist = new List ;
   lclistindex = new List ;
   nlclist = new List ;
   nlclistindex = new List ;
   ind=0;
   -- For each variable check if it increases the rank or not and add to the appropriate list
   for var in varlist do (
      i1 = ideal(lclist);
      i2 = ideal(var);
      jac = jacobian ( i1+i2+eqns );
      evjac = ptfunc( jac );
      newrank = rank( evjac );
      if newrank > oldrank then (
          lclist = append(lclist, var);
          lclistindex = append(lclistindex, ind);
      ) else (
          nlclist = append(nlclist, var);
          nlclistindex = append(nlclistindex, ind);
      );  
      oldrank = newrank;
      ind = ind+1;
   );
   { lclist, nlclist, lclistindex, nlclistindex }
)

-- returns primarily a matrix such that
--    the i,j entry describes the partial derivative of the -- implicit function by the -- local variable
--    where ... those are totally i and j
-- NOTE - this code is not optimal
--        We resubstitute values on each loop iteration while selecting the minimal system
getImpVarPartials = ( eqns , locvarinf, ptfunc ) -> (
    R = ring eqns;
    use R;
    jac = jacobian( eqns );
    evjac = ptfunc( jac );

    -- smat is the portion of the jacobian matrix which contains the derivatives with respect
    --   to the implicit functions
    --   with respect one of the actual local variables.
    --   
    --   We had g_k(x_i, y_j(x_i)) and we want to solve for
    --      del y_j/ del x_i  for all j i.
    --   we do this "one i at a time" (not really, but pretend)
    --   to compute del g_k / del x_i you should apply chain rule to the y_j coords
    --        in addition to including the "del g_k/del x_i" from the x_i coord
    --   the k-th column of smat contains the partial derivatives del g_k/dey y_j
    --   so you should imagine we are setting up the matrix equation
    --   using del g_k / del x_i = 0 
    --    - "del g_k/del x_i"  =  ( del y_j/ del x_i ) * smat
    --   wheere ( del y_j/ del x_i ) is a row vector

    smat = submatrix( evjac, locvarinf#3, );
    
    eqset = new List;
    oldrank = 0;

    -- we are trying to locally invert smat
    --   We want to select enough equations so that we can solve for each
    --   (Their are extra equations so we can drop some of the g_k
    --   
    for i in 0..((numgens(eqns))-1) do (
       tryeqset = append( eqset, i );
       trymat = submatrix( smat,, tryeqset );
       newrank = rank( trymat );
       if newrank > oldrank then (
          eqset = tryeqset;
          oldrank = newrank;
       );
    );

    invertablemat = submatrix( jac, locvarinf#3, eqset );
 
    --  each row of partsub corresponds to the vector  "del g_k/del x_i" for a fixed i, as we run over k.
    partsub = submatrix( jac, locvarinf#2, eqset );

    --  We know that trymat is locally invertible at the point of interest, so though we are working
    --  in the fraction field we know the inverse lifts to the local ring at the point.
    F = frac( R );
    tmi = (sub(invertablemat,F))^(-1);

    --   We have solved
    --       - "del g_k/del x_i"  =  ( del y_j/ del x_i ) * smat
    --   so that each row gives the solution for a different i
    --       ( del y_j/ del x_i )  = -  "del g_k/del x_i" * smat^-1
    --   where  "del g_k/del x_i" was the partial with respect to the x_i coord (knowing that the total partial wrt x_i is 0)
    implicitpart = -partsub * tmi;
    implicitpart
)
 
-- This function implements the quotient rule to take the jacobian of the function f
--   (provided f is a rational function) otherwise it just uses standard Jacobian Function.
jacobianWithQuotientRule = ( f ) -> (
    R = ring f;
    j = 0;
    if instance(R,FractionField) then (
      num = numerator f;
      den = denominator f;

      jn = jacobian( ideal num );
      jd = jacobian( ideal den );
    
      j = sub((jn*den - jd*num),ring f)*(1/den^2);
    ) else (
      j = jacobian( ideal(f) );
    );

    j
)

-- returns the "jacobian" of partial derivatives of f with respect to the local variables
--     it uses the precomputed partials of the implicit functions with respect to the local variables
getPartial = (f,locvarinf,impvarpart) -> (
   jacf =  jacobianWithQuotientRule ( f );
   jacfloc = submatrix( jacf, locvarinf#2, );
   jacfnloc = submatrix( jacf, locvarinf#3, );

   partiallist = new List;

   numlocvar = #(locvarinf#2);

   -- I think there is a "Matrix way" to implement this loop, probably does not matter for performance
   for i in 0..(numlocvar-1) do (
        resi = submatrix( impvarpart, {i} , );

        -- This should be the contribution from the chain rule
        --   The ith row of  impvarpart should contain
        --        ( del y_j/ del x_i ) 
        --   The column jacfnloc should contain
        --        ( del f/del y_j)
        chain = resi * jacfnloc;
        
        partiallist = append(partiallist, chain_(0,0) + jacfloc_(i,0)) 
   );
   partiallist
)


-- returns the "Hessian" of the partial deriviatives of f with respect to the local variables.
getHessian = (f, locvarinf, impvarpart) -> (

   -- List of first order partial derivatives by the local variables
   firstpartials = getPartial(f,locvarinf,impvarpart);

   hesslist = new List;
   for p in firstpartials do (
      hesslist = append(hesslist, getPartial(p,locvarinf,impvarpart) );
   );
   hess = matrix( hesslist );

   hess
);



-- In principal this is supposed to find a maximal rank sub-block of the Hessian
--   We try to delete rows and columns that do not lower rank
getSubHessian = (hess, ptfunc) -> (
   locvset = new List;
   entryhess = ptfunc( hess );
   oldrank = rank(  entryhess );
   
   for i in 0..(numgens target hess - 1) do (
      trylocvset = append( locvset, i );
      trysubhess = submatrix'( entryhess, trylocvset, trylocvset );
      newrank = rank( trysubhess );
      if newrank == oldrank then (
         locvset = trylocvset;
      );
   );

   subhess = submatrix'( hess, locvset, locvset );
   subhessent = submatrix'( entryhess, locvset, locvset );
   { subhess, subhessent, locvset }
)


-- Tries to find a maximal isotropic subspace
--   Does not actually try that hard, but is relatively likely to succeed.
findIso = ( hess, ptfunc ) -> (

   elthess =  ptfunc( matrix( hess ) );

   n=numColumns hess;
   yesList = new List;
   noList = new List;

   for i in 0..(n-1) do (
      tryList = append( yesList, i );
      sh = submatrix( hess, tryList, tryList );
      if sh == 0 then (
         yesList = append( yesList, i );
      ) else ( 
         noList = append( noList, i );
      );
   );

   oldrank = 0;
   compList = new List;
 
   zerList = yesList;
   zer = submatrix( hess, yesList, yesList );

   for i in noList do (
      tryList = append( yesList, i );
      sh = submatrix( elthess, tryList, tryList );
      if rank(sh) > oldrank then (
         yesList = append( yesList, i );
         oldrank=rank(sh);
      ) else (
         compList = append( compList, i ); 
      );
   );


   hyper = submatrix( hess, yesList, yesList );
   comp = submatrix( hess, compList, compList );
   B = submatrix( hess, yesList, compList );

   return({ zer, hyper, comp, B, zerList});
)

