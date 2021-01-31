----------------------------------------------------------------------------------------------------------------
-- File: VanishingCycles.m2
-- Author: Reginald Lybbert
--
-- This file provides methods that may be used during the computation of 
-- the Ev functor.  Currently, the full computation of Ev is unfinished,
-- as such, more functionality may be added to this file in the future.
--
-- Currently, the following functions are implemented
--
-- buildF RankConditions
-- jacobianForC'DualTimesCCover
--
-----------------------------------------------------------
-- buildF RankConditions
-- 
-- This function provides the function f used in the computation
-- of Ev.  To be precise, f:V x V* -> K is defined by 
-- f(X,Y) = tr(XY) = tr(YX)
--
-- Example:
--
-- i1 : orbit = new RankConditions from ({1,3,1},{{1,1},{0}})
--
-- o1 = 1   3   1
--        1   1
--          0
--
-- o1 : RankConditions
--
-- i2 : buildF orbit
--
-- o2 = x         y          + x         y          + x         y          + x         y          + x         y          + x         y          
--       {0, 0, 0} {0, 0, 0}    {0, 0, 1} {0, 1, 0}    {0, 0, 2} {0, 2, 0}    {1, 0, 0} {1, 0, 0}    {1, 1, 0} {1, 0, 1}    {1, 2, 0} {1, 0, 2}
--
-- o2 : QQ[x         , x         , x         , x         , x         , x         , y         , y         , y         , y         , y         , y         ]
--          {0, 0, 0}   {0, 0, 1}   {0, 0, 2}   {1, 0, 0}   {1, 1, 0}   {1, 2, 0}   {0, 0, 0}   {0, 1, 0}   {0, 2, 0}   {1, 0, 0}   {1, 0, 1}   {1, 0, 2}
-----------------------------------------------------------
-- jacobianForC'DualTimesCCover
--
-- This was intended as a helper function, but this is all that has
-- currently been implemented.  Given two orbits C' and C, this function
-- returns the jacobian of the equations for f, the cover  of C, C itself
-- and C'*.  The future goal will be the find the singularity locus of this
-- Jacobian on the space ~C x C'*.
--
-- Example
--
-- i3 : suborbit = new RankConditions from ({1,3,1},{{0,1},{0}})
--
-- o3 = 1   3   1
--        0   1
--          0
--
-- o3 : RankConditions
--
-- i4 : jacobianForC'DualTimesCCover (suborbit, orbit)
--
-- o4 = {1, 0, 0} | y_{0, 0, 0} x_{1, 0, 0} 0 0 0 0        0           0           0            1           |
--      {1, 0, 0} | y_{0, 1, 0} x_{1, 1, 0} 0 0 0 0        0           0           0            cc_{1, 0}   |
--      {1, 0, 0} | y_{0, 2, 0} x_{1, 2, 0} 0 0 0 0        0           0           0            cc_{1, 1}   |
--      {1, 0, 0} | y_{1, 0, 0} x_{0, 0, 0} 0 0 0 0        cc_{1,0}    c_{1,1}     0            0           |
--      {1, 0, 0} | y_{1, 0, 1} x_{0, 0, 1} 0 0 0 0        -1          0           cc_{1,1}     0           |
--      {1, 0, 0} | y_{1, 0, 2} x_{0, 0, 2} 0 0 0 0        0           -1          -cc_{1,0}    0           |
--      {0, 1, 0} | x_{0, 0, 0} 0           0 0 0 0        0           0           0            0           |
--      {0, 1, 0} | x_{0, 0, 1} 0           0 0 0 0        0           0           0            0           |
--      {0, 1, 0} | x_{0, 0, 2} 0           0 0 0 0        0           0           0            0           |
--      {0, 1, 0} | x_{1, 0, 0} 0           0 0 1 0        0           0           0            0           |
--      {0, 1, 0} | x_{1, 1, 0} 0           0 1 0 0        0           0           0            0           |
--      {0, 1, 0} | x_{1, 2, 0} 0           1 0 0 0        0           0           0            0           |
--      {0, 0, 1} | 0           0           0 0 0 cc_{1,3} x_{1, 0, 0} 0           -x_{1, 2, 0} x_{0, 0, 1} |
--      {0, 0, 1} | 0           0           0 0 0 -1       0           x_{1, 0, 0} x_{1, 1, 0}  x_{0, 0, 2} |
--      {0, 0, 1} | 0           0           0 0 0 1        0           0           0            0           |
--      {0, 0, 1} | 0           0           0 0 0 cc_{1,0} 0           0           0            0           |
--
--                                                                                                                                                                16                                                                                                                                                         10           
-- o4 : Matrix QQ[x         , x         , x         , x         , x         , x         , y         , y         , y         , y         , y         , y         ])   <--- (QQ[x         , x         , x         , x         , x         , x         , y         , y         , y         , y         , y         , y         ])
--                 {0, 0, 0}   {0, 0, 1}   {0, 0, 2}   {1, 0, 0}   {1, 1, 0}   {1, 2, 0}   {0, 0, 0}   {0, 1, 0}   {0, 2, 0}   {1, 0, 0}   {1, 0, 1}   {1, 0, 2}               {0, 0, 0}   {0, 0, 1}   {0, 0, 2}   {1, 0, 0}   {1, 1, 0}   {1, 2, 0}   {0, 0, 0}   {0, 1, 0}   {0, 2, 0}   {1, 0, 0}   {1, 0, 1}   {1, 0, 2}
--
----------------------------------------------------------------------------------------------------------------


needs "ComputeCover.m2"
needs "ComputeDuals.m2"

buildF = method()
buildF RankConditions := (R) -> (
   containingRing := (voganV R) ** (transposeVoganV R);
   use containingRing;
   k := #R.dims - 1;
   matrices := {};
   for i in 0..k-1 do (
       matrices = append(matrices, (map(containingRing^(R.dims#i),R.dims#(i+1),(m,n)->x_{i,m,n}),map(containingRing^(R.dims#(i+1)),R.dims#i,(m,n)->y_{i,m,n})))
   );
   sum apply(matrices, p -> trace(p#0*p#1))  
)

jacobianForC'DualTimesCCover = method()
jacobianForC'DualTimesCCover (RankConditions, RankConditions) := (C', C) -> (
    coverRing := getCoverRing C;
    bigRing := (voganV C) ** (transposeVoganV C') ** coverRing;
    generatingIdeal := sub(ideal buildF C,bigRing);
    generatingIdeal = generatingIdeal + sub(getEquations C,bigRing);
    generatingIdeal = generatingIdeal + sub(getTransposeEquations (computeDual C'),bigRing);
    generatingIdeal = generatingIdeal + sub(computeCover C,bigRing);
    jacobian generatingIdeal
    --generatingIdeal
)



