----------------------------------------------------------------------------------------------------------------
-- File: ComputeDuals.m2
-- Author: Reginald Lybbert
--
-- This file provides the functions related to the computation of dual 
-- orbits in a Vogan variety.  The following functions are intended for
-- use.
-- 
--
-- computeDualForm RankConditions
-- computeDual RankConditions
-- findPathologies List
--
--
-- Explanations for each function are given below 
--
-- The code for these functions, along with some helper functions is
-- found at the end of the file.  (Also, there may be some other useful
-- stuff hidden in the code.)
--
-----------------------------------------------------------
--computeDualForm RankConditions
--
-- This function takes a provided orbit, and produces a list of
-- lists, which can be seen as matrices, of 1s and 0s such that
-- any list of matrices Y which have 0s in the same place will 
-- satisfy the condition [X,Y] = 0, where X is the representative
-- of the provided orbit given by getRep2
--
-- While primarily a helper function for computeDual, the output
-- of this function may be useful.
--
-- Example:
-- 
-- i1 : orbit = new RankConditions from ({1,2,2,1},{{1,2,1},{1,1},{0}})
--
-- o1 = 1   2   2   1
--        1   2   1
--          1   1
--            0
--
-- o1 : RankConditions
--
-- i2 : computeDualForm orbit
--
-- o2 = {{{1}, {0}}, {{0, 1}, {0, 0}}, {{0, 1}}}
--
-- o2 : List
-----------------------------------------------------------
-- computeDual RankConditions
--
-- This function takes in an orbit, and returns the dual.
-- It uses the algorithm described in my earlier document
-- regarding dual compuations.
--
-- Example:
--
-- i3 : computeDual orbit
-- 
-- o3 = 1   2   2   1
--        1   1   1
--          0   0
--            0
-- 
-- o3 : RankConditions
--
-----------------------------------------------------------
-- findPathologies List
--
-- Given a list of dimensions that defines a Vogan variety V,
-- this method returns a list of all pairs of orbits (C', C)
-- such that C' < C and C'* < C*
--
-- Example:
--
-- i4 : findPathologies ({1,3,1})
--
-- o4 = {(1   3   1, 1   3   1), (1   3   1, 1   3   1)} 
--          0   1      1   1        1   0      1   1
--            0          0            0          0
--
-- o4 : List
--
-------------------------compucomp---------------------------------------------------------------------------------------


needs "VoganV.m2"
needs "NetworkFlow.m2"

isDoublyNestedRectangular = method()
isDoublyNestedRectangular List := (ls) -> (
      if #ls == 0 then false
      else if not instance(ls#0,List) then false
      else ( 
         n := #(ls#0);
         result := true;
         for x in ls do ( 
             result = result and instance(x,List) and #x == n;
         );
         result
      )
)

MatrixForm = new Type of MutableList
new MatrixForm from List := (MatrixForm, ls) -> ( if #ls == 0 then new MutableList from {}
                                                  else if isDoublyNestedRectangular(ls) then new MutableList from apply(ls, x -> new MutableList from x)
                                                  else error "expected rectangular doubly nested list" )

net MatrixForm := (m) -> net toList(apply(m,toList))

MatrixForm * MatrixForm := (m,n) -> (
    result := new MutableList from {};
    if #(m#0) == #n then (
        for i in 0..#m-1 do (
            result = append(result,new MutableList from {});
            for j in 0..#(n#0)-1 do (
                entry := 0;
                for k in 0..#n-1 do (
                     if (m#i#k != 0 and n#k#j != 0) then entry = 1; 
                );
                result#i = append(result#i, entry);
            );
        );
    );
    new MatrixForm from result
)

findRow = method()
findRow (List,ZZ) := (s,n) -> ( 
    result := false;
    for j in s do (if j#1 == n then result = j#0);
    result
)

findCol = method()
findCol (List,ZZ) := (s,n) -> ( 
    result := false;
    for j in s do (if j#0 == n then result = j#1);
    result
)

computeDualForm = method()
computeDualForm RankConditions := (S) -> (
     rep := getRep2 S;
     k := #S.dims - 1;
     forms := {};
     for i in 0..k-1 do (
          forms = append(forms, apply(apply(toList(0..(S.dims#i-1)),m -> (apply(toList(0..(S.dims#(i+1)-1)),p -> (m,p)))),n -> apply(n, q -> if member(q,rep#i) then 1 else 0)))    
     );
     forms = apply(forms,x-> new MatrixForm from x);
     
     duals = apply(toList(0..k-1), i -> apply(toList(0..(S.dims#(i+1)-1)),m -> (apply(toList(0..(S.dims#i-1)),p -> 1))));
     duals = apply(duals, x -> new MatrixForm from x);

     connections := new MutableHashTable from {};

     --AA' = 0
     for x in rep#0 do duals#0#(x#1) = new MutableList from toList(apply(0..S.dims#0-1, x -> 0));

     --NN' = M'M
     for i in 0..k-2 do (
        mm := duals#i * forms#i;
        nn := forms#(i+1) * duals#(i+1); 
        size := S.dims#(i+1)-1;

        for a in 0..size do (
            for b in 0..size do (
                if mm#a#b == 1 and nn#a#b == 1 then connections#(i+1,findCol(rep#(i+1),a),b) = (i,a,findRow(rep#i,b)) 
                else if mm#a#b == 0 and nn#a#b == 1 then duals#(i+1)#(findCol(rep#(i+1),a))#b = 0
                else if mm#a#b == 1 and nn#a#b == 0 then (
                    duals#i#a#(findRow(rep#i,b)) = 0;
                    x := (i,a,findRow(rep#i,b));
                    while connections#?x do (
                         y := connections#x;
                         duals#(y#0)#(y#1)#(y#2) = 0;
                         x = y
                    );
                );
            );
        );
     );

     --0 = Z'Z
     for x in rep#-1 do (
         for y in 0..S.dims#(-1)-1 do (
             duals#-1#y#(x#0) = 0;
             z := (#duals - 1,y,x#0);
             while connections#?z do (
                zz := connections#z;
                duals#(zz#0)#(zz#1)#(zz#2) = 0;
                z = zz;
             );
         );
     );

     duals    
)

computeDualOrbitFromForm = method()
computeDualOrbitFromForm List := (duals) -> (
     G := new GraphList from {new MutableHashTable from {},new MutableHashTable from {}};
     a := #(duals#0#0);
     dimensions := {a};
     edges := {};
     for i in 1..a do edges = join(edges, {{i,i+a},{-1,i},{i+a,-2}});
     i := 2*a+1;
     j := -3;
     
     for dual in duals do (
         a = #dual;
         dimensions = append(dimensions,a);
         b = #(dual#0);
         edges = join (edges, toList ({j,i}..{j,i+a-1}));
         edges = join (edges, toList ({i+a,j-1}..{i+2*a-1,j-1}));
         edges = join (edges, apply(toList (0..a-1), x -> {i+x,i+a+x}));
         for row in 0..a-1 do (
             for col in 0..b-1 do (
                 if dual#row#col == 1 then edges = append(edges, {i-b+col,i+row});
             );
         );
         i = i + 2*a;
         j = j-2;
     );
     
     vertices = join(toList(j+1..-1),toList(1..i-1));

     for edge in edges do (
        if G#children#?(edge#0) then G#children#(edge#0) = append(G#children#(edge#0),edge#1) else G#children#(edge#0) = {edge#1};
        if G#parents#?(edge#1) then G#parents#(edge#1) = append(G#parents#(edge#1),edge#0) else G#parents#(edge#1) = {edge#0};
     );

     levels := apply(toList (1..#duals), i -> apply(toList (0..#duals-i), j -> fordFulkersonUnitCapacity(G,-1-2*j, -2-2*(i+j))));

     new RankConditions from (dimensions, levels)
)


computeDual = method()
computeDual RankConditions := computeDualOrbitFromForm @@ computeDualForm


findPathologies = method()
findPathologies List := (dimensions) -> (
  listOfOrbits := getAllStrataR dimensions;
  pathologies := {};
  for c1 in listOfOrbits do (
        for c2 in listOfOrbits do (
            if c1 < c2 and (computeDual c1) < (computeDual c2) then pathologies = append(pathologies, (c1,c2))
        );
  );
  pathologies
)



