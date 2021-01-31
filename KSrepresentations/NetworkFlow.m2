----------------------------------------------------------------------------------------------------------------
-- 
-- File: NetworkFlow.m2
-- Author: Reginald Lybbert
--
-- This file provides an implementation of the Ford-Fulkerson algorithm
-- for finding maximum flow through a network in the special case that
-- each edge in the network has unit capacity.  This is used in the
-- computation of dual orbits.
--
-- This implementation is not intended for general use, only as a part
-- of the compute dual function, although other uses may be found.
--
----------------------------------------------------------------------------------------------------------------


GraphList = new Type of MutableHashTable
new GraphList from List := {GraphList, ls} -> hashTable {children => ls#0, parents => ls#1}

edges = method();
edges GraphList := (G) -> (
    result = {};
    for x in keys (G#children) do (
        result = join(result, apply(G#children#x, i -> {x,i}))
    );
    result
);

children = method()
children (GraphList, Thing) := (G,x) -> if (G#children)#?x then (G#children)#x else {};

parents = method()
parents (GraphList, Thing) := (G,x) -> if (G#parents)#?x then (G#parents)#x else {};

fordFulkersonUnitCapacity = method()
fordFulkersonUnitCapacity (GraphList, ZZ, ZZ) := (G,s,t) -> (
    f := 0;
    again := true;
    flow := new MutableHashTable from {};
    for edge in (edges G) do (
        flow#edge = 0;
    );
    while again do (
         again = false;
         labels := new MutableHashTable from {};
         ls := {s};
         while #ls > 0 do (
              x := ls#0;
              ls = drop(ls,1);
              for y in children (G,x) do (
                  if not labels#?y and flow#{x,y} == 0 then (
                     labels#y = (x,0);
                     ls = append(ls, y); 
                  );
              );
              for y in parents (G,x) do (
                  if not labels#?y and flow#{y,x} == 1 then (
                      labels#y = (x,1);
                      ls = append(ls, y);
                  );
              );
              if labels#?t then (
                   again = true;
                   f = f+1;
                   ls = {};
                   x = t;
                   while x != s do (
                       if labels#x#1 == 0 then (
                            flow#({labels#x#0,x}) = 1;
                            x = labels#x#0;
                       )
                       else if labels#x#1 == 1 then (
                            flow#({x,labels#x#0}) = 0;
                            x = labels#x#0;
                       );
                   );
              );
         );
    );
    f
)

