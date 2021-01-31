-------------------------------------------------------------------------------------------------------------------------
-- File: PSNF.m2
-- Author: Nicole Kitt
--
-- This file contains the function partialSmithNormalForm(Matrix).
--
-- The function partialSmithNormalForm(Matrix) inputs a matrix M over a polynomial ring, and outputs the matrix obtained
-- by performing row/column operations, excluding division, on M. More precisely, this function reduces M to an upper
-- left identity block, and lower right block whose entries consist of non-constant polynomials. 
--
-- The abbreviation PSNF stands for "partial Smith normal form". We refer to the process of using elementary row and
-- column operations to reduce a matrix to an upper left identity block and lower right block with no nonzero scalar 
-- entries, as computing partial Smith normal form. This is the definition of "partial Smith form", which can be found
-- in https://math.usask.ca/~bremner/research/publications/Slides-for-my-LinearAlgebra-Talk%20copy.pdf 
--
-- Note: This function will not work over non-polynomial rings, such as QQ. We can probably adapt it so that it works 
-- over QQ, but this extra functionality is not necessary for the project. 
-------------------------------------------------------------------------------------------------------------------------

partialSmithNormalForm = method()
partialSmithNormalForm(Matrix) := (M) -> (  
    numrow := numrows M;
    numcol := numcols M;
    R := ring M;
    R = frac R;
    M = sub(M,R);
    M = mutableMatrix M;
    i=0;
--- big outside loop going through all columns
    for j in 0..(numcol - 1) do (
-- finding most left column that contains next pivot entry
        L :={};
        br = j;
        for c in j..(numcol-1) do (
          for r in 0..(numrow -1) do ( if (M_(r,c)!=0 and #(set support M_(r,c))==0) then L=append(L,r) );
          if (#L!=0 or br == (numcol-1)) then break L;
          br = br + 1;
        );
-- moving that entry to make pivot
        for r in 0..(numrow -1) do ( 
          if (M_(r,br)!=0 and #(set support M_(r,br))==0) then (rowMult(M,r,1/M_(r,br)) and rowSwap(M,i,r) and columnSwap(M,j,br) and break);
        );
-- using pivot to remove all entries in that colummn
        for k in (i+1)..(numrow -1) do ( if (M_(i,j)!=0 and #(set support M_(i,j))==0) then rowAdd(M,k,-M_(k,j),j) );
-- using pivot to remove all entries in that row
        for k in (i+1)..(numcol -1) do ( if (M_(i,j)!=0 and #(set support M_(i,j))==0) then columnAdd(M,k,-M_(i,k),i) );
    if i == (numrow-1) then break M;
    i=i+1;
    );
    M = matrix M;
    M = compress M;
    M = transpose compress transpose M;
    M = mutableMatrix M;
    M
)


