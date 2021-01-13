BeginPackage["WeylMultiwaySystem`"]

(* PackageImport["importpackages"] *) 

PackageExport["WeylMultiwayFinite"]

SetUsage @ "
WeylMultiwayFinite[FiniteAlgebra, r, n, Options] gives the Multiway System of a
Finite Algebra {An,Bn,Cn,Dn,E6,E7,E8,F4,G2} of rank, r, after n iterations, with Options
matching those of MultiwayFunctionSystem.   
";

A::"usage"="represents the simple root systems of series A, for r>0";
B::"usage"="represents the simple root systems of series B, for r>2";
C::"usage"="represents the simple root systems of series C, for r>2";
D::"usage"="represents the simple root systems of series D, for r>3";
E::"usage"="represents the simple root systems of series E, for r=6,7,8";
F::"usage"="represents the simple root systems of series F, for r=4";
G::"usage"="represents the simple root systems of series G, for r=2";
 

(*)
SyntaxInformation[WeylMultiwayFinite] = {"ArgumentsPattern" -> {_}};
*)

 (* Inner product *) 
InnerProduct[x_, y_] := Total[x y];
 
(* Defining functions to give all the semi-simple algebra's simple \
roots *) 
SimpleRoots["A", r_] := 
  Table[UnitVector[r + 1 , i] - UnitVector[r + 1 , i + 1], {i, 1, r}]; 

SimpleRoots["B", r_] := 
  Join[Table[
    UnitVector[r , i] - UnitVector[r , i + 1], {i, 1, r - 1}], {UnitVector[r, r]}];

SimpleRoots["C", r_] := 
  Join[Table[
    UnitVector[r , i] - UnitVector[r , i + 1], {i, 1, r - 1}], {2 UnitVector[r, r]}]; 

SimpleRoots["D", r_] := 
  Join[Table[
    UnitVector[r , i] - UnitVector[r , i + 1], {i, 1, r - 1}], {UnitVector[r, r - 1] + UnitVector[r, r]}];

SimpleRoots["E", 6] := {{1, -1, -1, -1, -1, -1, -1, 1},
    {2, 2, 0, 0, 0, 0, 0, 0},
    {-2, 2, 0, 0, 0, 0, 0, 0},
    {0, -2, 2, 0, 0, 0, 0, 0},
    {0, 0, -2, 2, 0, 0, 0, 0},
    {0, 0, 0, -2, 2, 0, 0, 0}};

SimpleRoots["E", 7] := {{1, -1, -1, -1, -1, -1, -1, 1},
    {2, 2, 0, 0, 0, 0, 0, 0},
    {-2, 2, 0, 0, 0, 0, 0, 0},
    {0, -2, 2, 0, 0, 0, 0, 0},
    {0, 0, -2, 2, 0, 0, 0, 0},
    {0, 0, 0, -2, 2, 0, 0, 0},
    {0, 0, 0, 0, -2, 2, 0, 0}};

SimpleRoots["E", 8] := {{1, -1, -1, -1, -1, -1, -1, 1}, 
    {2, 2, 0, 0, 0, 0, 0, 0},
    {-2, 2, 0, 0, 0, 0, 0, 0},
    {0, -2, 2, 0, 0, 0, 0, 0},
    {0, 0, -2, 2, 0, 0, 0, 0},
    {0, 0, 0, -2, 2, 0, 0, 0},
    {0, 0, 0, 0, -2, 2, 0, 0},
    {0, 0, 0, 0, 0, -2, 2, 0}};

SimpleRoots["F", 4] := {{0, 2, -2, 0}, {0, 0, 2, -2},{0, 0, 0, 2}, {1, -1, -1, -1}};

SimpleRoots["G", 2] := {{1, -1, 0}, {-2, 1, 1}};


(* Function defining the action of the Weyl group *)
W[ai_, j_] := j - 2 (InnerProduct[ai, j]/InnerProduct[ai, ai]) ai ;


(* Function that gives a Multiway System for the sucessive action of the Weyl
group as defined above  *)
 WeylMultiwayFinite[Alg_, r_, n_, options___] := 
 With[{simpleRoots = root @@@ SimpleRoots[Alg, r]},
  ResourceFunction["MultiwayFunctionSystem"][
   Function[reflectedRoot, 
    root @@ W[List @@ #, List @@ reflectedRoot] & /@ simpleRoots], 
   simpleRoots, n, options]]

EndPackage[];