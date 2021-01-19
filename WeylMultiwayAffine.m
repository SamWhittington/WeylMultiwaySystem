Package["WeylMultiwaySystem`"]

PackageImport["GeneralUtilities`"]

PackageExport["WeylMultiwayAffine"]


SetUsage @ "
WeylMultiwayAffine[AffineAlgebra, r, m, n, Options] gives the Multiway System of a 
Affine Algebra {An,Bn,Cn,Dn,E6,E7,E8,F4,G2} of rank, r, after n iterations, 
and reflection in m-plane.Options are matching those of MultiwayFunctionSystem.   
";


(*)
SyntaxInformation[WeylMultiwayFinite] = {"ArgumentsPattern" -> {_}};
*)

 (* Inner product *) 
InnerProduct[x_, y_] := Total[x y];
 
(* Defining functions to give all the semi-simple algebra's simple \
roots as private functions *) 

RootsA[r_] := 
  Table[UnitVector[r + 1 , i] - UnitVector[r + 1 , i + 1], {i, 1, 
    r}]; (* Simple roots *) 
SimpleRootsAffine["A", r_] := 
 Join[RootsA[r], {-Sum[RootsA[r][[i]], {i, 1, r}]}]


RootsB[r_] := 
  Join[Table[
    UnitVector[r , i] - UnitVector[r , i + 1], {i, 1, 
     r - 1}], {UnitVector[r, r]}]; (* Simple roots *) 
SimpleRootsAffine["B", r_] := 
  Join[RootsB[r], -{UnitVector[r, 1] + UnitVector[r, 2]}];

RootsC[r_] := 
  Join[Table[
    UnitVector[r , i] - UnitVector[r , i + 1], {i, 1, 
     r - 1}], {2 UnitVector[r, r]}]; (* Simple roots *) 
SimpleRootsAffine["C", r_] := Join[RootsC[r], -2 {UnitVector[r, 1]}];

SimpleRootsAffine["D", r_] := 
  Join[Table[
    UnitVector[r , i] - UnitVector[r , i + 1], {i, 1, 
     r - 1}], {UnitVector[r, r - 1] + UnitVector[r, r]}];

SimpleRootsAffine["E", 
   6] := {2 {1/
     2, -(1/2), -(1/2), -(1/2), \[Minus](1/2), \[Minus](1/
      2), \[Minus](1/2), 1/2},
   2 {1, 1, 0, 0, 0, 0, 0, 0},
   2 {-1, 1, 0, 0, 0, 0, 0, 0},
   2 {0, -1, 1, 0, 0, 0, 0, 0},
   2 {0, 0, -1, 1, 0, 0, 0, 0},
   2 {0, 0, 0, -1, 1, 0, 0, 0},
   2 {-(1/2), -(1/2), -(1/2), -(1/2), -(1/2), 1/2, 1/2, -(1/2)}};

SimpleRootsAffine["E", 
   7] :=  {2 {1/
     2, -(1/2), -(1/2), -(1/2), \[Minus](1/2), \[Minus](1/
      2), \[Minus](1/2), 1/2},
   2 {1, 1, 0, 0, 0, 0, 0, 0},
   2 {-1, 1, 0, 0, 0, 0, 0, 0},
   2 {0, -1, 1, 0, 0, 0, 0, 0},
   2 {0, 0, -1, 1, 0, 0, 0, 0},
   2 {0, 0, 0, -1, 1, 0, 0, 0},
   2 {0, 0, 0, 0, -1, 1, 0, 0},
   2 {0, 0, 0, 0, 0, 0, 1, -1}};

SimpleRootsAffine["E", 
   8] := {{1, -1, -1, -1, \[Minus]1, \[Minus]1, \[Minus]1, 1},
   {2, 2, 0, 0, 0, 0, 0, 0},
   {-2, 2, 0, 0, 0, 0, 0, 0},
   {0, -2, 2, 0, 0, 0, 0, 0},
   {0, 0, -2, 2, 0, 0, 0, 0},
   {0, 0, 0, -2, 2, 0, 0, 0},
   {0, 0, 0, 0, -2, 2, 0, 0},
   {0, 0, 0, 0, 0, -2, 2, 0},
   {0, 0, 0, 0, 0, 0, -2, -2}};

SimpleRootsAffine["F", 4] := {2 {0, 1, -1, 0}, 2 {0, 0, 1, -1}, 
   2 {0, 0, 0, 1}, 2 {1/2, -(1/2), -(1/2), -(1/2)}, 2 {-1, -1, 0, 0}};

SimpleRootsAffine["G", 2] := {{1, -1, 0}, {-2, 1, 1}, {1, 1, -2}};


 (* Reflections in a hyperplane orthogonal to the roots - for a rank, \
r algebra *)
WAffine[ai_, j_, m_] := 
  j - 2 (InnerProduct[ai, j]/InnerProduct[ai, ai]) ai  + 
   2 (m/InnerProduct[ai , ai ]) ai;

 WeylMultiwayAffine[Alg_, r_, m_, n_, options___] := 
  With[{simpleRoots = root @@@ SimpleRootsAffine[Alg, r]},
   ResourceFunction["MultiwayFunctionSystem"][
    Function[reflectedRoot, 
     root @@ WAffine[List @@ #, List @@ reflectedRoot, m] & /@ 
      simpleRoots], simpleRoots, n, options]];

AffineWeylReflectionsGraph[Alg_, r_, m_, n_, options___] := 
  With[{simpleRoots = root @@@ SimpleRootsAffine[Alg, r]},
    NestGraph[
    Function[reflectedRoot, 
     root @@ WAffine[List @@ #, List @@ reflectedRoot, m] & /@ 
      simpleRoots], simpleRoots, n, options]];