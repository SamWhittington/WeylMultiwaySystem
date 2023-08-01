(* ::Package:: *)

Package["WeylMultiwaySystem`"]

PackageImport["GeneralUtilities`"]

PackageExport["WeylMultiwayHyperbolic"]


SetUsage @ "
WeylMultiwayHyperbolic[HyperbolicAlgebra, r, n, m, x, Options] gives the Multiway System of a 
Hyperbolic Kac-Moody Algebra {An,Bn,Cn,Dn,E6,E7,E8,F4,G2} of finite rank, r, with an n-extension, after x iterations, 
and reflection in m-plane. Options are matching those of MultiwayFunctionSystem.   
";


(*
SyntaxInformation[WeylMultiwayFinite] = {"ArgumentsPattern" -> {_}};
*)

 
 (* Inner product finite *) 
InnerProductFinite[x_, y_] := Total[x y];
 
   (* Inner product 
P[x_,y_]:= \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(r + 1\)]\(x[\([i]\)]y[\([i]\)]\)\)-\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]\((x[\([r + 1 + 2 i - 1]\)]y[\([r + 1 + 2 i\ ]\)] + y[\([r + 1 + 2 i\  - \ 1]\)]x[\([r + 1 + 2 i]\)])\)\);*) 
InnerProduct[x_,y_,r_,n_]:= Sum[x[[i]]y[[i]],{i,1,r+1}]- Sum[x[[r+1+2i-1]]y[[r+1+2i]]+y[[r+1+2i-1]]x[[r+1+2i]],{i,n}]
(* Defining functions to give all the semi-simple algebra's simple roots as private functions *) 

RootsAFinite[r_,n_] := 
  Table[UnitVector[r+1 +2n,i] - UnitVector[r+1 +2n,i+1], {i,1,r}]; (* Simple roots *) 

ExtendedPartK[r_,n_] := Table[Join[ ConstantArray[0,(r+1)], UnitVector[2n ,(2*i+1)]], {i,0,n-1}];

ExtendedPartKBar[r_,n_] := Table[Join[ ConstantArray[0,(r+1)],-UnitVector[2n ,(2*j)]], {j,1,n}];

RootsAAffineRoot[r_,n_]:={-Sum[RootsAFinite[r,n][[i]], {i, 1, r}]+ExtendedPartK[r,n][[1]]};

RootsAExtendedParts[r_,n_] := Join[{-(ExtendedPartK[r,n][[1]] + ExtendedPartKBar[r,n][[1]])},Table[ExtendedPartK[r,n][[j-1]] - (ExtendedPartK[r,n][[j]] + ExtendedPartKBar[r,n][[j]]) ,{j,2,n}]];

(* RootsAExtendedPartHighestRoot[] := *)

(* Join all the roots with finite and extended parts *)

SimpleRootsHyperbolic["A", r_,n_] := 
 Join[RootsAFinite[r,n], RootsAAffineRoot[r,n],RootsAExtendedParts[r,n]];

CartanA["A",r_,n_]:= Table[InnerProduct[SimpleRootsHyperbolic["A", r,n][[i]],SimpleRootsHyperbolic["A", r,n][[j]]  ,r,n],{i,1,r+n+1},{j,1,r+n+1}]//MatrixForm;


RootsB[r_] := 
  Join[Table[
    UnitVector[r , i] - UnitVector[r , i + 1], {i, 1, 
     r - 1}], {UnitVector[r, r]}];(* Simple roots *) 
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


(*(* Reflections in a hyperplane orthogonal to the roots - for a rank, r algebra *)

HyperbolicW[ai_,j_,m_] := j - 2(InnerProduct[ai,j]/InnerProduct[ai,ai])ai  +2(m/InnerProduct[ai ,ai ])ai;

WeylMultiwayHyperbolic[Alg_,r_,m_,n_,options___]:=With[{simpleRoots = root@@@SimpleRootsHyperbolic[Alg,r,n]},ResourceFunction["MultiwayFunctionSystem"][Function[reflectedRoot,root@@HyperbolicW[List@@#,List@@reflectedRoot,m]&/@simpleRoots],simpleRoots,n, options]];

HyperbolicWeylReflectionsGraph[Alg_,r_,m_,n_,options___]:=With[{simpleRoots = root@@@ SimpleRootsHyperbolic[Alg,r]},NestGraph[Function[reflectedRoot,root@@HyperbolicW[List@@#,List@@reflectedRoot,m]&/@simpleRoots],simpleRoots,n, options]]; *)






(* Reflections in a hyperplane orthogonal to the roots - for a rank, r algebra *)

HyperbolicW[ai_,j_,r_,n_,m_] := j - 2(InnerProduct[ai,j,r,n]/InnerProduct[ai,ai,r,n])ai+2(m/InnerProduct[ai,ai,r,n])ai;

WeylMultiwayHyperbolic[Alg_,r_,n_,m_,iterations_,options___]:=With[{simpleRoots = root@@@SimpleRootsHyperbolic[Alg,r,n]},ResourceFunction["MultiwayFunctionSystem"][Function[reflectedRoot,root@@HyperbolicW[List@@#,List@@reflectedRoot,r,n,m]&/@simpleRoots],simpleRoots,iterations, options]];

HyperbolicWeylReflectionsGraph[Alg_,r_,n_,m_,iterations_,options___]:=With[{simpleRoots = root@@@ SimpleRootsHyperbolic[Alg,r,n]},NestGraph[Function[reflectedRoot,root@@HyperbolicW[List@@#,List@@reflectedRoot,r,n,m]&/@simpleRoots],simpleRoots,iterations, options]];
