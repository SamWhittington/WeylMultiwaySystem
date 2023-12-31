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
 
   (* Inner product - Different for A and others are the same due to their dimension of representation due to Bourbaki
P[x_,y_]:= \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(r + 1\)]\(x[\([i]\)]y[\([i]\)]\)\)-\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]\((x[\([r + 1 + 2 i - 1]\)]y[\([r + 1 + 2 i\ ]\)] + y[\([r + 1 + 2 i\  - \ 1]\)]x[\([r + 1 + 2 i]\)])\)\);*) 
InnerProduct["A",x_,y_,r_,n_]:= Sum[x[[i]]y[[i]],{i,1,r+1}]- Sum[x[[r+1+2i-1]]y[[r+1+2i]]+y[[r+1+2i-1]]x[[r+1+2i]],{i,n}];
InnerProduct["B",x_,y_,r_,n_]:= Sum[x[[i]]y[[i]],{i,1,r}]- Sum[x[[r+2i-1]]y[[r+2i]]+y[[r+2i-1]]x[[r+2i]],{i,n}];
InnerProduct["C",x_,y_,r_,n_]:= Sum[x[[i]]y[[i]],{i,1,r}]- Sum[x[[r+2i-1]]y[[r+2i]]+y[[r+2i-1]]x[[r+2i]],{i,n}];
InnerProduct["D",x_,y_,r_,n_]:= Sum[x[[i]]y[[i]],{i,1,r}]- Sum[x[[r+2i-1]]y[[r+2i]]+y[[r+2i-1]]x[[r+2i]],{i,n}];
InnerProduct["D",x_,y_,r_,n_]:= Sum[x[[i]]y[[i]],{i,1,r}]- Sum[x[[r+2i-1]]y[[r+2i]]+y[[r+2i-1]]x[[r+2i]],{i,n}];
InnerProduct["E",x_,y_,r_,n_]:= Sum[x[[i]]y[[i]],{i,1,r}]- Sum[x[[r+2i-1]]y[[r+2i]]+y[[r+2i-1]]x[[r+2i]],{i,n}];
InnerProduct["F",x_,y_,r_,n_]:= Sum[x[[i]]y[[i]],{i,1,r}]- Sum[x[[r+2i-1]]y[[r+2i]]+y[[r+2i-1]]x[[r+2i]],{i,n}];
InnerProduct["G",x_,y_,r_,n_]:= Sum[x[[i]]y[[i]],{i,1,r}]- Sum[x[[r+2i-1]]y[[r+2i]]+y[[r+2i-1]]x[[r+2i]],{i,n}];

(* For the n-extended E series which has a finite basis of dimension 8*)
InnerProduct["En",x_,y_,r_,n_]:= Sum[x[[i]]y[[i]],{i,1,8}]- Sum[x[[(r+2i-1)+(8-r)]]y[[r+2i+(8-r)]]+y[[(r+2i-1)+(8-r)]]x[[r+2i+(8-r)]],{i,n}];


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

Cartan[alg_,r_,n_]:= Table[InnerProduct[alg,SimpleRootsHyperbolic[alg, r,n][[i]],SimpleRootsHyperbolic[alg, r,n][[j]]  ,r,n],{i,1,r+n+1},{j,1,r+n+1}]//MatrixForm;


(* n-extended simple root systems *)

(* A Series is given in detail above *)

(* General n-extended parts of the root structure - in r+2n diemensions - for all except the representation of the A series *)
RootsExtendedParts[r_,n_] := Join[{-(ExtendedPartK[r-1,n][[1]] + ExtendedPartKBar[r-1,n][[1]])},Table[ExtendedPartK[r-1,n][[j-1]] - (ExtendedPartK[r-1,n][[j]] + ExtendedPartKBar[r-1,n][[j]]) ,{j,2,n}]];
 
(* B *)
RootsBFinite[r_,n_] := 
  Join[Table[
    UnitVector[r +2n , i] - UnitVector[r + 2n , i + 1], {i, 1, 
     r - 1}], {UnitVector[r+2n, r]}];(* Simple roots *)
     
RootsBAffineRoot[r_,n_]:={-Sum[RootsBFinite[r,n][[i]], {i, 1, r}]+ExtendedPartK[r-1,n][[1]]};

SimpleRootsHyperbolic["B", r_,n_] := 
 Join[RootsBFinite[r,n], RootsBAffineRoot[r,n],RootsExtendedParts[r,n]];
 
 
(* C *)
RootsCFinite[r_,n_] := 
  Join[Table[
    UnitVector[r+2n , i] - UnitVector[r+2n , i + 1], {i, 1, 
     r - 1}], {2 UnitVector[r+2n, r]}]; (* Simple roots *) 
     
RootsCAffineRoot[r_,n_]:={-Sum[RootsCFinite[r,n][[i]], {i, 1, r}]+ExtendedPartK[r-1,n][[1]]}
SimpleRootsHyperbolic["C", r_,n_] := 
 Join[RootsCFinite[r,n], RootsCAffineRoot[r,n],RootsExtendedParts[r,n]];
 

(* D *)
SimpleRootsAffine["D", r_,n_] := 
  Join[Table[
    UnitVector[r+2n, i] - UnitVector[r+2n , i + 1], {i, 1, 
     r+2n - 1}], {UnitVector[r+2n, r - 1] + UnitVector[r+2n, r]+ExtendedPartK[r-1,n][[1]]}];

(* Extend directly from the affine set, more direct than the way done for A,B and C above*)
SimpleRootsHyperbolic["D", r_,n_] := 
 Join[SimpleRootsAffine["D",r,n],RootsExtendedParts[r,n]];
 
 
 
 (* TODO: NEXT: Work out how to pad/ get the E-series arrays in the correct positions for the n-extensions -->> Probably best to look in the code
 / notes for my old work because will have certainly done this somewhere else before... -->> Look where I have generated the 
 n-exteneded roots for the E-series and do this here below! *)  
       
(* E Series' *)
(* N.B. For the E-series its easiest to build from the finite basis and add the highest root as in our paper (not as in Bourbaki) *)

SimpleRoots["E", 6] := {{1, -1, -1, -1, -1, -1, -1, 1},
    {2, 2, 0, 0, 0, 0, 0, 0},
    {-2, 2, 0, 0, 0, 0, 0, 0},
    {0, -2, 2, 0, 0, 0, 0, 0},
    {0, 0, -2, 2, 0, 0, 0, 0},
    {0, 0, 0, -2, 2, 0, 0, 0}};
    
    
SimpleRootsAffine["E", 
   6,n_] := ArrayPad[{{0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5,0,0},
   {1, 1, 0, 0, 0, 0, 0, 0,0,0},
    {-1, 1, 0, 0, 0, 0, 0, 0,0,0},
    {0, -1, 1, 0, 0, 0, 0, 0,0,0},
    {0, 0, -1, 1, 0, 0, 0, 0,0,0},
    {0, 0, 0, -1, 1, 0, 0, 0,0,0},
    {0, 0, 0, 0, -1, 1, 0, 0,1,0}},{{0,0},{0,(2n-2)}}];
   
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

(* Should work this way *)
  SimpleRootsAffine["E", 8] :={{1/2,-(1/2),-(1/2),-(1/2),-(1/2),-(1/2),-(1/2),1/2,0,0},     (* Base rep for E9*)
{1,1,0,0,0,0,0,0,0,0},
{-1,1,0,0,0,0,0,0,0,0},
{0,-1,1,0,0,0,0,0,0,0},
{0,0,-1,1,0,0,0,0,0,0},
{0,0,0,-1,1,0,0,0,0,0},
{0,0,0,0,-1,1,0,0,0,0}, 
{0,0,0,0,0,-1,1,0,0,0},
{0,0,0,0,0,0,-1,-1,0,1}}; 
(* (* Way that I don't think works... *)
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
   *)
   
(* General hyperbolic roots for the E series' *)
SimpleRootsHyperbolic["E", r_,n_] := 
 Join[SimpleRootsAffine["E",r,n],RootsExtendedParts[r,n]];
 
 (* Trying the above extension again now we have the correct affine system - 8-dimensional basis *)
SimpleRootsHyperbolic["En", r_,n_] := 
 Join[SimpleRootsAffine["E",r,n],RootsExtendedParts[8,n]];
 
(* F *) 
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
HyperbolicW[alg_,ai_,j_,r_,n_,m_] := j - 2(InnerProduct[alg,ai,j,r,n]/InnerProduct[alg,ai,ai,r,n])ai+2(m/InnerProduct[alg,ai,ai,r,n])ai;

WeylMultiwayHyperbolic[Alg_,r_,n_,m_,iterations_,options___]:=With[{simpleRoots = root@@@SimpleRootsHyperbolic[Alg,r,n]},ResourceFunction["MultiwayFunctionSystem"][Function[reflectedRoot,root@@HyperbolicW[Alg,List@@#,List@@reflectedRoot,r,n,m]&/@simpleRoots],simpleRoots,iterations, options]];

HyperbolicWeylReflectionsGraph[Alg_,r_,n_,m_,iterations_,options___]:=With[{simpleRoots = root@@@ SimpleRootsHyperbolic[Alg,r,n]},NestGraph[Function[reflectedRoot,root@@HyperbolicW[Alg,List@@#,List@@reflectedRoot,r,n,m]&/@simpleRoots],simpleRoots,iterations, options]];
