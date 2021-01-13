Unprotect["WeylMultiwaySystem`*"];

ClearAll @@ (# <> "*" & /@ Contexts["WeylMultiwaySystem`*"]);

BeginPackage["WeylMultiwaySystem`"];

(* Unsure what else I need in here *) 

EndPackage[];

SetAttributes[#, {Protected, ReadProtected}] & /@ Evaluate @ Names @ "WeylMultiwaySystem`*";