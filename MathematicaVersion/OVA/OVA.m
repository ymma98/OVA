(* ::Package:: *)

                       
(****************************************************************)
(*                                                              *)
(*       Vector Analysis in Orthogonal Coordinates System      *)
(*                                                              *)
(*       Version 1.0                                            *)
(*                                                              *)
(****************************************************************)



(*
This is a modified package based on Prof. Qin's GVA(General Vector Analysis)
package Version1.0 in 1997. The detailed copyright of Qin's GVA Version 1.0
are shown below in its original state.

Usually we want to normalize the basis vectors in orthogal curvilinear 
coordinates (for example the vector analysis formulas in cylindrical and
spherical coordinates as what we usually see in many textbook appendices)
and keep the basis vectors unnormalized in non-orthogonal coordinates 
(It seems meaningless to normalize the basis vectors in non-orthogonal 
coordinates because the co-variant and contra-variant vectors are still 
different with each other after the normalization, while the normalization
 process will make the formula derivation much more tedious). 
 
In the OVA, we only perform the symbolic vector analysis in the orthogonal
coordinate system. All functions in OVA` context start with 
lowercase letters to avoid conflicts with symbols in System`.

ymma98,
ymma98@qq.com
*)


(***************************************************************************)    
(*                                                                         *)
(*       Copyright: Hong Qin, January 1, 1997                              *) 
(*                                                                         *)
(*              All rights reserved                                        *)
(*                                                                         *)
(*                                                                         *)
(* This software package and its accompanying documentation are provided   *)
(* as is, without guarantee of support or maintenance.  The copyright      *)
(* holder makes no express or implied warranty of any kind with respect    *)
(* to this software, including implied warranties of merchantability or    *)
(* fitness for a particular purpose, and is not liable for any damages     *)
(* resulting in any way from its use.                                      *)
(*                                                                         *)
(* Everyone is granted permission to copy, modify and redistribute this    *)
(* software package and its accompanying documentation, provided that:     *)
(*  1. All copies contain this notice in the main program file and in the  *)
(*     supporting documentation.                                           *)
(*  2. All modified copies carry a prominent notice stating who made the   *)
(*     last modification and the date of such modification.                *)
(*  3. No charge is made for this software or works derived from it, with  *)
(*     the exception of a distribution fee to cover the cost of materials  *)
(*     and/or transmission.                                                *)
(*                                                                         *)
(***************************************************************************)




(****************************************************************)
(*                                                              *)
(* Basic Functions:                                             *)
(*                                                              *)
(*  (1) Vector analysis in Orthogonal coordinates               *)
(*  (2) Define your own coordinates to work on                  *)
(*  (3) Perturbative vector analysis with respect to any small  *)
(*             parameter, e.g. inverse aspect ratio in tokamak, *)
(*             to any order                                     *)
(*  (4) Vector analysis in Shafranov coordinate,                *)
(*             flux coordinate, straight tokamak coordinate,    *)
(*             circular concentric tokamak coordinate, ......   *)  
(*                                                              *)   
(****************************************************************)



(************ Begin OVA` Package *****)
BeginPackage["OVA`"]

OVA`setCoordinateSystem::usage=
"setCoordinateSystem[$CoordinateSystem,p1,p2][c1,c2,c3] will set up a 
$coordinateSystem coordinate with coordinate c1,c2,c3 and parameter p1,p2."

OVA`defineCoordinateSystem::usage=
"defineCoordinateSystem[$CoordinateSystem,p1,p2,...][c1,c2,c3][MetricMatrix] will
add the $coordinateSystem with parameter p1,p2,... to the package. The metric 
tensor as a function of coordinate c1,c2,c3 is given by metricMatrix."


OVA`metric::usage=
"metric[1][c1,c2,c3] is the metric matrix g_{i,j},
 metric[2][c1,c2,c3] is the metric matrix g^{i,j}."

OVA`jacobian::usage="jacobian[c1,c2,c3] gives the jacobian.  "

OVA`christoffel::usage=
"GChristoffe1[1][c1,c2,c3] is the 1st kind Christoffel symbol,
 GChristoffe1[2][c1,c2,c3] is the 2nd kind Christoffel symbol."

OVA`defineVector::usage=
"defineVector[A1[c1,c2,c3],A2[c1,c2,c3],A3[c1,c2,c3]]
 creats a vector with A1,A2,A3 as its components. "

OVA`dot::usage=
"dot[A[c1,c2,c3],B[c1,c2,c3]] gives the dot product. "

OVA`cross::usage=
"cross[A[c1,c2,c3],B[c1,c2,c3]] gives the cross product."

OVA`grad::usage=
"grad[f[c1,c2,c3]] gives the gradiant of f."

OVA`div::usage=
"div[A[c1,c2,c3]] gives the divergence of A."

OVA`curl::usage=
"curl[A[c1,c2,c3]] gives the curl of A."

OVA`covariantDerivative::usage=
"covariantDerivative[A[c1,c2,c3],x[c1,c2,c3]] gives the directional
 derivative of x in A direction, 
 x can be either a vector or a scalar."

OVA`laplacian::usage=
"laplacian[f[c1,c2,c3]] gives the laplacian of scalar f."

OVA`absoluteValue::usage=
"absoluteValue[A[c1,c2,c3]] gives the absolute
 value of vector A[c1,c2,c3]. "

OVA`unitVector::usage=
"unitVector[A[c1,c2,c3]] gives the normalized 
 vector in the direction of A[c1,c2,c3]. "

OVA`parallel::usage=
"parallel[A[c1,c2,c3],B[c1,c2,c3]] gives the 
 projection of A[c1,c2,c3] in the direction of B[c1,c2,c3].
 It returns a scalar."

OVA`parallelUnitVector::usage=
"parallel[A[c1,c2,c3],B[c1,c2,c3]] gives the 
 the projection of A[c1,c2,c3] in the direction of B[c1,c2,c3].
 which is normalized. It returns a scalar."

OVA`perp::usage=
"perp[A[c1,c2,c3],B[c1,c2,c3]] gives the perpendicular
 components of A[c1,c2,c3] with respect to B[c1,c2,c3].
 It returns a vector."

OVA`perpUnitVector::usage=
"perpunitVector[A[c1,c2,c3],B[c1,c2,c3]] gives the perpendicular
 components of A[c1,c2,c3] with respect to B[c1,c2,c3], which is normalized.
 It returns a vector."

OVA`declareVector::usage=
"declareVector[A,B,....] declear A, B,... to be vectors. "

OVA`declareScalar::usage=
"declareScalar[a,b,...] declear a,b,... to be scalars. "

OVA`vector::usage=
"vector[{c1,c2,c3},{c4,c5,c6}] is a vector in a chosen coordinate system."

OVA`vectorQ::usage=
"vectorQ[A] test if A is a vector."

OVA`scalarQ::usage=
"scalarQ[f] test if f is a scalar."


OVA`vectorExpand::usage=
"Expand vector expression to the default canonical form."

OVA`covariantComponent::usage=
"covariantComponent[A,i] gives the ith covariant component of A."

OVA`contravariantComponent::usage=
"contravariantComponent[A,i] gives the ith contravariant component of A."

OVA`about::usage=
"Information about this package"

(**************************** Private ***************************)


Begin["`Private`"]


(* ===================VectorNotation[ ] ============================= *)
vectorNotation[]:=
Module[{},
Unprotect[jacobian] ;
(*===Make Notation for============dot  =============================*)

MakeExpression[RowBox[{a_,"\[CenterDot]",b_}],StandardForm]:=MakeExpression[RowBox[{"dot","[",a,",",b,"]"}],StandardForm] ;

dot/: MakeBoxes[dot[a_,b_],StandardForm]:=RowBox[{MakeBoxes[a,StandardForm],"\[CenterDot]",MakeBoxes[b,StandardForm]}] /; Length[a]==0 && Length[b]==0 ;

dot/: MakeBoxes[dot[a_,b_],StandardForm]:=RowBox[{"(",MakeBoxes[a,StandardForm],")","\[CenterDot]","(",MakeBoxes[b,StandardForm],")" }] /; Length[a]>=1 && Length[b]>=1 ;

dot/: MakeBoxes[dot[a_,b_],StandardForm]:=RowBox[{"(",MakeBoxes[a,StandardForm],")","\[CenterDot]",MakeBoxes[b,StandardForm] }] /; Length[a]>=1 && Length[b]==0 ;

dot/: MakeBoxes[dot[a_,b_],StandardForm]:=RowBox[{MakeBoxes[a,StandardForm],"\[CenterDot]","(",MakeBoxes[b,StandardForm],")" }] /; Length[a]==0 && Length[b]>=1 ;


(*===Make Notation for=============cross ========================*)

MakeExpression[RowBox[{a_,"\[Cross]",b_}],StandardForm]:=MakeExpression[RowBox[{"cross","[",a,",",b,"]"}],StandardForm] ;

cross/: MakeBoxes[cross[a_,b_],StandardForm]:=RowBox[{MakeBoxes[a,StandardForm],"\[Cross]",MakeBoxes[b,StandardForm] }] /; Length[a]==0 && Length[b]==0 ;

cross/: MakeBoxes[cross[a_,b_],StandardForm]:=RowBox[{"(",MakeBoxes[a,StandardForm],")","\[Cross]","(",MakeBoxes[b,StandardForm],")" }] /; Length[a]>=1 && Length[b]>=1 ;

cross/: MakeBoxes[cross[a_,b_],StandardForm]:=RowBox[{"(",MakeBoxes[a,StandardForm],")","\[Cross]",MakeBoxes[b,StandardForm] }] /; Length[a]>=1 && Length[b]==0 ;

cross/: MakeBoxes[cross[a_,b_],StandardForm]:=RowBox[{MakeBoxes[a,StandardForm],"\[Cross]","(",MakeBoxes[b,StandardForm],")" }] /; Length[a]==0 && Length[b]>=1 ;



(*===Make Notation for==============div======================*)

(*MakeExpression[RowBox[{"\[Del]", RowBox[{"\[CenterDot]",a_}]}],StandardForm]:=MakeExpression[RowBox[{"div","[",a,"]"}],StandardForm] ;*)
MakeExpression[RowBox[{"\[Del]","\[CenterDot]",a_}],StandardForm]:=MakeExpression[RowBox[{"div","[",a,"]"}],StandardForm] ;
div/: MakeBoxes[div[a_],StandardForm]:=RowBox[{"\[Del]","\[CenterDot]",MakeBoxes[a,StandardForm]}] /; Length[a]==0 ;


div/: MakeBoxes[div[a_],StandardForm]:=RowBox[{"\[Del]","\[CenterDot]","(",MakeBoxes[a,StandardForm],")"}] /; Length[a]>=1 ;


(*===Make Notation for==============curl====================*)
(*MakeExpression[RowBox[{"\[Del]", RowBox[{"\[Cross]",a_}]}],StandardForm]:=MakeExpression[RowBox[{"curl","[",a,"]"}],StandardForm] ;*)
MakeExpression[RowBox[{"\[Del]","\[Cross]",a_}],StandardForm]:=MakeExpression[RowBox[{"curl","[",a,"]"}],StandardForm] ;

curl/: MakeBoxes[curl[a_],StandardForm]:=RowBox[{"\[Del]","\[Cross]",MakeBoxes[a,StandardForm]}] /; Length[a]==0 ;

curl/: MakeBoxes[curl[a_],StandardForm]:=RowBox[{"\[Del]","\[Cross]","(",MakeBoxes[a,StandardForm],")"}] /; Length[a]>=1 ;



(*===Make Notation for==============CovariantDerivative======*)

MakeExpression[RowBox[{RowBox[{"(",RowBox[{a_, "\[CenterDot]", "\[Del]"}], ")"}], b_}],StandardForm]:=
MakeExpression[RowBox[{"covariantDerivative","[",a, ",",b,"]"}],StandardForm] ;

covariantDerivative/: MakeBoxes[covariantDerivative[a_,b_],StandardForm]:=RowBox[{"(",MakeBoxes[a,StandardForm],"\[CenterDot]", "\[Del]",")",MakeBoxes[b,StandardForm]}] /;Length[a]==0 && Length[b]==0 ;

covariantDerivative/: MakeBoxes[covariantDerivative[a_,b_],StandardForm]:=RowBox[{"(","(",MakeBoxes[a,StandardForm],")","\[CenterDot]", "\[Del]",")","(",MakeBoxes[b,StandardForm],")"}] /; Length[a]>=1 && Length[b]>=1 ;

covariantDerivative/: MakeBoxes[covariantDerivative[a_,b_],StandardForm]:=RowBox[{"(","(",MakeBoxes[a,StandardForm],")","\[CenterDot]", "\[Del]",")",MakeBoxes[b,StandardForm]}] /; Length[a]>=1 && Length[b]==0 ;


covariantDerivative/: MakeBoxes[covariantDerivative[a_,b_],StandardForm]:=RowBox[{"(",MakeBoxes[a,StandardForm],"\[CenterDot]", "\[Del]",")","(",MakeBoxes[b,StandardForm],")"}] /; Length[a]==0 && Length[b]>=1 ;





(*===Make Notation for==============laplacian================*)

MakeExpression[RowBox[{SuperscriptBox["\[Del]", "2"], f_}],StandardForm]:=MakeExpression[RowBox[{"laplacian","[",f,"]"}],StandardForm] ;

laplacian/: MakeBoxes[laplacian[f_],StandardForm]:=RowBox[{SuperscriptBox["\[Del]", "2"], MakeBoxes[f,StandardForm]}] /; Length[f]==0 ;

laplacian/: MakeBoxes[laplacian[f_],StandardForm]:=RowBox[{SuperscriptBox["\[Del]", "2"],"(", MakeBoxes[f,StandardForm],")"}] /; Length[f]>=1 ;


(*===Make Notation for==============grad====================*)
MakeExpression[RowBox[{"\[Del]",a_ }],StandardForm]:=MakeExpression[RowBox[{"grad","[",a,"]"}],StandardForm]  (* /; scalarQ[ToExpression[a]]*) ;

grad/: MakeBoxes[grad[a_],StandardForm]:=RowBox[{"\[Del]",MakeBoxes[a,StandardForm]}] /; Length[a]==0 ;

grad/: MakeBoxes[grad[a_],StandardForm]:=RowBox[{"\[Del]","(",MakeBoxes[a,StandardForm],")"}] /; Length[a]>=1 ;


(*===Make Notation for==============Superscript==============*)

MakeExpression[SuperscriptBox[a_,i_],StandardForm]:=MakeExpression[RowBox[{"Superscript","[",a,",",i,"]"}],StandardForm] ;

Superscript/:MakeBoxes[Superscript[a_,i_],StandardForm]:=SuperscriptBox[MakeBoxes[a,StandardForm],MakeBoxes[i,StandardForm]] /; Length[a]==0 && Length[i]==0 ;

Superscript/:MakeBoxes[Superscript[a_,i_],StandardForm]:=SuperscriptBox[RowBox[{"(",MakeBoxes[a,StandardForm],")"}],RowBox[{"(",MakeBoxes[i,StandardForm],")"}] ] /; Length[a]>=1 && Length[i]>=1 ;

Superscript/:MakeBoxes[Superscript[a_,i_],StandardForm]:=SuperscriptBox[MakeBoxes[a,StandardForm],RowBox[{"(",MakeBoxes[i,StandardForm],")"}] ] /; Length[a]==0&& Length[i]>=1 ;

Superscript/:MakeBoxes[Superscript[a_,i_],StandardForm]:=SuperscriptBox[RowBox[{"(",MakeBoxes[a,StandardForm],")"}],MakeBoxes[i,StandardForm] ] /; Length[a]>=1 && Length[i]==0 ;


(*===Make Notation for==============Subscript==============*)

MakeExpression[SubscriptBox[a_,i_],StandardForm]:=MakeExpression[RowBox[{"Subscript","[",a,",",i,"]"}],StandardForm] ;

Subscript/:MakeBoxes[Subscript[a_,i_],StandardForm]:=SubscriptBox[MakeBoxes[a,StandardForm],MakeBoxes[i,StandardForm]] /; Length[a]==0 && Length[i]==0 ;

Subscript/:MakeBoxes[Subscript[a_,i_],StandardForm]:=SubscriptBox[RowBox[{"(",MakeBoxes[a,StandardForm],")"}],RowBox[{"(",MakeBoxes[i,StandardForm],")"}] ] /; Length[a]>=1 && Length[i]>=1 ;

Subscript/:MakeBoxes[Subscript[a_,i_],StandardForm]:=SubscriptBox[MakeBoxes[a,StandardForm],RowBox[{"(",MakeBoxes[i,StandardForm],")"}] ] /; Length[a]==0&& Length[i]>=1 ;

Subscript/:MakeBoxes[Subscript[a_,i_],StandardForm]:=SubscriptBox[RowBox[{"(",MakeBoxes[a,StandardForm],")"}],MakeBoxes[i,StandardForm] ] /; Length[a]>=1 && Length[i]==0 ;


(*===Make Notation for==============unitVector==============*)

MakeExpression[OverscriptBox[a_, "^"],StandardForm]:=MakeExpression[RowBox[{"unitVector","[",a,"]"}],StandardForm] ;

unitVector/: MakeBoxes[unitVector[a_],StandardForm]:=OverscriptBox[MakeBoxes[a,StandardForm], "^"] /; Length[a]==0 ;

unitVector/: MakeBoxes[unitVector[a_],StandardForm]:=OverscriptBox[RowBox[{"(",MakeBoxes[a,StandardForm],")"}], "^"] /; Length[a]>=1 ;

(*===Make Notation for ============AbsoluteValue===========*)

MakeExpression[RowBox[{"|",a__,"|"}],StandardForm]:=MakeExpression[RowBox[{"absoluteValue","[",a,"]"}],StandardForm] ;

absoluteValue/: MakeBoxes[absoluteValue[a_],StandardForm]:=RowBox[{"|",MakeBoxes[a,StandardForm],"|"}] ;



(*===Make Notation for============= jacobian ==============*)

MakeExpression[RowBox[{"\[ScriptCapitalJ]", "[",a__,"]"}],StandardForm]:=MakeExpression[RowBox[{"jacobian","[",a,"]"}],StandardForm] ;

jacobian/: MakeBoxes[jacobian[a_,b_,c_],StandardForm]:=RowBox[{"\[ScriptCapitalJ]", "[",MakeBoxes[a,StandardForm],",",MakeBoxes[b,StandardForm],",",MakeBoxes[c,StandardForm],"]"}] ;

]


(*============= Open the input palette ================*)
about[]:=Module[{},
Print["OVA packege(Vector Analysis in Orthogonal Coordinates System)"];
Print["Version 1.0"];
Print["This is a modified package based on Prof. Qin's GVA(General Vector Analysis)
 package Version 1.0 in 1997. The detailed copyright of Qin's GVA Version 1.0 are
 shown in OVA.m in its original state."];
Print["
In the OVA, we only perform the symbolic vector analysis in the orthogonal 
coordinate system. All functions in OVA` context start with lowercase letters 
to avoid conflicts with symbols in System`.
"];
Print["author: ymma98"];
Print["e-mail: ymma98@qq.com"];
];


(* ===================SetCoordinateSystem for "None" ================= *)


NotebookPut[Notebook[{
Cell[BoxData[{GridBox[{
        {
          ButtonBox[\(\[SelectionPlaceholder]\[CenterDot]\[Placeholder]\)], 
          ButtonBox[\(\[SelectionPlaceholder]\[Cross]\[Placeholder]\)], 
          
          ButtonBox[
            \(\(\[SelectionPlaceholder]\[CenterDot]\) 
              \((\[Placeholder]\[Cross]\[Placeholder])\)\)], 
          ButtonBox[\(\[Del]\[SelectionPlaceholder]\)]},
        {
          ButtonBox[\(\[Del]\(\[CenterDot]\[SelectionPlaceholder]\)\)], 
          ButtonBox[\(\[Del]\(\[Cross]\[Placeholder]\)\)], 
          
          ButtonBox[
            \(\((\(\[SelectionPlaceholder]\[CenterDot]\) \[Del])\) 
              \[Placeholder]\)], 
          ButtonBox[\(\[Del]\^2 \[SelectionPlaceholder]\)]}
        },
      RowSpacings->0,
      ColumnSpacings->0,
      GridDefaultElement:>ButtonBox[ "\\[Placeholder]"]], 
    ButtonBox[\(\ \ QuickInput\ \ \ \ \),
      ButtonStyle->"Evaluate"]}], NotebookDefault,
  CellMargins->{{Inherited, Inherited}, {5, Inherited}},
  Evaluatable->True,
  CellGroupingRules->"InputGrouping",
  CellHorizontalScrolling->True,
  PageBreakAbove->True,
  PageBreakWithin->False,
  GroupPageBreakWithin->False,
  CellLabelMargins->{{11, Inherited}, {Inherited, Inherited}},
  DefaultFormatType->DefaultInputFormatType,
  LineSpacing->{1.25, 0},
  AutoItalicWords->{},
  FormatType->InputForm,
  ScriptMinSize->9,
  ShowStringCharacters->True,
  NumberMarks->True,
  CounterIncrements->"Input",
  StyleMenuListing->None,
  FontFamily->"Courier",
  FontWeight->"Bold"]
},
FrontEndVersion->"X 3.0",
ScreenRectangle->{{0, 1280}, {0, 1024}},
Editable->False,
WindowToolbars->{},
PageWidth->387,
WindowSize->{180, 80},
WindowMargins->{{58, Automatic}, {Automatic, 404}},
WindowFrame->"Palette",
WindowElements->{},
WindowFrameElements->"CloseBox",
WindowClickSelect->False,
ScrollingOptions->{"PagewiseScrolling"->True},
ShowCellBracket->False,
CellMargins->{{0, 0}, {Inherited, 0}},
Active->True,
CellOpen->True,
ShowCellLabel->False,
ShowCellTags->False,
ImageMargins->{{0, Inherited}, {Inherited, 0}},
Magnification->1
]
];


setCoordinateSystem["None"]:=
Module[{},
(*Unprotect[vectorQ];*)
Unprotect[metric, jacobian, christoffel,
defineVector, dot, cross, grad, div, curl, covariantDerivative, laplacian, absoluteValue,
unitVector, parallel, parallelUnitVector, perp, perpUnitVector, declareVector,
declareScalar, vector, vectorQ, scalarQ, vectorExpand, covariantComponent,
contravariantComponent];
(* Unprotect[curl, div, grad]; *)

ClearAll[vectorQ, scalarQ, dot, cross, grad, div, curl, covariantDerivative];

$coordinateSystem="None";

vectorNotation[];

vectorQ[A_+B_]:=True /; (vectorQ[A] && vectorQ[B]) ;
vectorQ[(x_/; (NumberQ[x] || scalarQ[x])) A_]:=True /; (vectorQ[A]) ;

vectorQ[cross[A_,B_]]:=True /; (vectorQ[A] && vectorQ[B]) ;
vectorQ[curl[A_]]:=True /;(vectorQ[A])  ;
vectorQ[grad[f_]]:=True /;(scalarQ[f])  ;
vectorQ[covariantDerivative[A_,B_]]:=True /; (vectorQ[A] && vectorQ[B]) ;

scalarQ[div[A_]]:=True /;(vectorQ[A]) ;
scalarQ[dot[A_,B_]]:=True /; (vectorQ[A] && vectorQ[B]) ;
scalarQ[f_+g_]:=True /; (scalarQ[f] && scalarQ[g]) ;
scalarQ[(x_/; NumberQ[x]) f_]:=True /; (scalarQ[f]) ;

(*==========Subscript[a,i] and Superscript[a,i] ==========*)
Superscript[a_,i_]:=Power[a,i];


declareVector[a__]:=
Module[{tmp},
tmp={a};
Unprotect[vectorQ];
Do[vectorQ[tmp[[n]]]=True,{n,Length[tmp]}];
Protect[vectorQ];
Print[tmp," are Vectors now."]
];

declareScalar[a__]:=
Module[{tmp},
tmp={a};
Unprotect[scalarQ];
Do[scalarQ[tmp[[n]]]=True,{n,Length[tmp]}];
Protect[scalarQ];
Print[tmp, " are Scalars now."]
];

dot[A_,cross[A_,C_]]:=0 /; (vectorQ[A] && vectorQ[C]) ;

dot[A_,cross[B_,C_]]:=
Module[{tmp},
	tmp=Sort[{A,B,C}];
	Signature[{A,B,C}]dot[tmp[[1]],cross[tmp[[2]],tmp[[3]]]]
] /;  (vectorQ[A] && vectorQ[B] && vectorQ[C] && OrderedQ[{A,B,C}]!=True ) ;


cross[A_,A_]:=0 /; (vectorQ[A]) ;

HoldPattern[curl[grad[f_]]]:=0 /; scalarQ[f] ;

HoldPattern[div[curl[A_]]]:=0 /; vectorQ[A] ;


(*==============  Pull the constant out ========================== *)
HoldPattern[cross[(x_/; NumberQ[x]) A_, B_] ]:= 
x  cross[A,B] /; (vectorQ[A] && vectorQ[B]) ;

HoldPattern[cross[A_, (x_/; NumberQ[x]) B_] ]:= 
x  cross[A,B] /; (vectorQ[A] && vectorQ[B]) ;

HoldPattern[dot[(x_/; NumberQ[x]) A_, B_] ]:= 
x dot[A,B]      /; (vectorQ[A] && vectorQ[B]) ;

HoldPattern[dot[A_, (x_/; NumberQ[x]) B_] ]:= 
x dot[A,B]      /; (vectorQ[A] && vectorQ[B]) ;

HoldPattern[grad[(x_/; NumberQ[x]) f_]]:=x grad[f] /;(scalarQ[f]) ;

HoldPattern[div[(x_/; NumberQ[x]) A_]]:=x div[A] /;(vectorQ[A]) ;

HoldPattern[curl[(x_/; NumberQ[x]) A_]]:=x curl[A] /;(vectorQ[A]) ;

HoldPattern[CovariantDerivative[(x_/; NumberQ[x]) A_,B_]]:=
x CovariantDerivative[A,B] /;(vectorQ[A] && vectorQ[B]);

HoldPattern[CovariantDerivative[A_,(x_/; NumberQ[x]) B_]]:=
x CovariantDerivative[A,B] /;(vectorQ[A] && vectorQ[B]);

SetAttributes[dot, Orderless];

cross[A_,B_]:=-cross[B,A] /; (vectorQ[A] && vectorQ[B] && Order[A,B]==-1);


vectorExpandDispatch=Dispatch[{

(*================== 4th product ==================================*)

HoldPattern[dot[cross[A_,B_],cross[C_,D_]]] :> 
dot[A,C] dot[B,D]-dot[A,D] dot[B,C] /;
(vectorQ[A] && vectorQ[B] && vectorQ[C] && vectorQ[D] ),


(*================== Tripple product ==============================*)
HoldPattern[cross[A_,cross[B_,C_]]]:>dot[A,C] B-dot[A,B] C /; 
(vectorQ[A] && vectorQ[B] && vectorQ[C])  ,


(*================== Differential rules ==================================*)
HoldPattern[grad[f_ g_]] :> f grad[g] + g grad[f] /; (scalarQ[f] && scalarQ[g]),

HoldPattern[div[f_ A_]] :> f div[A]+ dot[A,grad[f]] /; (scalarQ[f] && vectorQ[A]),

HoldPattern[div[A_ f_]] :> f div[A]+ dot[A,grad[f]] /; (scalarQ[f] && vectorQ[A]),

HoldPattern[curl[f_ A_]] :> f curl[A]+ cross[grad[f],A] /; (scalarQ[f] && vectorQ[A]),

HoldPattern[curl[A_ f_]] :> f curl[A]+ cross[grad[f],A] /; (scalarQ[f] && vectorQ[A]),

HoldPattern[div[cross[A_,B_]]] :> dot[B,curl[A]]-dot[A,curl[B]] 
/; (vectorQ[A] && vectorQ[B]),

HoldPattern[curl[cross[A_,B_]]] :> A div[B] - B div[A] + covariantDerivative[B,A]-
covariantDerivative[A,B ] /; (vectorQ[A] && vectorQ[B]),

HoldPattern[grad[dot[A_,B_]]] :> cross[A,curl[B]]+ cross[B,curl[A]] +
covariantDerivative[A,B] +covariantDerivative[B,A] /; (vectorQ[A] && vectorQ[B]),

HoldPattern[div[grad[f_]]] :> laplacian[f] /; scalarQ[f],




(* ================== Bilinearity ==================================*)
HoldPattern[cross[A_+B_,C_]]:>cross[A,C]+cross[B,C] /; 
(vectorQ[A] && vectorQ[B] && vectorQ[C]),

HoldPattern[cross[A_,B_+C_]]:>cross[A,C]+cross[A,B] /; 
(vectorQ[A] && vectorQ[B] && vectorQ[C]),

HoldPattern[dot[A_+B_,C_]]:>dot[A,C]+dot[B,C] /; 
(vectorQ[A] && vectorQ[B] && vectorQ[C]),

HoldPattern[dot[A_,B_+C_]]:>dot[A,C]+dot[A,B] /; 
(vectorQ[A] && vectorQ[B] && vectorQ[C]),

HoldPattern[grad[f_+g_]]:>grad[f]+grad[g] /; (scalarQ[f] && scalarQ[g]),

HoldPattern[div[A_+B_]]:>div[A]+div[B] /; (vectorQ[A] && vectorQ[B]),

HoldPattern[curl[A_+B_]]:>curl[A]+curl[B] /; (vectorQ[A] && vectorQ[B])


}];


vectorExpand[a_]:=
Module[{tmp},
tmp=a //. vectorExpandDispatch;
Expand[tmp]
];

$coordinateSystem;

Protect[metric, jacobian, christoffel,
defineVector, dot, cross, grad, div, curl, covariantDerivative, laplacian, absoluteValue,
unitVector, parallel, parallelUnitVector, perp, perpUnitVector, declareVector,
declareScalar, vector, vectorQ, scalarQ, vectorExpand, covariantComponent,
contravariantComponent];
]


(*=============SetCoordinateSysterm[ ][ ] for a special coordinate system ==============*)

setCoordinateSystem[CoordinateSystem__][x_,y_,z_]:=
Module[{},

(*Unprotect[vectorQ,jacobian];*)
Unprotect[metric, jacobian, christoffel,u,
defineVector, dot, cross, grad, div, curl, covariantDerivative, laplacian, absoluteValue,
unitVector, parallel, parallelUnitVector, perp, perpUnitVector, declareVector,
declareScalar, vector, vectorQ, scalarQ, vectorExpand, covariantComponent,
contravariantComponent];
(* Unprotect[curl, div, grad]; *)

ClearAll[christoffel,u,metric,jacobian];
ClearAll[vectorQ, scalarQ, dot, cross, grad, div, curl, covariantDerivative];

$coordinateSystem=CoordinateSystem;

vectorNotation[];

vectorQ[a_]:= (Head[a]===vector) ;



(* ============ Set up coordinates =================*)
u[1]=x;
u[2]=y;
u[3]=z;

(* ============ Metric Matrix for e^{i}, Metric[2] ===== *)

metric[2][r_,t_,p_]=metrictmp[CoordinateSystem][r,t,p];


(* ============ Metric Matrix for e_{i}, Metric[1] ============ *)
metric[1][r_,t_,p_]=Inverse[metric[2][r,t,p]];


(* ========== jacobian==================== *)
jacobian[r_,t_,p_]=Simplify[PowerExpand[Sqrt[Det[metric[1][r,t,p]]]]];



(* ========== 1st and 2nd kind Christoffel symbol ============*)

christoffel[1][j_,i_,k_,r_,t_,p_]:=
Module[{v},
 v[1]=r;
 v[2]=t;
 v[3]=p;
 1/2(D[metric[1][r,t,p][[j,i]],v[k]]+
     D[metric[1][r,t,p][[j,k]],v[i]]-
     D[metric[1][r,t,p][[i,k]],v[j]])
]; 

christoffel[2][i_,j_,k_,r_,t_,p_]:=
(
Sum[metric[2][r,t,p][[i,n]] christoffel[1][n,j,k,r,t,p],{n,1,3}]
);

christoffel[1][r_,t_,p_]=
             Table[Simplify[christoffel[1][i,j,k,r,t,p],Trig->False],{i,3},{j,3},{k,3}];
christoffel[2][r_,t_,p_]=
             Table[Simplify[christoffel[2][i,j,k,r,t,p],Trig->False],{i,3},{j,3},{k,3}];

Print["============The coordinate system is set up============="];
Print[" The Metric Matrix g^{i,j} is:"];
Print[ metric[2][u[1],u[2],u[3]]  ];
Print[ "======================================================="];
Print[" The Metric Matrix g_{i,j} is:"];
Print[ metric[1][u[1],u[2],u[3]] ];
Print[ "======================================================="];
Print[" The 1st kind Christoffel symbol matrix is:"];
Print[  christoffel[1][u[1],u[2],u[3]]  ];
Print[ "======================================================="];
Print[" The 2st kind Christoffel symbol matrix is:"];
Print[ christoffel[2][u[1],u[2],u[3]]  ];
Print[ "======================================================="];

(* ========== Levi-Civita symbol ====================*)

lc[i_,j_,k_]:=Signature[{i,j,k}];

(* ========= reverse between 1 and 2 =================*)
rev[i_ /; i==2]:=1 ;
rev[i_ /; i==1]:=2 ;

(* ======================================== OVA ======================================== *)
(* ========== Lame coefficients {h1, h2, h3} ========== *)
lame=Table[Sqrt[ metric[1][u[1],u[2],u[3]][[i,i]] ],{i,1,3}];

(* =========== Creat a Vector ============ *)
defineVector[f1_,f2_,f3_]:=
Module[{tmp},
tmp={f1,f2,f3};
vector[tmp]
];

(* ========== Times a constant and Plus ===============*)
vector/: Plus[vector[a__],vector[b__]]:=vector@@Plus[List@@vector[a]+List@@vector[b]];

vector/: Times[q_,vector[a__]] := vector@@Times[q, List@@vector[a] ];

vector/: D[vector[a__],b_]:=vector@@D[List@@vector[a],b];

(* =========== dot[a,b] ============== *)
dot[a_,b_]:=Sum[a[[1,n]] b[[1,n]],{n,1,3}] ;

(* =========== cross[a,b] ======== *)
cross[a_,b_]:= 
Module[{tmp1,tmp},
tmp1 = Table[Sum[lc[i,j,k] a[[1,i]] / lame[[i]] b[[1,j]] / lame[[j]],{i,1,3}, {j,1,3}],{k,1,3}];
tmp = Table[tmp1[[i]] /lame[[i]],{i,3}];
vector[tmp]
];


(* =========== grad[f]===========*)
grad[f_]:=
Module[{tmp1,  tmp},
If[vectorQ[f],
(
tmp1=Table[D[f[[1,j]] lame[[j]],u[i]]-Sum[christoffel[2][u[1],u[2],u[3]][[k,j,i]] f[[1,k]] lame[[k]],{k,3}],{i,3},{j,3}];
tmp=Table[tmp1[[i,j]]/lame[[i]]/lame[[j]],{i,3},{j,3}];
Return[tmp]
)
,
(
tmp1=Table[D[f,u[i]],{i,3}];
tmp = Table[tmp1[[i]] /lame[[i]],{i,3}];
Return[vector[tmp]]
)
]
];


(* =========== div[a]============*)
div[a_]:=
Module[{},
Sum[D[jacobian[u[1],u[2],u[3]] a[[1,i]] / lame[[i]], u[i]],{i,1,3}] / jacobian[u[1],u[2],u[3]]
];


(* ============ curl[a] ==========*)
curl[a_]:=
Module[{tmp2, tmp},
tmp2 = Table[Sum[lc[i,j,k] D[a[[1,j]] lame[[j]], u[i]] ,{i,1,3},{j,1,3}]
                                           ,{k,3}]/jacobian[u[1],u[2],u[3]];
tmp = Table[tmp2[[i]] lame[[i]],{i,3}];
vector[tmp]
];


(* ============ CovariantDerivative[a,b]========*)
covariantDerivative[a_,b_]:=
Module[{tmp, tmp2},
If[vectorQ[b],
(tmp2=Table[Sum[a[[1,i]] / lame[[i]] D[b[[1,j]] / lame[[j]],u[i]],{i,1,3}]+
   Sum[a[[1,i]]/lame[[i]] b[[1,k]]/lame[[k]] christoffel[2][u[1],u[2],u[3]][[j,k,i]]
    ,{i,1,3},{k,1,3}],{j,3}];
    tmp = Table[tmp2[[i]] lame[[i]],{i,3}];
    Return[ vector[tmp] ]
    ),
   Return [ Sum[a[[1,i]] / lame[[i]] D[b,u[i]] ,{i,1,3}] ]
]
];


(* =========== laplacian[a] ==============*)
laplacian[a_]:=
If[MatrixQ[a],
   Return[grad[div[a]]-curl[curl[a]]],
   Return[div[grad[a]]]
  ] ;
(* =========== AbsoluteValue[a] ===========*)
absoluteValue[a_]:=
PowerExpand[Sqrt[dot[a,a]]] ;

(* =========== unitVector[a] ==============*)
unitVector[a_]:=
a/PowerExpand[Sqrt[dot[a,a]]] ;

(* =========== Parallel[a,b] ============*)
parallel[a_,b_]:= dot[a,unitVector[b]] ;

parallelunitVector[a_,b_]:=dot[a,b];

(* =========== Perp[a,b] ================*)
perp[a_,b_]:=
Module[{tmp1,tmp2},
tmp1=unitVector[b];
 Simplify[cross[cross[tmp1,a],tmp1]]
] ;

perpUnitVector[a_,b_]:=
cross[cross[b,a],b] ;

(* ========== Subscript[a,i]  and CovariantComponent[a,i] =================*)
Subscript[a_,i_]:=covariantComponent[a,i] /; (vectorQ[a] && IntegerQ[i]) ;
covariantComponent[a_,i_]:=a[[1,i]] /; (vectorQ[a] && IntegerQ[i]);

(* ========== Superscript[a,i] and ContravariantComponent[a,i]================*)

Superscript[a_,i_]:=contravariantComponent[a,i] /; (vectorQ[a] && IntegerQ[i]) ;
Superscript[a_,i_]:=Power[a,i] /; !(vectorQ[a] && IntegerQ[i]);
contravariantComponent[a_,i_]:=a[[1,i]] /; (vectorQ[a]==True && IntegerQ[i]==True);

Print[{CoordinateSystem}, " is set up."];

Protect[metric, jacobian, christoffel,u,
defineVector, dot, cross, grad, div, curl, covariantDerivative, laplacian, absoluteValue,
unitVector, parallel, parallelUnitVector, perp, perpUnitVector, declareVector,
declareScalar, vector, vectorQ, scalarQ, vectorExpand, covariantComponent,
contravariantComponent, about];

]




(* =================DefineCoordinateSystem ========================*)
defineCoordinateSystem[CoordinateSystem__][c1_,c2_,c3_][matrix__]:=
Module[{r,t,p},
metrictmp[CoordinateSystem][r_,t_,p_]=(matrix /. {c1->r,c2->t,c3->p});
Print[{CoordinateSystem}," is added to current Mathematica session."];
Print["The Riemann Metric Tensor is:"];
Print["============================="];
Print[metrictmp[CoordinateSystem][r,t,p]];
Print["Use SetCoordinateSystem to set up a coordinate system"];
]


(* ================= Coordinate Systems ========================*)

(* ----------------Cylindrical ------------------------ *)
metrictmp["Cylindrical"][r_,t_,z_]:={{1,0,0},{0,1/r^2,0},{0,0,1}}
 

(* ----------------Cartesian---------------------------*)
metrictmp["Cartesian"][x_,y_,z_]:={{1,0,0},{0,1,0},{0,0,1}}


(* ----------------Spherical --------------------------*)
metrictmp["Spherical"][r_,t_,p_]:={{1,0,0},{0,1/r^2,0},
                               {0,0,1/(r^2 Sin[t]^2)}}

(* --------------STKCircular-----------------------------*)
(* --------------No small parameter here     ------------*)
metrictmp["STKCircular",R0_,ep_][r_,t_,p_]:=
(
{
 {
  1,
  0,
  0
 },
 {
  0,
  1/r^2,
  0
 },
 {
  0,
  0,
  1/R0^2
 }
}
)

(*----------------TKGeneralPerp------------------------*)
(*----------------TKGeneralPerp------------------------*)
metrictmp["TKGeneralPerp",g1_, g2_,g3_][r_,t_,p_]:=
(
{
 {
  g1[r,t]^2,
  0,
  0
 },
 {
  0,
  g2[r,t]^2,
  0
 },
 {
  0,
  0,
  g3[r,t]^2
 }
}
)

(* ----------------TKCircularFlux ----------------------*)
(* --------------To the leading order of toroidicity ---*)

metrictmp["TKCircularFlux",R0_,ept_][r_,t_,p_]:=
(
{
 {
  1,
  -(r/R0 ept  Sin[t])/r+Series[ept^2,{ept,0,1}],
  0
 },
 {
  -(r/R0 ept  Sin[t])/r+Series[ept^2,{ept,0,1}],
  (1-2r/R0 ept  Cos[t])/r^2+Series[ept^2,{ept,0,1}],
  0
 },
 {
  0,
  0,
  (1-2r/R0 ept Cos[t])/R0^2 + Series[ept^2,{ept,0,1}]
 }
}
)


(* ----------------TKCircularFlux2----------------------*)
(* --------------To the 2nd order of toroidicity ---*)

metrictmp["TKCircularFlux2",R0_,ept_][r_,t_,p_]:=
(
{
 {
  1,
  -(r/R0 ept  Sin[t])/r+ r t/R0^2 ept^2 +Series[ept^3,{ept,0,2}],
  0
 },
 {
  -(r/R0 ept  Sin[t])/r+ r t/R0^2 ept^2 +Series[ept^3,{ept,0,2}],
  (1-2r/R0 ept  Cos[t])/r^2+ 3/R0^2 ept^2 +Series[ept^3,{ept,0,2}],
  0
 },
 {
  0,
  0,
  (1-2r/R0 ept Cos[t])/R0^2 + r^2/R0^4 (5+Cos[2 t])/2 ept^2 +Series[ept^3,{ept,0,2}]
 }
}
)
(* ---------------TKCircular-----------------------------*)
(* ---------------don't use the ep expansion --------*)
(* --------Thought, ep is kept there for later use ------*)
metrictmp["TKCircular",R0_,ep_][r_,t_,p_]:=
(
{
 {
  1,
  0,
  0
 },
 {
  0,
  1/r^2,
  0
 },
 {
  0,
  0,
  1/R0^2 /(1+r/R0 ep Cos[t])^2
 }
}
)
(* ---------------TKCircular6-----------------------------*)
(* ---------------the odering of toroidicity is 6 --------*)

metrictmp["TKCircular6",R0_,ep_][r_,t_,p_]:=
(
{
 {
  1,
  0,
  0
 },
 {
  0,
  1/r^2,
  0
 },
 {
  0,
  0,
  Series[1/R0^2 /(1+r/R0 ep Cos[t])^2,{ep,0,6}]
 }
}
)
(* ---------------TKCircular5-----------------------------*)
(* ---------------the odering of toroidicity is 5 --------*)

metrictmp["TKCircular5",R0_,ep_][r_,t_,p_]:=
(
{
 {
  1,
  0,
  0
 },
 {
  0,
  1/r^2,
  0
 },
 {
  0,
  0,
  Series[1/R0^2 /(1+r/R0 ep Cos[t])^2,{ep,0,5}]
 }
}
)
(* ---------------TKCircular4-----------------------------*)
(* ---------------the odering of toroidicity is 4 --------*)

metrictmp["TKCircular4",R0_,ep_][r_,t_,p_]:=
(
{
 {
  1,
  0,
  0
 },
 {
  0,
  1/r^2,
  0
 },
 {
  0,
  0,
  Series[1/R0^2 /(1+r/R0 ep Cos[t])^2,{ep,0,4}]
 }
}
)
(* ---------------TKShiftedCircular4-----------------------------*)
(* ---------------the odering of toroidicity is 4 --------*)

metrictmp["TKShiftedCircular4",R0_,ep_,R0dp_][r_,t_,p_]:=
(
{
 {
  1+2 R0dp[r] Cos[t]/R0 ep,
  -R0dp[r] Sin[t] /r /R0 ep,
  0
 },
 {
  -R0dp[r] Sin[t] /r /R0 ep,
  1/r^2,
  0
 },
 {
  0,
  0,
  Series[1/R0^2 /(1+r/R0 ep Cos[t])^2,{ep,0,4}]
 }
}
)



(* ---------------TKCircular3-----------------------------*)
(* ---------------the odering of toroidicity is 3 --------*)

metrictmp["TKCircular3",R0_,ep_][r_,t_,p_]:=
(
{
 {
  1,
  0,
  0
 },
 {
  0,
  1/r^2,
  0
 },
 {
  0,
  0,
  Series[1/R0^2 /(1+r/R0 ep Cos[t])^2,{ep,0,3}]
 }
}
)
(* ---------------TKCircular2-----------------------------*)
(* ---------------the odering of toroidicity is 2 --------*)

metrictmp["TKCircular2",R0_,ep_][r_,t_,p_]:=
(
{
 {
  1,
  0,
  0
 },
 {
  0,
  1/r^2,
  0
 },
 {
  0,
  0,
  Series[1/R0^2 /(1+r/R0 ep Cos[t])^2,{ep,0,2}]
 }
}
)




(*--put the metric matrix of your own coordinate system here ----*)



Print["The Vector Calculus on Orthogonal Coordinate is loaded in"]
Print["Use SetCoordinateSystem[ ]  to set up a coordinate system"]
Print["The default CoordinateSystem is  None "]
$coordinateSystem="None"
setCoordinateSystem["None"]


End[]
Protect[setCoordinateSystem, defineCoordinateSystem, about];
EndPackage[]



(*-- End --*)
