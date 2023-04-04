(* ::Package:: *)

                       
(****************************************************************)
(*                                                              *)
(*       Vector Analysis in General Coordinates  System         *)
(*                                                              *)
(*       Version 1.0                                            *)
(*                                                              *)
(****************************************************************)



(****************************************************************)   
(*                                                              *)
(*       Hong Qin                                               *)
(*       Princeton Plasma Physics Laboratory                    *)
(*       Princeton University, Princeton, NJ, 08540             *)
(*       Phone:(609) 243-3653, Fax:(609) 243-2662               *)
(*       Email: HongQin@Princeton.edu                           *)
(*       WWW:http://w3.pppl.gov/~hongqin                        *)
(*                                                              *)
(****************************************************************)



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
(*  (1) Vector analysis in general coordinates                  *)
(*  (2) Define your own coordinates to work on                  *)
(*  (3) Perturbative vector analysis with respect to any small  *)
(*             parameter, e.g. inverse aspect ratio in tokamak, *)
(*             to any order                                     *)
(*  (4) Vector analysis in Shafranov coordinate,                *)
(*             flux coordinate, straight tokamak coordinate,    *)
(*             circular concentric tokamak coordinate, ......   *)  
(*                                                              *)   
(****************************************************************)



(************ Release Notes *************************************)
(*
For a complete report,see the paper,

automatic symbolic vector analysis in general coordinate systems
and its application to plasma physics,

and other related documents 
available at http://w3.lyman.pppl.gov/~hongqin/ .


This package is designed to perform automatic symbolic computation of  
3D vector analysis in general coordinate systems. The general coordinate
system could be any possible mathematically well defined coordinate 
system ( its Jacobian is non-vanishing everywhere.) It includes these
well known orthogonal coordinate systems like Cartesian coordinat,      
cylindrical coordinate, and spherical coordinate,  and also other      
non-orthogonal coordinate systems such as the flux coordinate used    
extensively in thermal fusion devices. 


The modern viewpoint of 3D vector calculus is used to simplify and
unify the algorithm. A general coordinate is defined by the Reimann  
metric matrix  expressed in its own coordinates. A vector field is     
isotopic to a one form ( covariant components) in the 3D manifolds and 
its image under the  Hodger star operator (contravariant components).  
All vector calculus operators can be expressed in the language of       
differential forms in a simple and unified manner. To define a new     
coordinate system of yourself, you only need to define the corresponding 
Reimann metric matrix. 


Users of this package need not to familarize themselves to this
mathematical theory, unless they want to set up new coordinate systems 
of their owns. The only thing need to realize here is the fact that in 
this package a vector is a  2x3 list, the [[1]] part of which is the
covariant component and [[2]] part contravariant component.

*)

(************ Begin Calculus`GeneralVectorAnalysis` Package *****)



BeginPackage["Calculus`GeneralVectorAnalysis`"]

Calculus`GeneralVectorAnalysis`SetCoordinateSystem::usage=
"SetCoordinateSystem[$CoordinateSystem,p1,p2][c1,c2,c3] will set up a 
$CoordinateSystem coordinate with coordinate c1,c2,c3 and parameter p1,p2."

Calculus`GeneralVectorAnalysis`DefineCoordinateSystem::usage=
"DefineCoordinateSystem[$CoordinateSystem,p1,p2,...][c1,c2,c3][MetricMatrix] will
add the $CoordinateSystem with parameter p1,p2,... to the package. The Metric 
tensor as a function of coordinate c1,c2,c3 is given by MetricMatrix."


Calculus`GeneralVectorAnalysis`Metric::usage=
"Metric[1][c1,c2,c3] is the metric matrix g_{i,j},
 Metric[2][c1,c2,c3] is the metric matrix g^{i,j}."

Calculus`GeneralVectorAnalysis`Jacobian::usage="Jacobian[c1,c2,c3] gives the Jacobian.  "

Calculus`GeneralVectorAnalysis`Christoffel::usage=
"GChristoffe1[1][c1,c2,c3] is the 1st kind Christoffel symbol,
 GChristoffe1[2][c1,c2,c3] is the 2nd kind Christoffel symbol."

Calculus`GeneralVectorAnalysis`DefineVector::usage=
"DefineVector[1,A1[c1,c2,c3],A2[c1,c2,c3],A3[c1,c2,c3]]
 creats a vector with A1,A2,A3 as its covariant component; 
 DefineVector[2,A1[c1,c2,c3],A2[c1,c2,c3],A3[c1,c2,c3]]
 creats a vector with A1,A2,A3 as its comtravariant component. "

Calculus`GeneralVectorAnalysis`DotProduct::usage=
"DotProduct[A[c1,c2,c3],B[c1,c2,c3]] gives the dot product. "

Calculus`GeneralVectorAnalysis`CrossProduct::usage=
"CrossProduct[A[c1,c2,c3],B[c1,c2,c3]] gives the cross product."

Calculus`GeneralVectorAnalysis`Grad::usage=
"Grad[f[c1,c2,c3]] gives the gradiant of f."

Calculus`GeneralVectorAnalysis`Div::usage=
"Div[A[c1,c2,c3]] gives the divergence of A."

Calculus`GeneralVectorAnalysis`Curl::usage=
"Curl[A[c1,c2,c3]] gives the curl of A."

Calculus`GeneralVectorAnalysis`CovariantDerivative::usage=
"CovariantDerivative[A[c1,c2,c3],x[c1,c2,c3]] gives the directional
 derivative of x in A direction, 
 x can be either a vector or a scalar."

Calculus`GeneralVectorAnalysis`Laplacian::usage=
"Laplacian[f[c1,c2,c3]] gives the Laplacian of scalar f."

Calculus`GeneralVectorAnalysis`AbsoluteValue::usage=
"AbsoluteValue[A[c1,c2,c3]] gives the absolute
 value of vector A[c1,c2,c3]. "

Calculus`GeneralVectorAnalysis`UnitVector::usage=
"UnitVector[A[c1,c2,c3]] gives the normalized 
 vector in the direction of A[c1,c2,c3]. "

Calculus`GeneralVectorAnalysis`Parallel::usage=
"Parallel[A[c1,c2,c3],B[c1,c2,c3]] gives the 
 projection of A[c1,c2,c3] in the direction of B[c1,c2,c3].
 It returns a scalar."

Calculus`GeneralVectorAnalysis`ParallelUnitVector::usage=
"Parallel[A[c1,c2,c3],B[c1,c2,c3]] gives the 
 the projection of A[c1,c2,c3] in the direction of B[c1,c2,c3].
 which is normalized. It returns a scalar."

Calculus`GeneralVectorAnalysis`Perp::usage=
"Perp[A[c1,c2,c3],B[c1,c2,c3]] gives the perpendicular
 components of A[c1,c2,c3] with respect to B[c1,c2,c3].
 It returns a vector."

Calculus`GeneralVectorAnalysis`PerpUnitVector::usage=
"PerpUnitVector[A[c1,c2,c3],B[c1,c2,c3]] gives the perpendicular
 components of A[c1,c2,c3] with respect to B[c1,c2,c3], which is normalized.
 It returns a vector."

Calculus`GeneralVectorAnalysis`DeclareVector::usage=
"DeclareVector[A,B,....] declear A, B,... to be vectors. "

Calculus`GeneralVectorAnalysis`DeclareScalar::usage=
"DeclareScalar[a,b,...] declear a,b,... to be scalars. "

Calculus`GeneralVectorAnalysis`Vector::usage=
"Vector[{c1,c2,c3},{c4,c5,c6}] is a vector in a chosen coordinate system."

Calculus`GeneralVectorAnalysis`VectorQ::usage=
"VectorQ[A] test if A is a vector."

Calculus`GeneralVectorAnalysis`ScalarQ::usage=
"ScalarQ[f] test if f is a scalar."


Calculus`GeneralVectorAnalysis`VectorExpand::usage=
"Expand vector expression to the default canonical form."

Calculus`GeneralVectorAnalysis`CovariantComponent::usage=
"CovariantComponent[A,i] gives the ith covariant component of A."

Calculus`GeneralVectorAnalysis`ContravariantComponent::usage=
"ContravariantComponent[A,i] gives the ith contravariant component of A."

Calculus`GeneralVectorAnalysis`About::usage=
"Information about this package"

(**************************** Private ***************************)





Begin["`Private`"]


(* ===================VectorNotation[ ] ============================= *)
VectorNotation[]:=
Module[{},
Unprotect[Jacobian] ;
(*===Make Notation for============DotProduct  =============================*)

MakeExpression[RowBox[{a_,"\[CenterDot]",b_}],StandardForm]:=MakeExpression[RowBox[{"DotProduct","[",a,",",b,"]"}],StandardForm] ;

DotProduct/: MakeBoxes[DotProduct[a_,b_],StandardForm]:=RowBox[{MakeBoxes[a,StandardForm],"\[CenterDot]",MakeBoxes[b,StandardForm]}] /; Length[a]==0 && Length[b]==0 ;

DotProduct/: MakeBoxes[DotProduct[a_,b_],StandardForm]:=RowBox[{"(",MakeBoxes[a,StandardForm],")","\[CenterDot]","(",MakeBoxes[b,StandardForm],")" }] /; Length[a]>=1 && Length[b]>=1 ;

DotProduct/: MakeBoxes[DotProduct[a_,b_],StandardForm]:=RowBox[{"(",MakeBoxes[a,StandardForm],")","\[CenterDot]",MakeBoxes[b,StandardForm] }] /; Length[a]>=1 && Length[b]==0 ;

DotProduct/: MakeBoxes[DotProduct[a_,b_],StandardForm]:=RowBox[{MakeBoxes[a,StandardForm],"\[CenterDot]","(",MakeBoxes[b,StandardForm],")" }] /; Length[a]==0 && Length[b]>=1 ;


(*===Make Notation for=============CrossProduct ========================*)

MakeExpression[RowBox[{a_,"\[Cross]",b_}],StandardForm]:=MakeExpression[RowBox[{"CrossProduct","[",a,",",b,"]"}],StandardForm] ;

CrossProduct/: MakeBoxes[CrossProduct[a_,b_],StandardForm]:=RowBox[{MakeBoxes[a,StandardForm],"\[Cross]",MakeBoxes[b,StandardForm] }] /; Length[a]==0 && Length[b]==0 ;

CrossProduct/: MakeBoxes[CrossProduct[a_,b_],StandardForm]:=RowBox[{"(",MakeBoxes[a,StandardForm],")","\[Cross]","(",MakeBoxes[b,StandardForm],")" }] /; Length[a]>=1 && Length[b]>=1 ;

CrossProduct/: MakeBoxes[CrossProduct[a_,b_],StandardForm]:=RowBox[{"(",MakeBoxes[a,StandardForm],")","\[Cross]",MakeBoxes[b,StandardForm] }] /; Length[a]>=1 && Length[b]==0 ;

CrossProduct/: MakeBoxes[CrossProduct[a_,b_],StandardForm]:=RowBox[{MakeBoxes[a,StandardForm],"\[Cross]","(",MakeBoxes[b,StandardForm],")" }] /; Length[a]==0 && Length[b]>=1 ;



(*===Make Notation for==============Div======================*)

MakeExpression[RowBox[{"\[Del]", RowBox[{"\[CenterDot]",a_}]}],StandardForm]:=MakeExpression[RowBox[{"Div","[",a,"]"}],StandardForm] ;

Div/: MakeBoxes[Div[a_],StandardForm]:=RowBox[{"\[Del]","\[CenterDot]",MakeBoxes[a,StandardForm]}] /; Length[a]==0 ;


Div/: MakeBoxes[Div[a_],StandardForm]:=RowBox[{"\[Del]","\[CenterDot]","(",MakeBoxes[a,StandardForm],")"}] /; Length[a]>=1 ;


(*===Make Notation for==============Curl====================*)
MakeExpression[RowBox[{"\[Del]", RowBox[{"\[Cross]",a_}]}],StandardForm]:=MakeExpression[RowBox[{"Curl","[",a,"]"}],StandardForm] ;

Curl/: MakeBoxes[Curl[a_],StandardForm]:=RowBox[{"\[Del]","\[Cross]",MakeBoxes[a,StandardForm]}] /; Length[a]==0 ;

Curl/: MakeBoxes[Curl[a_],StandardForm]:=RowBox[{"\[Del]","\[Cross]","(",MakeBoxes[a,StandardForm],")"}] /; Length[a]>=1 ;



(*===Make Notation for==============CovariantDerivative======*)

MakeExpression[RowBox[{RowBox[{"(",RowBox[{a_, "\[CenterDot]", "\[Del]"}], ")"}], b_}],StandardForm]:=
MakeExpression[RowBox[{"CovariantDerivative","[",a, ",",b,"]"}],StandardForm] ;

CovariantDerivative/: MakeBoxes[CovariantDerivative[a_,b_],StandardForm]:=RowBox[{"(",MakeBoxes[a,StandardForm],"\[CenterDot]", "\[Del]",")",MakeBoxes[b,StandardForm]}] /;Length[a]==0 && Length[b]==0 ;

CovariantDerivative/: MakeBoxes[CovariantDerivative[a_,b_],StandardForm]:=RowBox[{"(","(",MakeBoxes[a,StandardForm],")","\[CenterDot]", "\[Del]",")","(",MakeBoxes[b,StandardForm],")"}] /; Length[a]>=1 && Length[b]>=1 ;

CovariantDerivative/: MakeBoxes[CovariantDerivative[a_,b_],StandardForm]:=RowBox[{"(","(",MakeBoxes[a,StandardForm],")","\[CenterDot]", "\[Del]",")",MakeBoxes[b,StandardForm]}] /; Length[a]>=1 && Length[b]==0 ;


CovariantDerivative/: MakeBoxes[CovariantDerivative[a_,b_],StandardForm]:=RowBox[{"(",MakeBoxes[a,StandardForm],"\[CenterDot]", "\[Del]",")","(",MakeBoxes[b,StandardForm],")"}] /; Length[a]==0 && Length[b]>=1 ;





(*===Make Notation for==============Laplacian================*)

MakeExpression[RowBox[{SuperscriptBox["\[Del]", "2"], f_}],StandardForm]:=MakeExpression[RowBox[{"Laplacian","[",f,"]"}],StandardForm] ;

Laplacian/: MakeBoxes[Laplacian[f_],StandardForm]:=RowBox[{SuperscriptBox["\[Del]", "2"], MakeBoxes[f,StandardForm]}] /; Length[f]==0 ;

Laplacian/: MakeBoxes[Laplacian[f_],StandardForm]:=RowBox[{SuperscriptBox["\[Del]", "2"],"(", MakeBoxes[f,StandardForm],")"}] /; Length[f]>=1 ;


(*===Make Notation for==============Grad====================*)
MakeExpression[RowBox[{"\[Del]",a_ }],StandardForm]:=MakeExpression[RowBox[{"Grad","[",a,"]"}],StandardForm]  (* /; ScalarQ[ToExpression[a]] *)  ;

Grad/: MakeBoxes[Grad[a_],StandardForm]:=RowBox[{"\[Del]",MakeBoxes[a,StandardForm]}] /; Length[a]==0 ;

Grad/: MakeBoxes[Grad[a_],StandardForm]:=RowBox[{"\[Del]","(",MakeBoxes[a,StandardForm],")"}] /; Length[a]>=1 ;


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


(*===Make Notation for==============UnitVector==============*)

MakeExpression[OverscriptBox[a_, "^"],StandardForm]:=MakeExpression[RowBox[{"UnitVector","[",a,"]"}],StandardForm] ;

UnitVector/: MakeBoxes[UnitVector[a_],StandardForm]:=OverscriptBox[MakeBoxes[a,StandardForm], "^"] /; Length[a]==0 ;

UnitVector/: MakeBoxes[UnitVector[a_],StandardForm]:=OverscriptBox[RowBox[{"(",MakeBoxes[a,StandardForm],")"}], "^"] /; Length[a]>=1 ;

(*===Make Notation for ============AbsoluteValue===========*)

MakeExpression[RowBox[{"|",a__,"|"}],StandardForm]:=MakeExpression[RowBox[{"AbsoluteValue","[",a,"]"}],StandardForm] ;

AbsoluteValue/: MakeBoxes[AbsoluteValue[a_],StandardForm]:=RowBox[{"|",MakeBoxes[a,StandardForm],"|"}] ;



(*===Make Notation for============= Jacobian ==============*)

MakeExpression[RowBox[{"\[ScriptCapitalJ]", "[",a__,"]"}],StandardForm]:=MakeExpression[RowBox[{"Jacobian","[",a,"]"}],StandardForm] ;

Jacobian/: MakeBoxes[Jacobian[a_,b_,c_],StandardForm]:=RowBox[{"\[ScriptCapitalJ]", "[",MakeBoxes[a,StandardForm],",",MakeBoxes[b,StandardForm],",",MakeBoxes[c,StandardForm],"]"}] ;

]

(*============= Open the input palette ================*)
About[]:=NotebookPut[NotebookPut[Notebook[{
Cell[BoxData[GridBox[{
        {\(GeneralVectorAnalysis  Version  1.0\)},
        {\(Copy   Right :  Hong  Qin 1997 \)},
        {\(Send Comments  To :  HongQin@Princeton.Edu\)}
        }]], NotebookDefault,
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
PageWidth->528,
WindowSize->{283, 50},
WindowMargins->{{Automatic, 444}, {Automatic, 57}},
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
Magnification->1 ]] ]


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
    ButtonBox[\(\ \ GeneralVectorAnalysis\ \ \ \ \),
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



(* ===================SetCoordinateSystem for "None" ================= *)


SetCoordinateSystem["None"]:=
Module[{},
Unprotect[VectorQ];
(* Unprotect[Curl, Div, Grad]; *)

ClearAll[VectorQ, ScalarQ, DotProduct, CrossProduct, Grad, Div, Curl, CovariantDerivative];

$CoordinateSystem="None";

VectorNotation[];

VectorQ[A_+B_]:=True /; (VectorQ[A] && VectorQ[B]) ;
VectorQ[(x_/; (NumberQ[x] || ScalarQ[x])) A_]:=True /; (VectorQ[A]) ;

VectorQ[CrossProduct[A_,B_]]:=True /; (VectorQ[A] && VectorQ[B]) ;
VectorQ[Curl[A_]]:=True /;(VectorQ[A])  ;
VectorQ[Grad[f_]]:=True /;(ScalarQ[f])  ;
VectorQ[CovariantDerivative[A_,B_]]:=True /; (VectorQ[A] && VectorQ[B]) ;

ScalarQ[Div[A_]]:=True /;(VectorQ[A]) ;
ScalarQ[DotProduct[A_,B_]]:=True /; (VectorQ[A] && VectorQ[B]) ;
ScalarQ[f_+g_]:=True /; (ScalarQ[f] && ScalarQ[g]) ;

ScalarQ[(x_/; NumberQ[x]) f_]:=True /; (ScalarQ[f]) ;

(*==========Subscript[a,i] and Superscript[a,i] ==========*)
Superscript[a_,i_]:=Power[a,i];


DeclareVector[a__]:=
Module[{tmp},
tmp={a};
Do[VectorQ[tmp[[n]]]=True,{n,Length[tmp]}];
Print[tmp," are Vectors now."]
];

DeclareScalar[a__]:=
Module[{tmp},
tmp={a};
Do[ScalarQ[tmp[[n]]]=True,{n,Length[tmp]}];
Print[tmp, " are Scalars now."]
];

DotProduct[A_,CrossProduct[A_,C_]]:=0 /; (VectorQ[A] && VectorQ[C]) ;

DotProduct[A_,CrossProduct[B_,C_]]:=
Module[{tmp},
	tmp=Sort[{A,B,C}];
	Signature[{A,B,C}]DotProduct[tmp[[1]],CrossProduct[tmp[[2]],tmp[[3]]]]
] /;  (VectorQ[A] && VectorQ[B] && VectorQ[C] && OrderedQ[{A,B,C}]!=True ) ;


CrossProduct[A_,A_]:=0 /; (VectorQ[A]) ;

HoldPattern[Curl[Grad[f_]]]:=0 /; ScalarQ[f] ;

HoldPattern[Div[Curl[A_]]]:=0 /; VectorQ[A] ;


(*==============  Pull the constant out ========================== *)
HoldPattern[CrossProduct[(x_/; NumberQ[x]) A_, B_] ]:= 
x  CrossProduct[A,B] /; (VectorQ[A] && VectorQ[B]) ;

HoldPattern[CrossProduct[A_, (x_/; NumberQ[x]) B_] ]:= 
x  CrossProduct[A,B] /; (VectorQ[A] && VectorQ[B]) ;

HoldPattern[DotProduct[(x_/; NumberQ[x]) A_, B_] ]:= 
x DotProduct[A,B]      /; (VectorQ[A] && VectorQ[B]) ;

HoldPattern[DotProduct[A_, (x_/; NumberQ[x]) B_] ]:= 
x DotProduct[A,B]      /; (VectorQ[A] && VectorQ[B]) ;

HoldPattern[Grad[(x_/; NumberQ[x]) f_]]:=x Grad[f] /;(ScalarQ[f]) ;

HoldPattern[Div[(x_/; NumberQ[x]) A_]]:=x Div[A] /;(VectorQ[A]) ;

HoldPattern[Curl[(x_/; NumberQ[x]) A_]]:=x Curl[A] /;(VectorQ[A]) ;

HoldPattern[CovariantDerivative[(x_/; NumberQ[x]) A_,B_]]:=
x CovariantDerivative[A,B] /;(VectorQ[A] && VectorQ[B]);

HoldPattern[CovariantDerivative[A_,(x_/; NumberQ[x]) B_]]:=
x CovariantDerivative[A,B] /;(VectorQ[A] && VectorQ[B]);

SetAttributes[DotProduct, Orderless];

CrossProduct[A_,B_]:=-CrossProduct[B,A] /; (VectorQ[A] && VectorQ[B] && Order[A,B]==-1);


VectorExpandDispatch=Dispatch[{

(*================== 4th product ==================================*)

HoldPattern[DotProduct[CrossProduct[A_,B_],CrossProduct[C_,D_]]] :> 
DotProduct[A,C] DotProduct[B,D]-DotProduct[A,D] DotProduct[B,C] /;
(VectorQ[A] && VectorQ[B] && VectorQ[C] && VectorQ[D] ),


(*================== Tripple product ==============================*)
HoldPattern[CrossProduct[A_,CrossProduct[B_,C_]]]:>DotProduct[A,C] B-DotProduct[A,B] C /; 
(VectorQ[A] && VectorQ[B] && VectorQ[C])  ,


(*================== Differential rules ==================================*)
HoldPattern[Grad[f_ g_]] :> f Grad[g] + g Grad[f] /; (ScalarQ[f] && ScalarQ[g]),

HoldPattern[Div[f_ A_]] :> f Div[A]+ DotProduct[A,Grad[f]] /; (ScalarQ[f] && VectorQ[A]),

HoldPattern[Div[A_ f_]] :> f Div[A]+ DotProduct[A,Grad[f]] /; (ScalarQ[f] && VectorQ[A]),

HoldPattern[Curl[f_ A_]] :> f Curl[A]+ CrossProduct[Grad[f],A] /; (ScalarQ[f] && VectorQ[A]),

HoldPattern[Curl[A_ f_]] :> f Curl[A]+ CrossProduct[Grad[f],A] /; (ScalarQ[f] && VectorQ[A]),

HoldPattern[Div[CrossProduct[A_,B_]]] :> DotProduct[B,Curl[A]]-DotProduct[A,Curl[B]] 
/; (VectorQ[A] && VectorQ[B]),

HoldPattern[Curl[CrossProduct[A_,B_]]] :> A Div[B] - B Div[A] + CovariantDerivative[B,A]-
CovariantDerivative[A,B ] /; (VectorQ[A] && VectorQ[B]),

HoldPattern[Grad[DotProduct[A_,B_]]] :> CrossProduct[A,Curl[B]]+ CrossProduct[B,Curl[A]] +
CovariantDerivative[A,B] +CovariantDerivative[B,A] /; (VectorQ[A] && VectorQ[B]),

HoldPattern[Div[Grad[f_]]] :> Laplacian[f] /; ScalarQ[f],




(* ================== Bilinearity ==================================*)
HoldPattern[CrossProduct[A_+B_,C_]]:>CrossProduct[A,C]+CrossProduct[B,C] /; 
(VectorQ[A] && VectorQ[B] && VectorQ[C]),

HoldPattern[CrossProduct[A_,B_+C_]]:>CrossProduct[A,C]+CrossProduct[A,B] /; 
(VectorQ[A] && VectorQ[B] && VectorQ[C]),

HoldPattern[DotProduct[A_+B_,C_]]:>DotProduct[A,C]+DotProduct[B,C] /; 
(VectorQ[A] && VectorQ[B] && VectorQ[C]),

HoldPattern[DotProduct[A_,B_+C_]]:>DotProduct[A,C]+DotProduct[A,B] /; 
(VectorQ[A] && VectorQ[B] && VectorQ[C]),

HoldPattern[Grad[f_+g_]]:>Grad[f]+Grad[g] /; (ScalarQ[f] && ScalarQ[g]),

HoldPattern[Div[A_+B_]]:>Div[A]+Div[B] /; (VectorQ[A] && VectorQ[B]),

HoldPattern[Curl[A_+B_]]:>Curl[A]+Curl[B] /; (VectorQ[A] && VectorQ[B])


}];


VectorExpand[a_]:=
Module[{tmp},
tmp=a //. VectorExpandDispatch;
Expand[tmp]
];

$CoordinateSystem

]


(*=============SetCoordinateSysterm[ ][ ] for a special coordinate system ==============*)

SetCoordinateSystem[CoordinateSystem__][x_,y_,z_]:=
Module[{},

Unprotect[VectorQ,Jacobian];
(* Unprotect[Curl, Div, Grad]; *)

ClearAll[Christoffel,u,Metric,Jacobian];
ClearAll[VectorQ, ScalarQ, DotProduct, CrossProduct, Grad, Div, Curl, CovariantDerivative];

$CoordinateSystem=CoordinateSystem;

VectorNotation[];

VectorQ[a_]:= (Head[a]===Vector) ;



(* ============ Set up coordinates =================*)
u[1]=x;
u[2]=y;
u[3]=z;

(* ============ Metric Matrix for e^{i}, Metric[2] ===== *)

Metric[2][r_,t_,p_]=Metrictmp[CoordinateSystem][r,t,p];


(* ============ Metric Matrix for e_{i}, Metric[1] ============ *)
Metric[1][r_,t_,p_]=Inverse[Metric[2][r,t,p]];


(* ========== Jacobian==================== *)
Jacobian[r_,t_,p_]=Simplify[PowerExpand[Sqrt[Det[Metric[1][r,t,p]]]]];



(* ========== 1st and 2nd kind Christoffel symbol ============*)

Christoffel[1][j_,i_,k_,r_,t_,p_]:=
Module[{v},
 v[1]=r;
 v[2]=t;
 v[3]=p;
 1/2(D[Metric[1][r,t,p][[j,i]],v[k]]+
     D[Metric[1][r,t,p][[j,k]],v[i]]-
     D[Metric[1][r,t,p][[i,k]],v[j]])
]; 

Christoffel[2][i_,j_,k_,r_,t_,p_]:=
(
Sum[Metric[2][r,t,p][[i,n]] Christoffel[1][n,j,k,r,t,p],{n,1,3}]
);

Christoffel[1][r_,t_,p_]=
             Table[Simplify[Christoffel[1][i,j,k,r,t,p],Trig->False],{i,3},{j,3},{k,3}];
Christoffel[2][r_,t_,p_]=
             Table[Simplify[Christoffel[2][i,j,k,r,t,p],Trig->False],{i,3},{j,3},{k,3}];

Print["============The coordinate system is set up============="];
Print[" The Metric Matrix g^{i,j} is:"];
Print[ Metric[2][u[1],u[2],u[3]]  ];
Print[ "======================================================="];
Print[" The Metric Matrix g_{i,j} is:"];
Print[ Metric[1][u[1],u[2],u[3]] ];
Print[ "======================================================="];
Print[" The 1st kind Christoffel symbol matrix is:"];
Print[  Christoffel[1][u[1],u[2],u[3]]  ];
Print[ "======================================================="];
Print[" The 2st kind Christoffel symbol matrix is:"];
Print[ Christoffel[2][u[1],u[2],u[3]]  ];
Print[ "======================================================="];


(* ========== Levi-Civita symbol ====================*)

lc[i_,j_,k_]:=Signature[{i,j,k}];

(* ========= reverse between 1 and 2 =================*)
rev[i_ /; i==2]:=1 ;
rev[i_ /; i==1]:=2 ;


(* =========== Creat a Vector ============ *)
DefineVector[i_,f1_,f2_,f3_]:=
Module[{tmp1,tmp2},
If[i==1, tmp1={f1,f2,f3};
   tmp2=Metric[2][u[1],u[2],u[3]].tmp1,         
   tmp2={f1,f2,f3};
   tmp1=Metric[1][u[1],u[2],u[3]].tmp2  ];
   Vector[tmp1,tmp2] 
];

(* ========== Times a constant and Plus ===============*)

Vector/: Plus[Vector[a__],Vector[b__]]:=Vector@@Plus[List@@Vector[a]+List@@Vector[b]];

Vector/: Times[q_,Vector[a__]] := Vector@@Times[q, List@@Vector[a] ];

Vector/: D[Vector[a__],b_]:=Vector@@D[List@@Vector[a],b];

(* =========== DotProduct[a,b] ============== *)
DotProduct[a_,b_]:=Sum[a[[1,n]] b[[2,n]],{n,1,3}] ;

DotProduct2[a_,b_]:=Sum[a[[2,n]] b[[1,n]],{n,1,3}] ;


(* =========== CrossProduct[a,b] ======== *)

CrossProduct[a_,b_]:= 
Module[{tmp1,tmp2},
tmp1=Table[Sum[lc[i,j,k]a[[2,i]] b[[2,j]],{i,1,3},{j,1,3}]
                                     ,{k,3}] Jacobian[u[1],u[2],u[3]];
tmp2=Table[Sum[lc[i,j,k]a[[1,i]] b[[1,j]],{i,1,3},{j,1,3}]
                                     ,{k,3}]/Jacobian[u[1],u[2],u[3]];
 Vector[tmp1, tmp2]
] ;


(* =========== Grad[f]===========*)
Grad[f_]:=
Module[{tmp1,tmp2},
 tmp1=Table[D[f,u[i]],{i,3}];
 tmp2=Metric[2][u[1],u[2],u[3]].tmp1;
 Vector[tmp1,tmp2] 
] ;

(* =========== Div[a]============*)
Div[a_]:=
Module[{}, 
Sum[D[Jacobian[u[1],u[2],u[3]] a[[2,i]], u[i]],
   {i,3}]/Jacobian[u[1],u[2],u[3]]
] ;




(* ============ Curl[a] ==========*)
Curl[a_]:=
Module[{tmp1,tmp2},
tmp2=Table[Sum[lc[i,j,k] D[a[[1,j]], u[i]] ,{i,1,3},{j,1,3}]
                                           ,{k,3}]/Jacobian[u[1],u[2],u[3]];
 tmp1=Metric[1][u[1],u[2],u[3]].tmp2;
 Vector[tmp1,tmp2]
] ;

(* ============ CovariantDerivative[a,b]========*)
CovariantDerivative[a_,b_]:=
Module[{tmp1,tmp2},
If[VectorQ[b],
  (tmp2=Table[Sum[a[[2,i]] D[b[[2,j]],u[i]],{i,1,3}]+
   Sum[a[[2,i]] b[[2,k]] Christoffel[2][u[1],u[2],u[3]][[j,k,i]]
    ,{i,1,3},{k,1,3}],{j,3}];
    tmp1=Metric[1][u[1],u[2],u[3]].tmp2;
    Return[ Vector[tmp1,tmp2] ]
  ),
  Return [Sum[a[[2,i]] D[b,u[i]],{i,1,3}]]
  ]
] ;


(* =========== Laplacian[a] ==============*)
Laplacian[a_]:=
If[MatrixQ[a],
   Return[Grad[Div[a]]-Curl[Curl[a]]],
   Return[Div[Grad[a]]]
  ] ;

(* =========== AbsoluteValue[a] ===========*)
AbsoluteValue[a_]:=
PowerExpand[Sqrt[DotProduct[a,a]]] ;

(* =========== UnitVector[a] ==============*)
UnitVector[a_]:=
a/PowerExpand[Sqrt[DotProduct[a,a]]] ;

(* =========== Parallel[a,b] ============*)
Parallel[a_,b_]:=
DotProduct[a,UnitVector[b]] ;

ParallelUnitVector[a_,b_]:=
DotProduct[a,b] ;


(* =========== Perp[a,b] ================*)
Perp[a_,b_]:=
Module[{tmp1,tmp2},
tmp1=UnitVector[b];
 Simplify[CrossProduct[CrossProduct[tmp1,a],tmp1]]
] ;

PerpUnitVector[a_,b_]:=
CrossProduct[CrossProduct[b,a],b] ;
 
(* ========== Subscript[a,i]  and CovariantComponent[a,i] =================*)

Subscript[a_,i_]:=CovariantComponent[a,i] /; (VectorQ[a] && IntegerQ[i]) ;

CovariantComponent[a_,i_]:=a[[1,i]] /; (VectorQ[a] && IntegerQ[i]);

(* ========== Superscript[a,i] and ContravariantComponent[a,i]================*)

Superscript[a_,i_]:=ContravariantComponent[a,i] /; (VectorQ[a] && IntegerQ[i]) ;

Superscript[a_,i_]:=Power[a,i] /; !(VectorQ[a] && IntegerQ[i]);

ContravariantComponent[a_,i_]:=a[[2,i]] /; (VectorQ[a]==True && IntegerQ[i]==True);



Print[{CoordinateSystem}, " is set up."];

]




(* =================DefineCoordinateSystem ========================*)
DefineCoordinateSystem[CoordinateSystem__][c1_,c2_,c3_][matrix__]:=
Module[{r,t,p},
Metrictmp[CoordinateSystem][r_,t_,p_]=(matrix /. {c1->r,c2->t,c3->p});
Print[{CoordinateSystem}," is added to current Mathematica session."];
Print["The Riemann Metric Tensor is:"];
Print["============================="];
Print[Metrictmp[CoordinateSystem][r,t,p]];
Print["Use SetCoordinateSystem to set up a coordinate system"];
]


(* ================= Coordinate Systems ========================*)

(* ----------------Cylindrical ------------------------ *)
Metrictmp["Cylindrical"][r_,t_,z_]:={{1,0,0},{0,1/r^2,0},{0,0,1}}
 

(* ----------------Cartesian---------------------------*)
Metrictmp["Cartesian"][x_,y_,z_]:={{1,0,0},{0,1,0},{0,0,1}}


(* ----------------Spherical --------------------------*)
Metrictmp["Spherical"][r_,t_,p_]:={{1,0,0},{0,1/r^2,0},
                               {0,0,1/(r^2 Sin[t]^2)}}

(* --------------STKCircular-----------------------------*)
(* --------------No small parameter here     ------------*)
Metrictmp["STKCircular",R0_,ep_][r_,t_,p_]:=
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
Metrictmp["TKGeneralPerp",g1_, g2_,g3_][r_,t_,p_]:=
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

Metrictmp["TKCircularFlux",R0_,ept_][r_,t_,p_]:=
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

Metrictmp["TKCircularFlux2",R0_,ept_][r_,t_,p_]:=
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
Metrictmp["TKCircular",R0_,ep_][r_,t_,p_]:=
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

Metrictmp["TKCircular6",R0_,ep_][r_,t_,p_]:=
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

Metrictmp["TKCircular5",R0_,ep_][r_,t_,p_]:=
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

Metrictmp["TKCircular4",R0_,ep_][r_,t_,p_]:=
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

Metrictmp["TKShiftedCircular4",R0_,ep_,R0dp_][r_,t_,p_]:=
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

Metrictmp["TKCircular3",R0_,ep_][r_,t_,p_]:=
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

Metrictmp["TKCircular2",R0_,ep_][r_,t_,p_]:=
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



Print["The Vector Calculus on General Coordinate is loaded in"]
Print["Use SetCoordinateSystem[ ]  to set up a coordinate system"]
Print["The default CoordinateSystem is  None "]
$CoordinateSystem="None"
SetCoordinateSystem["None"]


End[]

EndPackage[]



-- End --
