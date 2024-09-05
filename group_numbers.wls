#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory[Directory[]]
symbols=ReadList[File["selected_elements.txt"],String]
elements=Map[Entity["Element",#1]&,symbols]
groupNumbers=Map[#1[EntityProperty["Element","Group"]]&,elements]
outTable=Transpose[{Prepend[symbols,"symbol"],Prepend[groupNumbers,"group"]}]
Export["group_numbers.csv",outTable]
