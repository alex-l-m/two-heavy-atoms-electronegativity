#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory[Directory[]]
symbols=ReadList[File["selected_elements.txt"],String]
elements=Map[Entity["Element",#1]&,symbols]
atomicNumbers=Map[#1[EntityProperty["Element","AtomicNumber"]]&,elements]
outTable=Transpose[{Prepend[symbols,"symbol"],Prepend[atomicNumbers,"atomic_number"]}]
Export["atomic_numbers.csv",outTable]
