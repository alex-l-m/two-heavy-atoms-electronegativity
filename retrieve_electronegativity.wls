#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory[Directory[]]
targetElementTable = Import["target_elements.csv"]
columnIndex = Position[targetElementTable[[1]], "symbol"][[1]][[1]]
targetElementSymbols = Drop[Transpose[targetElementTable][[columnIndex]], 1]
elements = Map[Entity["Element", #] &, targetElementSymbols]
electronegativities = Map[#[EntityProperty["Element", "Electronegativity"]] &, elements]
symbols = Map[#[EntityProperty["Element", "AtomicSymbol"]] &, elements]
table = Transpose[{Prepend[symbols, "symbol"], Prepend[electronegativities, "pauling_electronegativity"]}]
Export["pauling_electronegativity.csv", table]
