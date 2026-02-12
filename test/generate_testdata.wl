(* generate_testdata.wl - Generate test/testdata_expected.hh *)
(* Run with: /home/ker/.local/bin/math -script test/generate_testdata.wl *)

outputFile = "test/testdata_expected.hh";

(* Helper: convert Mathematica polynomial to CLPoly C++ construction string *)
(* Uses polynomial_ZZ constructed via variable operations *)
polyToCSrc[poly_, vars_] := Module[
  {terms, result, termStrs},
  If[poly === 0, Return["polynomial_ZZ()"]];
  terms = MonomialList[poly, vars, "DegreeLexicographic"];
  termStrs = Map[
    Module[{coeff, monStr, varParts},
      coeff = # /. Thread[vars -> 1];
      varParts = {};
      Do[
        Module[{exp = Exponent[#, v]},
          If[exp > 0,
            If[exp == 1,
              AppendTo[varParts, ToString[v]],
              AppendTo[varParts, "pow(" <> ToString[v] <> "," <> ToString[exp] <> ")"]
            ]
          ]
        ],
        {v, vars}
      ];
      monStr = If[varParts === {}, "", StringRiffle[varParts, "*"]];
      Which[
        monStr === "" && coeff === 1, "1",
        monStr === "" && coeff === -1, "-1",
        monStr === "", "polynomial_ZZ({{{}," <> ToString[coeff] <> "}})",
        coeff === 1, monStr,
        coeff === -1, "-" <> monStr,
        True, ToString[coeff] <> "*" <> monStr
      ]
    ] &,
    terms
  ];
  (* wrap integer constants properly *)
  result = StringRiffle[termStrs, "+"];
  (* fix +-  to - *)
  result = StringReplace[result, "+-" -> "-"];
  result
];

(* Helper: convert integer to ZZ string *)
zzStr[n_] := "ZZ(\"" <> ToString[n] <> "\")";

(* Open file *)
stream = OpenWrite[outputFile];

WriteString[stream, "/**\n"];
WriteString[stream, " * @file testdata_expected.hh\n"];
WriteString[stream, " * @brief Auto-generated test data from Mathematica\n"];
WriteString[stream, " * DO NOT EDIT - regenerate with: math -script test/generate_testdata.wl\n"];
WriteString[stream, " */\n"];
WriteString[stream, "#ifndef CLPOLY_TESTDATA_EXPECTED_HH\n"];
WriteString[stream, "#define CLPOLY_TESTDATA_EXPECTED_HH\n\n"];
WriteString[stream, "#include <clpoly/clpoly.hh>\n"];
WriteString[stream, "#include <vector>\n"];
WriteString[stream, "#include <utility>\n\n"];
WriteString[stream, "namespace testdata {\n\n"];
WriteString[stream, "using namespace clpoly;\n\n"];

(* ============ ARITHMETIC TEST DATA ============ *)
WriteString[stream, "// ==================== Arithmetic ====================\n\n"];

(* Simple: single variable *)
SeedRandom[42];
vars1 = {x};
f1s = 3x^3 - 2x^2 + x - 5;
g1s = x^2 + 2x + 1;

WriteString[stream, "inline polynomial_ZZ arith_f1_simple() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[f1s, vars1] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ arith_g1_simple() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[g1s, vars1] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ arith_add_simple() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[Expand[f1s + g1s], vars1] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ arith_sub_simple() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[Expand[f1s - g1s], vars1] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ arith_mul_simple() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[Expand[f1s * g1s], vars1] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ arith_pow2_simple() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[Expand[f1s^2], vars1] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ arith_pow3_simple() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[Expand[f1s^3], vars1] <> ";\n"];
WriteString[stream, "}\n\n"];

(* Medium: two variables *)
f1m = 2x^3*y^2 - 3x*y^4 + x^2*y - 7;
g1m = x^2*y - y^3 + 2x + 1;

WriteString[stream, "inline polynomial_ZZ arith_f1_medium() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[f1m, {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ arith_g1_medium() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[g1m, {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ arith_mul_medium() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[Expand[f1m * g1m], {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];

(* Complex: three variables *)
f1c = x^3*y*z - 2*x^2*z^3 + 3*y^2*z^2 - x*y + 5;
g1c = x*y^2*z + z^3 - x^2 + y - 1;

WriteString[stream, "inline polynomial_ZZ arith_f1_complex() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\"), z(\"z\");\n"];
WriteString[stream, "    return " <> polyToCSrc[f1c, {x, y, z}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ arith_g1_complex() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\"), z(\"z\");\n"];
WriteString[stream, "    return " <> polyToCSrc[g1c, {x, y, z}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ arith_mul_complex() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\"), z(\"z\");\n"];
WriteString[stream, "    return " <> polyToCSrc[Expand[f1c * g1c], {x, y, z}] <> ";\n"];
WriteString[stream, "}\n\n"];

(* ============ PREM/PQUO TEST DATA ============ *)
WriteString[stream, "// ==================== Prem/Pquo ====================\n\n"];

(* Simple prem *)
fp = 2x^4 + 3x^3 - x + 7;
gp = x^2 + x - 1;
premResult = PolynomialRemainder[fp * LeadingCoefficient[gp, x]^(Exponent[fp, x] - Exponent[gp, x] + 1), gp, x];
(* Actually prem in Mathematica is PolynomialPseudoRemainder *)

(* Use proper pseudo-remainder computation *)
prem1f = x^3 + 2x^2 - x + 3;
prem1g = x^2 - x + 1;
lcg = Coefficient[prem1g, x, 2];  (* = 1 *)
degDiff = Exponent[prem1f, x] - Exponent[prem1g, x] + 1;  (* = 2 *)
prem1Result = PolynomialRemainder[lcg^degDiff * prem1f, prem1g, x];
pquo1Result = PolynomialQuotient[lcg^degDiff * prem1f, prem1g, x];

WriteString[stream, "inline polynomial_ZZ prem_f1() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[prem1f, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ prem_g1() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[prem1g, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ prem_result1() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[prem1Result, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ pquo_result1() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[pquo1Result, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

(* Medium: multivariate prem *)
prem2f = x^3 + y*x^2 + (y^2-1)*x + y^3;
prem2g = x^2 + y*x + y^2;
lc2 = Coefficient[prem2g, x, 2];
dd2 = Exponent[prem2f, x] - Exponent[prem2g, x] + 1;
prem2Result = PolynomialRemainder[lc2^dd2 * prem2f, prem2g, x];
pquo2Result = PolynomialQuotient[lc2^dd2 * prem2f, prem2g, x];

WriteString[stream, "inline polynomial_ZZ prem_f2() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[prem2f, {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ prem_g2() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[prem2g, {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ prem_result2() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[prem2Result, {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ pquo_result2() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[pquo2Result, {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];

(* ============ GCD TEST DATA ============ *)
WriteString[stream, "// ==================== GCD ====================\n\n"];

(* Simple: known GCD *)
gcdCommon1 = x^2 + 1;
gcdF1 = Expand[gcdCommon1 * (x^2 - 2x + 3)];
gcdG1 = Expand[gcdCommon1 * (x + 5)];
gcdResult1 = PolynomialGCD[gcdF1, gcdG1];

WriteString[stream, "inline polynomial_ZZ gcd_f1() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[gcdF1, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ gcd_g1() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[gcdG1, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ gcd_result1() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[gcdResult1, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

(* Coprime polynomials *)
gcdF2 = x^3 + x + 1;
gcdG2 = x^2 + x + 1;

WriteString[stream, "inline polynomial_ZZ gcd_f2_coprime() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[gcdF2, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ gcd_g2_coprime() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[gcdG2, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

(* Multivariate GCD *)
gcdCommon3 = x*y + 1;
gcdF3 = Expand[gcdCommon3 * (x^2 - y)];
gcdG3 = Expand[gcdCommon3 * (y^2 + x)];
gcdResult3 = PolynomialGCD[gcdF3, gcdG3];

WriteString[stream, "inline polynomial_ZZ gcd_f3_mv() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[gcdF3, {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ gcd_g3_mv() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[gcdG3, {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ gcd_result3_mv() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[gcdResult3, {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];

(* ============ RESULTANT / DISCRIMINANT TEST DATA ============ *)
WriteString[stream, "// ==================== Resultant / Discriminant ====================\n\n"];

(* Simple resultant *)
resF1 = x^2 + x + 1;
resG1 = x^3 - 1;
resResult1 = Resultant[resF1, resG1, x];

WriteString[stream, "inline polynomial_ZZ res_f1() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[resF1, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ res_g1() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[resG1, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline ZZ res_result1_val() { return " <> zzStr[resResult1] <> "; }\n\n"];

(* Multivariate resultant *)
resF2 = x^2 + y*x + y^2;
resG2 = x^3 - y^3;
resResult2 = Resultant[resF2, resG2, x];

WriteString[stream, "inline polynomial_ZZ res_f2() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[resF2, {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ res_g2() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[resG2, {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ res_result2() {\n"];
WriteString[stream, "    variable y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[Expand[resResult2], {y}] <> ";\n"];
WriteString[stream, "}\n\n"];

(* Resultant with common factor (should be 0) *)
resCommon = x + y;
resF3 = Expand[resCommon * (x - y)];
resG3 = Expand[resCommon * (x + 1)];
resResult3 = Resultant[resF3, resG3, x];

WriteString[stream, "inline polynomial_ZZ res_f3_common() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[resF3, {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ res_g3_common() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[resG3, {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];

(* Discriminant *)
discF1 = x^3 - 3x + 2;  (* has repeated root at x=1 *)
discResult1 = Discriminant[discF1, x];

WriteString[stream, "inline polynomial_ZZ disc_f1_repeated() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[discF1, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline ZZ disc_result1_val() { return " <> zzStr[discResult1] <> "; }\n\n"];

discF2 = x^3 + x + 1;
discResult2 = Discriminant[discF2, x];

WriteString[stream, "inline polynomial_ZZ disc_f2() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[discF2, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline ZZ disc_result2_val() { return " <> zzStr[discResult2] <> "; }\n\n"];

(* Multivariate discriminant *)
discF3 = x^2 + y*x + y^2 - 1;
discResult3 = Discriminant[discF3, x];

WriteString[stream, "inline polynomial_ZZ disc_f3_mv() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[discF3, {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ disc_result3_mv() {\n"];
WriteString[stream, "    variable y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[Expand[discResult3], {y}] <> ";\n"];
WriteString[stream, "}\n\n"];

(* ============ SQUAREFREE TEST DATA ============ *)
WriteString[stream, "// ==================== Squarefree ====================\n\n"];

(* Polynomial with square factors *)
sqfF1 = Expand[(x+1)^2 * (x-2)];

WriteString[stream, "inline polynomial_ZZ sqf_f1() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[sqfF1, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

sqfF2 = Expand[(x+1)^2 * (x-1)^3];

WriteString[stream, "inline polynomial_ZZ sqf_f2() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[sqfF2, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

(* Already squarefree *)
sqfF3 = x^3 + x + 1;

WriteString[stream, "inline polynomial_ZZ sqf_f3_already() {\n"];
WriteString[stream, "    variable x(\"x\");\n"];
WriteString[stream, "    return " <> polyToCSrc[sqfF3, {x}] <> ";\n"];
WriteString[stream, "}\n\n"];

(* Multivariate square factor *)
sqfF4 = Expand[(x + y)^2 * (x - y + 1)];

WriteString[stream, "inline polynomial_ZZ sqf_f4_mv() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\");\n"];
WriteString[stream, "    return " <> polyToCSrc[sqfF4, {x, y}] <> ";\n"];
WriteString[stream, "}\n\n"];


(* ============ COMPLEX / RANDOM TEST DATA ============ *)
WriteString[stream, "// ==================== Complex (random with fixed seed) ====================\n\n"];

SeedRandom[12345];

(* Generate complex 3-variable polynomials *)
randomPoly3[deg_, nterms_] := Module[{terms = {}, i, e1, e2, e3, c},
  While[Length[terms] < nterms,
    e1 = RandomInteger[{0, deg}];
    e2 = RandomInteger[{0, deg - e1}];
    e3 = RandomInteger[{0, deg - e1 - e2}];
    c = RandomInteger[{-10, 10}];
    If[c != 0, AppendTo[terms, c * x^e1 * y^e2 * z^e3]];
  ];
  Expand[Total[terms]]
];

rc1 = randomPoly3[5, 6];
rc2 = randomPoly3[5, 6];

WriteString[stream, "inline polynomial_ZZ complex_f1() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\"), z(\"z\");\n"];
WriteString[stream, "    return " <> polyToCSrc[rc1, {x, y, z}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ complex_g1() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\"), z(\"z\");\n"];
WriteString[stream, "    return " <> polyToCSrc[rc2, {x, y, z}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ complex_mul1() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\"), z(\"z\");\n"];
WriteString[stream, "    return " <> polyToCSrc[Expand[rc1 * rc2], {x, y, z}] <> ";\n"];
WriteString[stream, "}\n\n"];

(* GCD with known factor *)
rcCommon = x*y + z + 1;
rcGcdF = Expand[rcCommon * rc1];
rcGcdG = Expand[rcCommon * rc2];

WriteString[stream, "inline polynomial_ZZ complex_gcd_f() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\"), z(\"z\");\n"];
WriteString[stream, "    return " <> polyToCSrc[rcGcdF, {x, y, z}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ complex_gcd_g() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\"), z(\"z\");\n"];
WriteString[stream, "    return " <> polyToCSrc[rcGcdG, {x, y, z}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "inline polynomial_ZZ complex_gcd_common() {\n"];
WriteString[stream, "    variable x(\"x\"), y(\"y\"), z(\"z\");\n"];
WriteString[stream, "    return " <> polyToCSrc[rcCommon, {x, y, z}] <> ";\n"];
WriteString[stream, "}\n\n"];

(* Resultant *)
rcResResult = Resultant[rc1, rc2, x];

WriteString[stream, "inline polynomial_ZZ complex_res_result() {\n"];
WriteString[stream, "    variable y(\"y\"), z(\"z\");\n"];
WriteString[stream, "    return " <> polyToCSrc[Expand[rcResResult], {y, z}] <> ";\n"];
WriteString[stream, "}\n\n"];

WriteString[stream, "} // namespace testdata\n\n"];
WriteString[stream, "#endif\n"];

Close[stream];

Print["Generated " <> outputFile <> " successfully."];
