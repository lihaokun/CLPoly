"""
Pass 9: refine_gen — 从 MIRFunc + REFINEMENT_MAP 生成精化定理骨架

输入：MIRFunc 列表（Pass 8 已生成 L1 corpus）
输出：
  * proof/lean/CLPoly/Refinement/GeneratedSkeletons/*.lean（待证明骨架）
  * proof/lean/CLPoly/Refinement/Generated.lean（已证明公开契约的集中入口）

生成器绝不覆盖手写且已证明的 Refinement/*.lean。已证明公开契约保留原始
C++ 函数名（包括前导双下划线），并桥接到严格执行语义的内部证明。

四类定理形态（由 result_kind 决定）：

  direct_eq    : l1_ir x = l2 x                        （纯数值函数）
  map_eq       : toPolyList (l1_ir f) p = l2 (f.toPoly p)  （列表输出，Array→List）
  map_eq_pair  : toPoly (l1_ir f).snd = l2 (f.toPoly p)    （pair 输出，取第二个分量）
  pair_eq      : toPoly (l1_ir ...).1 = (l2 ...).1 ∧ ...  （多返回值）
  conditional  : l1_ir f = .ok result → toPoly result = ...（有 throw/Except）

参数信息的三个来源（优先级从高到低）：
  1. MIRFunc（来自 Pass 1-8 流水线）
  2. REFINEMENT_MAP 的 l1_params / l2_call 字段
  3. _FIRST_PARAM / _HAS_SPARSE_POLY 固定表
"""

from __future__ import annotations
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from class_map import REFINEMENT_MAP
from ir_types import (
    BaseType,
    MIRFunc,
    NamedType,
    ArrayType,
    PairType,
)


# ============================================================
# 输出目录
# ============================================================

REFINEMENT_DIR = Path(__file__).resolve().parent.parent.parent / "lean" / "CLPoly" / "Refinement"


# ============================================================
# 片段生成
# ============================================================

LEAN_KEYWORDS = frozenset({
    "let", "fun", "match", "with", "do", "if", "then", "else",
    "where", "by", "and", "or", "not", "in", "of", "from",
    "def", "theorem", "sorry",
})


def _safe_name(s: str) -> str:
    if s in LEAN_KEYWORDS:
        return f"«{s}»"
    return s


def _to_lean_type(ty) -> str:
    if isinstance(ty, BaseType):
        return ty.value
    if isinstance(ty, NamedType):
        return ty.name
    if isinstance(ty, ArrayType):
        return f"Array ({_to_lean_type(ty.elem)})"
    if isinstance(ty, PairType):
        return f"({_to_lean_type(ty.fst)} × {_to_lean_type(ty.snd)})"
    return " sorry /* type */ "


# 固定参数表（当 MIRFunc 和 REFINEMENT_MAP.l1_params 都不可用时使用）
_PARAM_TABLE = {
    "__make_zp":        [("val", "Int64"), ("modulus", "UInt64")],
    "__upoly_make_monic": [("f", "SparsePolyZp")],
    "__upoly_divmod":   [("q", "SparsePolyZp"), ("r", "SparsePolyZp"),
                         ("f", "SparsePolyZp"), ("g", "SparsePolyZp")],
    "__squarefree_Zp": [("f", "SparsePolyZp")],
    "__ddf_Zp": [("f", "SparsePolyZp")],
    "__symmetric_mod": [("a", "ZZ"), ("m", "ZZ")],
    "__binomial": [("n", "Int64"), ("k", "Int64")],
    "__isqrt_ceil": [("n", "ZZ")],
    "__upoly_mod":      [("f", "SparsePolyZp"), ("g", "SparsePolyZp")],
    "__edf_Zp":         [("result", "Array SparsePolyZp"), ("f", "SparsePolyZp"),
                         ("d", "UInt64"), ("rng", "Rng")],
    "__factor_Zp":      [("f", "SparsePolyZp")],
    "factorize":        [("F", "SparsePolyZZ")],
}

# 哪些 base_name 的参数包含 SparsePolyZp
_SPARSE_POLY_FUNCS = {
    "factorize", "__squarefree_Zp", "__ddf_Zp", "__edf_Zp",
    "__factor_Zp", "__upoly_make_monic", "__upoly_divmod",
}


def _get_params(base_name: str) -> list[tuple[str, str]]:
    return _PARAM_TABLE.get(base_name, [("x", "ZZ")])


def _has_sparse_poly(base_name: str) -> bool:
    return base_name in _SPARSE_POLY_FUNCS


def _first_sparse_poly_param(base_name: str) -> str:
    for name, ty in _PARAM_TABLE.get(base_name, []):
        if "SparsePolyZp" in ty:
            return name
    return "x"


def _emit_theorem_body(l1_name: str, info: dict, f: Optional[MIRFunc] = None,
                       base_name: str = "") -> str:
    """生成定理结论（= 号右侧的表达式）。"""
    kind = info["result_kind"]
    l2_name = info["l2_name"]
    l2_call = info.get("l2_call", l2_name)
    bn = base_name or l1_name.replace("_ir", "")

    # L1 调用参数（Lean 空格分隔函数参数）
    if f is not None:
        l1_call_params = " ".join(_safe_name(p.name) for p in f.params)
    else:
        l1_call_params = " ".join(_safe_name(p[0]) for p in _get_params(bn))
    l1_call = f"Generated.{l1_name} {l1_call_params}"

    if kind == "direct_eq":
        if _has_sparse_poly(bn):
            return f"  SparsePolyZp.toPoly p ({l1_call}) = {l2_call}"
        else:
            return f"  {l1_call} = {l2_call}"
    elif kind == "map_eq":
        return f"  toPolyList ({l1_call}) p = {l2_call}"
    elif kind == "map_eq_pair":
        return f"  SparsePolyZp.toPoly p ({l1_call}).snd = SparsePolyZp.toPoly p ({l2_call})"
    elif kind == "pair_eq":
        return (
            f"  SparsePolyZp.toPoly p ({l1_call}).1 = SparsePolyZp.toPoly p ({l2_call}).1 ∧\n"
            f"  SparsePolyZp.toPoly p ({l1_call}).2 = SparsePolyZp.toPoly p ({l2_call}).2"
        )
    elif kind == "conditional":
        return (
            f"  {l1_call} = .ok result →\n"
            f"  SparsePolyZp.toPoly p result = {l2_call}"
        )
    else:
        return f"  sorry  /* unknown result_kind: {kind} */"


def _emit_refinement_theorem(f: Optional[MIRFunc], info: dict, base_name: str = "") -> str:
    if f is not None:
        l1_lean_name = f.lean_name
        bn = f.base_name
    else:
        l1_lean_name = info.get("l1_name", f"{base_name}_ir")
        bn = base_name

    cpp_source = info.get("cpp_source", "unknown")
    doc = info.get("doc", "")
    l2_name = info["l2_name"]

    lines: list[str] = []
    # docstring
    lines.append("/--")
    lines.append(f"  L1 `{l1_lean_name}` (C++: `{cpp_source}`) → L2 `{l2_name}`")
    if doc:
        lines.append(f"")
        lines.append(f"  {doc}")
    lines.append(f"")
    lines.append(f"  C++ source: `{cpp_source}` — Clang AST → cpp2lean v2 Pass 1-8 → Corpus.lean")
    lines.append(f"  L2 model : Algorithm/{info.get('refinement_file', '?')}.lean — hand-written, proven correct")
    lines.append(f"  Bridge   : SparsePolyZp.toPoly (see Math/Univariate.lean)")
    lines.append(f"")
    lines.append(f"  Proof status: Skeleton (sorry) — fill to complete L1→L2 verification chain")
    lines.append("-/")

    # 参数列表
    if f is not None:
        params = [(p.name, _to_lean_type(p.ty)) for p in f.params]
    else:
        params = info.get("l1_params", _get_params(bn))

    # 定理签名
    lines.append(f"theorem {l1_lean_name}_refines (p : ℕ) [hp : Fact (Nat.Prime p)]")
    for pname, ptype in params:
        lines.append(f"    ({_safe_name(pname)} : {ptype})")

    # 前提条件（SparsePolyZp 参数）
    has_spzp = any("SparsePolyZp" in t for _, t in params)
    for pname, ptype in params:
        if "SparsePolyZp" in ptype:
            lines.append(f"    (hwf_{pname} : SparsePolyZp.WellFormed p {_safe_name(pname)})")
            lines.append(f"    (hred_{pname} : SparsePolyZp.AllReduced p {_safe_name(pname)}.toList)")
    if has_spzp:
        lines.append(f"    (hp_size : 2 * p ≤ UInt64.size)")

    # 结论
    body = _emit_theorem_body(l1_lean_name, info, f, base_name=bn)
    lines.append(f"    : {body} :=")
    lines.append(f"  by")
    lines.append(f"  sorry")
    lines.append("")

    return "\n".join(lines)


def _emit_refinement_file(funcs: list[Optional[MIRFunc]], base_names: list[str], module: str) -> str:
    lines: list[str] = []
    lines.append("-- Auto-generated by cpp2lean v2 Pass 9 (refine_gen)")
    lines.append("-- L1 → L2 refinement theorem skeletons")
    lines.append("")

    seen_imports: set[str] = set()
    for bn in base_names:
        imp = REFINEMENT_MAP.get(bn, {}).get("l2_import", "")
        if imp and imp not in seen_imports:
            lines.append(f"import {imp}")
            seen_imports.add(imp)
    lines.append("import CLPoly.Model")
    lines.append("import CLPoly.Generated.Corpus")
    lines.append("import CLPoly.Refinement.Basic")
    lines.append("import CLPoly.Math.Univariate")
    lines.append("")
    lines.append("set_option autoImplicit false")
    lines.append("")
    lines.append("open Polynomial")
    lines.append("open CLPoly.Math")
    lines.append("")
    lines.append("namespace Refinement")
    lines.append("")
    lines.append("variable {p : ℕ} [hp : Fact (Nat.Prime p)]")
    lines.append("")

    for i, bn in enumerate(base_names):
        f = funcs[i] if i < len(funcs) else None
        info = REFINEMENT_MAP.get(bn, {})
        theorem = _emit_refinement_theorem(f, info, base_name=bn)
        lines.append(theorem)

    lines.append("end Refinement")
    lines.append("")
    return "\n".join(lines)


def _emit_verified_contract(info: dict) -> str:
    """Emit one checked public wrapper around an already-proved strict theorem."""
    contract = info["verified_contract"]
    theorem_name = contract["theorem_name"]
    proof_theorem = contract["proof_theorem"]
    if contract["kind"] == "strict_select_prime":
        return f'''/-- Generated public contract for the original C++
`__select_prime` entry.  The strict L1 program uses well-founded machine-prime
enumeration and executes modular reduction, derivative/GCD, make-monic, DDF,
and EDF for every accepted candidate. -/
theorem {theorem_name}
    {{State : Type}} (engine : Generated.StrictEDF.RandomEngine State)
    (provider : StrictSelectPrime.CandidateRuntimeProvider engine)
    (initialRng : State) (useLargePrime : Bool) (f : SparsePolyZZ)
    (hinitialPrimeCorrect : Nat.Prime
      (if useLargePrime then
        ((18446744073709551615 : UInt64) - 58).toNat
      else (2 : UInt64).toNat))
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (hdegree : 2 ≤ (SparsePolyZZ.toPoly f).natDegree)
    (hdegreeBound : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 62)
    (hlcSemantic : ∀ p : UInt64, Nat.Prime p.toNat →
      ((SparsePolyZZ.front! f).2 : ZMod p.toNat) =
        ((SparsePolyZZ.toPoly f).leadingCoeff : ZMod p.toNat))
    (result : PrimeSelectionResult)
    (hrun : Generated.StrictSelectPrime.__select_prime_raw_ir
      (StrictSelectPrime.selectPrimeRawOps
        (StrictSelectPrime.concreteTryCandidate engine provider))
      (StrictSelectPrime.selectPrimeTermination
        (StrictSelectPrime.concreteTryCandidate engine provider))
      initialRng useLargePrime f = .ok result) :
    StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) result ∧
      StrictSelectPrime.SelectionPhysical result := by
  exact {proof_theorem} engine provider initialRng useLargePrime f
    hinitialPrimeCorrect hcanonical hnonempty hdegree hdegreeBound
      hlcSemantic result hrun
'''
    if contract["kind"] == "strict_factor_zp":
        return f"""/-- Generated public end-to-end contract for the original C++ `__factor_Zp`
entry.  Its executable side is the source-shaped strict L1 entry, including
make-monic, SQF, DDF, EDF, multiplicity attachment, RNG threading and the
`std::sort` permutation boundary. -/
theorem {theorem_name}
    {{State : Type}} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (sqfPhysical : StrictSquarefreeZp.YunRawGCDWorkspaceProvider this hcfg)
    (providers : StrictDDF.DDFRawProviders this)
    (termination : Generated.StrictEDF.EDFTermination
      (StrictEDF.strictEDFRawOps engine this providers))
    (sort : StrictFactorZp.SortByDegreeProvider) (initialRng : State)
    (f : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat f)
    (hnonempty : 0 < f.size)
    (hdegreePositive :
      0 < (SparsePolyZp.toPoly this._p.toNat f).natDegree)
    (hdegreeBound :
      (SparsePolyZp.toPoly this._p.toNat f).natDegree < 2 ^ 62) :
    let ops : Generated.StrictFactorZp.FactorZpRawOps State := {{
      makeMonic := StrictSquarefreeZp.upolyMakeMonicIR this
      squarefree := StrictFactorZp.strictSQFCall this hcfg sqfPhysical
      ddf := StrictFactorZp.strictDDFCall this providers
      edf := StrictFactorZp.strictEDFCall engine this providers termination
      sortByDegree := sort.run }}
    ∃ lc output,
      Generated.StrictFactorZp.__factor_Zp_raw_ir ops initialRng f =
        .ok (lc, output) ∧
      FactorZpCorrect (SparsePolyZp.toPoly this._p.toNat f)
        (Zp.toZMod this._p.toNat lc)
        (StrictFactorZp.factorResultToL2 this._p.toNat output) ∧
      (∀ item ∈ output.toList,
        SparsePolyZp.Canonical this._p.toNat item.1) := by
  exact {proof_theorem} engine this hcfg sqfPhysical providers termination
    sort initialRng f hcanonical hnonempty hdegreePositive hdegreeBound
"""
    if contract["kind"] == "strict_factor_zz":
        return f"""/-- Generated public end-to-end contract for the original
C++ `__factor_squarefree_primitive_ZZ` entry.  Its executable side fixes the
actual generated prime-selection, modular factorization, quadratic Hensel
lifting, van-Hoeij/Zassenhaus recombination, and source branch structure. -/
theorem {theorem_name}
    {{State : Type}} (engine : Generated.StrictEDF.RandomEngine State)
    (provider : StrictSelectPrime.CandidateRuntimeProvider engine)
    (initialRng : State) (useLargePrime : Bool) (f : SparsePolyZZ)
    (hinitialPrimeCorrect : Nat.Prime
      (if useLargePrime then
        ((18446744073709551615 : UInt64) - 58).toNat
      else (2 : UInt64).toNat))
    (hhenselInvariant : ∀ (selection : PrimeSelectionResult)
      (hp : Nat.Prime selection.prime.toNat) (aTarget : Int32),
      let candidate := provider.physical selection.prime hp
      @StrictFactorZZ.HenselLiftRuntimeReadiness candidate.dense
        StrictHensel.concreteDivmodTermination candidate.providers.mul
        f selection.factors aTarget)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hnonempty : 0 < f.size)
    (hdegree : 2 ≤ (SparsePolyZZ.toPoly f).natDegree)
    (hdegree62 : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 62)
    (hdegree63 : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (output : Array SparsePolyZZ)
    (hrun : Generated.StrictFactorZZ.__factor_squarefree_primitive_ZZ_raw_ir
      (StrictFactorZZ.concreteRecombineFactorZZRawOps
        (StrictFactorZZ.concreteSelectPrime engine provider initialRng)
        (StrictFactorZZ.concreteHenselLift engine provider))
      useLargePrime f = .ok output) :
    FactorZZCorrect (SparsePolyZZ.toPoly f)
      (output.toList.map SparsePolyZZ.toPoly) := by
  exact {proof_theorem} engine provider initialRng useLargePrime f
    hinitialPrimeCorrect hhenselInvariant hcanonical hprimitive hnonempty
      hdegree hdegree62 hdegree63 leading hleading output hrun
"""
    if contract["kind"] == "strict_squarefree_zp":
        return f"""/-- Generated public contract for the original C++ `__squarefree_Zp` entry.
The executable side is the strict, well-founded L1 semantics; the result is
identified with the L2 `sqfZp` model. -/
theorem {theorem_name}
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : StrictSquarefreeZp.YunRawGCDWorkspaceProvider this hcfg)
    (source : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat source)
    (hmonic : (SparsePolyZp.toPoly this._p.toNat source).Monic)
    (hnonempty : 0 < source.size)
    (hpositive : 0 < (SparsePolyZp.toPoly this._p.toNat source).natDegree)
    (hbound : CLPoly.Impl.StrictPolynomialGCDRefinement.sparseDenseLength source ≤ 2 ^ 63) :
    ∃ factors,
      Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir
          (StrictSquarefreeGenerated.strictSQFRawOps this hcfg physical)
          source (fun _ => ⟨hcanonical, hmonic, hnonempty, hpositive,
            hbound⟩) = .ok factors ∧
      toPolyList factors this._p.toNat =
        sqfZp (SparsePolyZp.toPoly this._p.toNat source) ∧
      ∀ item ∈ factors.toList,
        SparsePolyZp.Canonical this._p.toNat item.1 := by
  exact {proof_theorem} this hcfg physical source hcanonical hmonic
    hnonempty hpositive hbound
"""
    if contract["kind"] == "strict_ddf_zp":
        return f"""/-- Generated public contract for the original C++ `__ddf_Zp` entry.
The executable side is the strict, well-founded L1 semantics; the result is
identified with the L2 `ddf` model. -/
theorem {theorem_name}
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this) (f : SparsePolyZp)
    (hfPrime : f[0]!.2.prime = this._p)
    (hfCanonical : SparsePolyZp.Canonical this._p.toNat f)
    (hfDegree : ∀ term ∈ f.toList, term.1.deg < 2 ^ 62)
    (hfMonic : (SparsePolyZp.toPoly this._p.toNat f).Monic)
    (hfSquarefree : Squarefree (SparsePolyZp.toPoly this._p.toNat f)) :
    ∃ output,
      Generated.StrictDDF.__ddf_Zp_raw_ir
          (StrictDDF.strictDDFRawOps this providers
            (SparsePolyZp.toPoly this._p.toNat f)) f
          (fun _ => StrictDDF.DDFLoopInvariant.initial this f
            f[0]!.2.prime hfPrime
            hfCanonical hfDegree hfMonic hfSquarefree) = .ok output ∧
      StrictDDF.ddfResultToL2 this._p.toNat output =
        ddf (SparsePolyZp.toPoly this._p.toNat f) ∧
      (∀ item ∈ output.toList,
        SparsePolyZp.Canonical this._p.toNat item.1) ∧
      (∀ item ∈ output.toList,
        0 < (SparsePolyZp.toPoly this._p.toNat item.1).natDegree) ∧
      ∀ item ∈ output.toList, 0 < item.2.toNat := by
  exact {proof_theorem} this providers f hfPrime
    hfCanonical hfDegree hfMonic hfSquarefree
"""
    if contract["kind"] == "strict_edf_zp":
        return f"""/-- Generated public contract for the original C++ `__edf_Zp` entry.
The executable side is the strict, well-founded L1 semantics, including its
exact RNG transition; every newly appended concrete factor satisfies the L2
`EDFCorrect` specification. -/
theorem {theorem_name}
    {{State : Type}} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (termination : Generated.StrictEDF.EDFTermination
      (StrictEDF.strictEDFRawOps engine this providers))
    (result : Array SparsePolyZp) (f : SparsePolyZp) (d : UInt64)
    (rng : State) (hinvariant : StrictEDF.EDFEntryInvariant this f d)
    (hresultCanonical : ∀ factor ∈ result.toList,
      SparsePolyZp.Canonical this._p.toNat factor) :
    ∃ output rng' factors,
      Generated.StrictEDF.__edf_Zp_raw_ir
          (StrictEDF.strictEDFRawOps engine this providers)
          (StrictEDF.strictEDFSplitLaw engine this providers) termination
          result f d rng hinvariant = .ok (output, rng') ∧
      StrictEDF.edfResultToL2 this._p.toNat output =
        StrictEDF.edfResultToL2 this._p.toNat result ++ factors ∧
      EDFCorrect (SparsePolyZp.toPoly this._p.toNat f) d.toNat factors ∧
      ∀ factor ∈ output.toList,
        SparsePolyZp.Canonical this._p.toNat factor := by
  exact {proof_theorem} engine this providers termination result f d rng
    hinvariant hresultCanonical
"""
    if contract["kind"] == "strict_hensel_eea":
        return f"""/-- Generated public contract for the polynomial extended
Euclidean entry used by Hensel tree construction. The executable side is the
strict degree-measured L1 loop, including generated division, inverse, and
terminal coefficient scaling; its first result is the L2 monic gcd and all
three concrete results satisfy the Bézout contract. -/
theorem {theorem_name}
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (left right : SparsePolyZp)
    (hleftCanonical : SparsePolyZp.Canonical this._p.toNat left)
    (hrightCanonical : SparsePolyZp.Canonical this._p.toNat right)
    (hleftNonempty : 0 < left.size) :
    ∃ gcd s t,
      Generated.StrictHensel.__polynomial_GCD_eea_raw_ir
          (StrictHensel.strictHenselEEARawOps this)
          (StrictHensel.strictHenselEEATermination this)
          (Generated.StrictHensel.henselEEAInitialState this._p left right) =
        .ok (gcd, s, t) ∧
      SparsePolyZp.toPoly this._p.toNat gcd =
        SparsePolyZp.toPoly this._p.toNat s *
            SparsePolyZp.toPoly this._p.toNat left +
          SparsePolyZp.toPoly this._p.toNat t *
            SparsePolyZp.toPoly this._p.toNat right ∧
      (SparsePolyZp.toPoly this._p.toNat gcd).Monic ∧
      SparsePolyZp.toPoly this._p.toNat gcd ∣
        SparsePolyZp.toPoly this._p.toNat left ∧
      SparsePolyZp.toPoly this._p.toNat gcd ∣
        SparsePolyZp.toPoly this._p.toNat right ∧
      SparsePolyZp.toPoly this._p.toNat gcd =
        GCDMonoid.gcd
          (SparsePolyZp.toPoly this._p.toNat left)
          (SparsePolyZp.toPoly this._p.toNat right) := by
  exact {proof_theorem} this hcfg left right
    hleftCanonical hrightCanonical hleftNonempty
"""
    if contract["kind"] == "strict_hensel_step":
        return f"""/-- Generated public contract for the original C++ `__hensel_step` entry.
The executable side is the exact strict raw L1 program, and its output
satisfies the L2 quadratic Hensel-step invariant without a model fallback. -/
theorem {theorem_name}
    (node : HenselNode) (f : SparsePolyZZ) (m : Nat)
    (hinvariant : StrictHensel.HenselStepRefinementInvariant
      StrictHensel.concreteDivmodTermination node f m) :
    ∃ output,
      Generated.StrictHensel.__hensel_step_raw_ir
          (StrictHensel.strictHenselRawOps
            StrictHensel.concreteDivmodTermination)
          node f (m : Int) = .ok output ∧
      StrictHensel.HenselStepCorrect f m node output := by
  exact {proof_theorem} StrictHensel.concreteDivmodTermination node f m hinvariant
"""
    if contract["kind"] == "strict_hensel_lift_recursive":
        return f"""/-- Generated public contract for the original C++
`__hensel_lift_recursive` entry.  The executable side performs the exact
top-down raw tree traversal with structural well-founded recursion; its
semantic trace proves every concrete node update is a quadratic Hensel step. -/
theorem {theorem_name}
    (tree : Generated.StrictHensel.HenselLiftTree)
    (nodes : Array HenselNode) (target : SparsePolyZZ) (m : Nat)
    (hinvariant : StrictHensel.HenselLiftRecursiveRefinementInvariant
      StrictHensel.concreteDivmodTermination tree nodes target m) :
    ∃ output,
      Generated.StrictHensel.__hensel_lift_recursive_raw_ir
          (StrictHensel.strictHenselRawOps
            StrictHensel.concreteDivmodTermination)
          tree nodes target (m : Int) = .ok output ∧
      StrictHensel.HenselLiftRecursiveCorrect
        StrictHensel.concreteDivmodTermination m tree nodes target output := by
  exact {proof_theorem} StrictHensel.concreteDivmodTermination tree nodes target m hinvariant
"""
    if contract["kind"] == "strict_hensel_extract_factors":
        return f"""/-- Generated public contract for the original C++
`__hensel_extract_factors` entry.  The executable side performs the exact
left-before-right traversal of the concrete Hensel tree; its semantic trace
records precisely which leaf factors are appended to the input array. -/
theorem {theorem_name}
    (tree : Generated.StrictHensel.HenselLiftTree)
    (nodes : Array HenselNode) (factors : Array SparsePolyZZ)
    (hinvariant : StrictHensel.HenselExtractInvariant tree nodes) :
    ∃ output,
      Generated.StrictHensel.__hensel_extract_factors_raw_ir
          tree nodes factors = .ok output ∧
      StrictHensel.HenselExtractCorrect tree nodes factors output := by
  exact {proof_theorem} tree nodes factors hinvariant
"""
    if contract["kind"] == "strict_hensel_tree_build":
        return f"""/-- Generated public contract for the original C++
`__hensel_tree_build` entry. The actual strict raw builder allocates exactly
the canonical preorder topology and produces the invariant consumed by the
generated extraction traversal; no mathematical tree is supplied as an
oracle. -/
theorem {theorem_name}
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    (factors : Array SparsePolyZp)
    (hfactors : ∀ factor ∈ factors.toList,
      SparsePolyZp.Canonical this._p.toNat factor)
    (hfactorsNonempty : ∀ factor ∈ factors.toList, 0 < factor.size)
    (hfactorsMonicAfterZero : ∀ index (hindex : index < factors.size),
      0 < index → (SparsePolyZp.toPoly this._p.toNat
        (getElem factors index hindex)).Monic)
    (hpairwise : ∀ i j (hi : i < factors.size) (hj : j < factors.size),
      i < j → IsCoprime
        (SparsePolyZp.toPoly this._p.toNat (getElem factors i hi))
        (SparsePolyZp.toPoly this._p.toNat (getElem factors j hj)))
    (htwo : 2 ≤ factors.size)
    (hfactorFits : factors.size < 2 ^ 31)
    (hfitsInt32 : StrictHensel.henselTreeInternalNodeCount
      0 factors.size < 2 ^ 31) :
    let tree := StrictHensel.henselTreeBuildTopology 0 factors.size 0
    ∃ output,
      Generated.StrictHensel.__hensel_tree_build_raw_ir
          (StrictHensel.strictHenselTreeBuildRawOps this mulProvider)
          factors this._p = .ok output ∧
      ∃ hroot : 0 < output.size,
      output.size = tree.nodeCount ∧ tree.rootIndex = 0 ∧
      StrictHensel.liftChildMatches (getElem output 0 hroot).left
        (match tree with | .node _ left _ => left) ∧
      StrictHensel.liftChildMatches (getElem output 0 hroot).right
        (match tree with | .node _ _ right => right) ∧
      StrictHensel.HenselExtractInvariant tree output ∧
      StrictHensel.HenselTreeSemanticBuildCertificate this._p.toNat factors 0
        0 factors.size tree output ∧
      StrictHensel.HenselTreeNodeInitialInvariant this._p.toNat factors
        0 factors.size (getElem output 0 hroot) ∧
      StrictHensel.HenselArrayCanonical output ∧
      StrictHensel.HenselArrayHOneHead output := by
  exact {proof_theorem} this hcfg mulProvider factors hfactors
    hfactorsNonempty hfactorsMonicAfterZero hpairwise htwo hfactorFits
      hfitsInt32
"""
    if contract["kind"] == "strict_hensel_lift_upoly":
        return f"""/-- Generated public contract for the original C++
`__hensel_lift_upoly` entry.  The strict generated L1 program executes target
selection, coefficient adjustment, tree construction, quadratic lifting,
leaf extraction, and normalization; the L2 trace records every actual stage
result without an oracle, fallback, or fuel parameter. -/
theorem {theorem_name}
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    (f : SparsePolyZZ) (factors : Array SparsePolyZp) (aTarget : Int32)
    (hinvariant : StrictHensel.HenselLiftEntryInvariant this
      StrictHensel.concreteDivmodTermination mulProvider f factors aTarget)
    (hfactorCount : 2 ≤ factors.size)
    (hfactorFits : factors.size < 2 ^ 31) :
    ∃ output,
      Generated.StrictHensel.__hensel_lift_upoly_raw_ir
          (StrictHensel.strictHenselRawOps
            StrictHensel.concreteDivmodTermination)
          (StrictHensel.strictHenselTreeBuildRawOps this mulProvider)
          f factors this._p aTarget
            (Nat.Prime.two_le (Fact.out : Nat.Prime this._p.toNat)) =
              .ok output ∧
      StrictHensel.HenselLiftEntryCorrect
        StrictHensel.concreteDivmodTermination f factors this._p aTarget
        output := by
  exact {proof_theorem} this hcfg
    StrictHensel.concreteDivmodTermination mulProvider f factors aTarget
      hinvariant hfactorCount hfactorFits
"""
    raise ValueError(f"unsupported verified contract kind: {contract['kind']}")


def _emit_verified_contracts_file(infos: list[dict]) -> str:
    """Emit the single discoverable module containing every closed contract."""
    contracts = [info for info in infos if "verified_contract" in info]
    imports = sorted({info["verified_contract"]["proof_import"] for info in contracts})
    lines = [
        "-- Auto-generated by cpp2lean v2 Pass 9 (refine_gen)",
        "-- Do not edit: all completed public L1 → L2 contracts live here.",
        "",
        *(f"import {module}" for module in imports),
        "",
        "set_option autoImplicit false",
        "set_option linter.style.nameCheck false",
        "",
        "open Polynomial",
        "open CLPoly.Math",
        "",
        "namespace Refinement",
        "",
    ]
    for info in contracts:
        lines.extend((_emit_verified_contract(info).rstrip(), ""))
    lines.extend(("end Refinement", ""))
    return "\n".join(lines)


# ============================================================
# Pass 9 入口
# ============================================================

def refine_gen_pass(top_funcs: Optional[list[MIRFunc]] = None) -> None:
    mir_index: dict[str, MIRFunc] = {}
    if top_funcs is not None:
        for f in top_funcs:
            mir_index[f.base_name] = f

    # Public refinement output is deliberately centralized.  Incomplete
    # per-algorithm skeletons used to obscure which contracts were actually
    # proved and introduced `sorry` declarations beside the checked API.
    out_path = REFINEMENT_DIR / "Generated.lean"
    out_path.write_text(_emit_verified_contracts_file(list(REFINEMENT_MAP.values())))
    print(f"  Verified refinements: {out_path}", file=sys.stderr)


def main():
    import json
    from pass1_parse import parse_pass
    from pass2_ref_elim import ref_elim_pass
    from pass2b_callsite_ref_elim import callsite_ref_elim_pass
    from pass3_lambda_lift import lambda_lift_pass
    from pass3b_lambda_ref_elim import lambda_ref_elim_pass
    from pass4_iter_recognize import iter_recognize_pass
    from pass5_operator_resolve import operator_resolve_pass
    from pass6_ssa_build import ssa_build_pass
    from pass7_loop_lower import loop_lower_pass

    ast_cache = Path(__file__).resolve().parent.parent / "tests" / "ast_cache"
    targets = [k for k in REFINEMENT_MAP if not k.startswith("_loop_")]
    mirs = []
    for fn in targets:
        cache = ast_cache / f"{fn}.json"
        if not cache.exists():
            print(f"  skip: {fn} (no AST cache)", file=sys.stderr)
            continue
        try:
            with open(cache) as f:
                ast = json.load(f)
            h = lambda_ref_elim_pass(lambda_lift_pass(callsite_ref_elim_pass(
                ref_elim_pass(parse_pass(ast)))))
            h3 = iter_recognize_pass(h)
            h4, _ = operator_resolve_pass(h3)
            h4 = callsite_ref_elim_pass(h4, callee_filter={"Rng.next_advance"})
            mir0 = ssa_build_pass(h4)
            mir1 = loop_lower_pass(mir0)
            mirs.append(mir1)
        except Exception as e:
            print(f"  error: {fn} ({type(e).__name__}: {e})", file=sys.stderr)
    refine_gen_pass(mirs)


if __name__ == "__main__":
    main()
