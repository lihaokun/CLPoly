# Strict select-prime candidate proof audit

## Scope

This note covers the refinement of one successful execution of the generated
`__select_prime` candidate body.  It does not claim the outer prime-search
loop complete.

## Natural-language proof draft

1. Execute the generated leading-coefficient guard and strict integer
   polynomial reduction.  The reduction theorem identifies the resulting L1
   sparse polynomial with `Polynomial.map (Int.castRingHom (ZMod p)) f` and
   proves its canonical representation.
2. The generated degree and derivative guards rule out the empty and zero
   derivative cases.  Apply the existing strict derivative refinement and the
   strict heap-GCD refinement to the exact values produced by the raw body.
3. The successful GCD guard says the canonical normalized GCD does not have
   positive degree.  Its monicity and nonzeroness therefore force degree zero;
   hence the source polynomial and its derivative are coprime and the modular
   polynomial is squarefree.
4. Apply the strict make-monic refinement to the raw result, preserving
   canonicality, degree, and association with the pre-normalized polynomial.
5. Apply the already proved strict DDF call theorem, then induct over the
   generated DDF-component loop.  For each concrete component, invoke the
   strict EDF call theorem and induct over the generated append loop.  This
   proves that the actual returned array consists of monic irreducibles and
   that its product is associated with the mapped input polynomial.
6. Decode the constructor equality from the raw successful return; this
   identifies the exact L1 `fp`, factor array, and RNG state with the values in
   the proof.  Package `GoodPrime`, product association, and factor quality as
   `CandidateCorrect`.

## Type-level modulus issue

The original `CandidatePhysical p` stored a `DenseUPolyZp` plus an equality
`dense._p = p`.  Since `ZMod dense._p.toNat` and `ZMod p.toNat` carry dependent
field instances, repeated transport across that equality was not robust.  It
was replaced by `CandidateArithmetic p`, whose reconstruction assigns `_p :=
p` definitionally.  The derivative/GCD/DDF/EDF chain now stays in one field
without casts.

The first loop theorem also retained only irreducibility and monicity of EDF
outputs.  That was insufficient to recover the factor product.  Its invariant
now carries the exact product of the unvisited concrete DDF suffix, derived at
each step from the actual `EDFCorrect` result.

The outer loop additionally needs a prime-iterator invariant: the initial
word is prime, `nextPrime`/`prevPrime` preserve primality, every physical
candidate is requested only with that proof, and a positive `tried` count
implies that `best` satisfies `CandidateCorrect`.

## 度量

- 耗时：~1.5 小时（接口审计、形式化调试、证明草稿）
- 迭代：12 轮编译-修复循环
- Lean 新增/修改行数：约 190 行净修改（含加强的循环不变量）
- 对应 C++ 行数：约 55 行候选执行控制流
- 放弃的方案：在每个中间定理处用 `rw [primeField]` 搬运 `ZMod`
  类型；依赖的域实例使该方案脆弱且部分 motive 不合法
