# Strict polynomial symmetric-mod refinement

The van-Hoeij candidate validator no longer delegates C++
`__upoly_symmetric_mod` to a raw callback.  Its generator emits the concrete
range-for loop that symmetrically reduces each integer coefficient, omits
zeros, and advances with measure `input.size - index`.

The refinement proof follows that loop term by term.  It proves from integer
division/remainder identities that `ZZ.symmetricMod c m` casts to the same
element of `ZMod m` as `c`, handles the zero-omission branch explicitly, and
derives that every successful `symmetricModRaw` execution represents the same
polynomial modulo the positive modulus.

There is no fuel, partial definition, admitted lemma, semantic oracle, or
factorization witness in this boundary.

Estimated FactorZZ progress: about 93%.  The next concrete callees are
primitive-part normalization and integer sparse trial division, followed by
the full candidate/factor-product invariant and the public FactorZZ pipeline.
Estimated remaining focused time: 1--2.5 workdays.
