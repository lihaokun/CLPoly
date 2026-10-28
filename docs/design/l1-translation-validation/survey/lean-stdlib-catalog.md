# Lean 4 Stdlib & Mathlib API Reference for cpp2lean v2

**Version**: Lean 4.28.0 + Mathlib4 + Batteries  
**Date**: 2026-04-19  
**Purpose**: Exact API specification for C++ ↔ Lean 4 type and operation mapping in CLPoly

---

## Overview

This catalog lists exact type signatures, UB conditions, and usage notes for types/operations needed to translate CLPoly C++ code to Lean 4 L1 IR. Each entry includes:

- **Full Path**: Module path in stdlib/Mathlib/Batteries
- **Type Parameters**: Implicit/explicit universe levels
- **Core Operations**: Signature + whether computable
- **UB Conditions**: Panic/undefined behavior on edge cases
- **Frequency**: How often used in typical C++ translation

---

## 1. Basic Numeric Types

### 1.1 Nat (Natural Numbers)

| Field | Value |
|-------|-------|
| **Full Path** | `Init.Data.Nat` (stdlib core) |
| **Type Constructor** | `Nat : Type` |
| **Universe** | Sort 0 (not a Type u) |
| **Instances** | `BEq`, `Hashable`, `Ord`, `LinearOrder` |

**Core Operations:**

| Op | Signature | Computable | Notes |
|----|-----------|-----------|-------|
| `.add` | `Nat → Nat → Nat` | Yes | No overflow (unbounded) |
| `.sub` | `Nat → Nat → Nat` | Yes | Saturating: `a - b = 0` if `a < b` |
| `.mul` | `Nat → Nat → Nat` | Yes | No overflow |
| `.div` | `Nat → Nat → Nat` | Yes | **Divides by 0 → 0** (not panic) |
| `.mod` | `Nat → Nat → Nat` | Yes | `a % 0 = a` (not panic) |
| `.pow` | `Nat → Nat → Nat` | Yes | Exponentiation, no overflow |
| `.max`, `.min` | `Nat → Nat → Nat` | Yes | Binary max/min |
| `.log` | `Nat → Nat → Nat` | Yes | Integer log base (Mathlib) |
| `.sqrt` | `Nat → Nat` | Yes | Integer square root |

**UB Conditions:**
- **Division by zero**: `a / 0 = 0` (NOT panic; differs from C++ UB)
- **Modulo by zero**: `a % 0 = a` (not panic)
- **Subtraction underflow**: Saturates to 0 (not UB)
- **Overflow**: None; Nat is unbounded

**Comparison:**
- `.< , .<= , .> , .>= : Nat → Nat → Bool` (decidable)
- `.compare : Nat → Nat → Ordering` (total order)

**Mathlib-specific operations** (in `Data.Nat.*`):
- `Nat.Prime : Nat → Prop` (irreducible, for factorization)
- `Nat.factorization : Nat → Multiset ℕ` (prime factorization)
- `Nat.log b n : Nat` (integer log base; `Nat.log 2 8 = 3`)

---

### 1.2 UInt64 (Unsigned 64-bit Integer)

| Field | Value |
|-------|-------|
| **Full Path** | `Batteries.Data.UInt` (Batteries library) |
| **Type Constructor** | `UInt64 : Type` |
| **Universe** | Sort 0 |
| **Instances** | `BEq`, `Ord`, `Hashable` |
| **Representation** | Wrapper around `{ n : Nat // n < 2^64 }` |

**Core Operations:**

| Op | Signature | Computable | Notes |
|----|-----------|-----------|-------|
| `.toNat` | `UInt64 → Nat` | Yes | Unwrap to `Nat` (lossless) |
| `.ofNat` | `Nat → UInt64` | Yes | Wraps `n % 2^64` |
| `+` | `UInt64 → UInt64 → UInt64` | Yes | **Wraps on overflow** (mod 2^64) |
| `-` | `UInt64 → UInt64 → UInt64` | Yes | **Wraps on underflow** (mod 2^64) |
| `*` | `UInt64 → UInt64 → UInt64` | Yes | **Wraps on overflow** (mod 2^64) |
| `/` | `UInt64 → UInt64 → UInt64` | Yes | **Divides by 0 → 0** (not panic) |
| `%` | `UInt64 → UInt64 → UInt64` | Yes | `a % 0 = a` (not panic) |
| `.shiftLeft (<<)` | `UInt64 → Nat → UInt64` | Yes | **Shift ≥ 64 → 0** (wraps) |
| `.shiftRight (>>)` | `UInt64 → Nat → UInt64` | Yes | **Shift ≥ 64 → 0** (wraps) |

**UB Conditions:**
- **Division by zero**: `a / 0 = 0` (NOT panic; **differs from C++ undefined behavior**)
- **Modulo by zero**: `a % 0 = a`
- **Overflow/Underflow**: Wraps silently (mod 2^64)
- **Shift by ≥ 64**: Result is 0 (not panic)

**Equivalence to Nat:**
```lean
x.toNat < 2^64              -- Always true by construction
(UInt64.ofNat n).toNat = n % 2^64
UInt64.ext : x.toNat = y.toNat → x = y  -- Extensionality
```

**Comparison:**
- `x ≤ y ↔ x.toNat ≤ y.toNat` (`UInt64.le_iff_toNat_le_toNat`)
- `x < y ↔ x.toNat < y.toNat` (`UInt64.lt_iff_toNat_lt_toNat`)
- `.compare : UInt64 → UInt64 → Ordering`

**Conversion:**
- `.toInt64 : UInt64 → Int64` (via `.toNat`)
- `.toNat : UInt64 → Nat` (lossless)

**CRITICAL FOR C++ TRANSLATION:**
- Lean's `/` and `%` on `UInt64` do **NOT panic** on zero (C++ UB). To match C++ semantics, require `b ≠ 0` via explicit hypothesis or check in guards.
- Shifts ≥ 64 are well-defined in Lean (result = 0), not UB as in C++.

---

### 1.3 Int64 / Int (Signed Integers)

| Field | Value |
|-------|-------|
| **Full Path** | `Init.Data.Int` (stdlib) or `Batteries.Data.UInt` (signed variants) |
| **Type Constructor** | `Int : Type` (unbounded) or `Int64` (bounded, limited support) |
| **Universe** | Sort 0 |
| **Instances** | `BEq`, `Ord`, `Ring`, `LinearOrder` |

**Core Operations:**

| Op | Signature | Computable | Notes |
|----|-----------|-----------|-------|
| `.add` | `Int → Int → Int` | Yes | No overflow (unbounded) |
| `.sub` | `Int → Int → Int` | Yes | No overflow |
| `.mul` | `Int → Int → Int` | Yes | No overflow |
| `.div` | `Int → Int → Int` | Yes | **Divides by 0 → 0** (not panic; Euclidean division) |
| `.mod` | `Int → Int → Int` | Yes | Euclidean: `a % 0 = a` |

**UB Conditions:**
- **Division by zero**: `a / 0 = 0` (not panic)
- **No overflow**: `Int` is unbounded in Lean (unlike C++ `int64_t`)

**Comparison to C++ int64_t:**
- Lean's `Int` is arbitrary precision, not 2's complement
- For C++ semantics, use `UInt64` + bitwise ops + manual wrapping

---

### 1.4 Bool

| Field | Value |
|-------|-------|
| **Full Path** | `Init.Data.Bool` (stdlib) |
| **Type Constructor** | `Bool : Type` with `true`, `false` |
| **Universe** | Sort 0 |
| **Instances** | `BEq`, `Decidable` |

**Core Operations:**

| Op | Signature | Short-Circuit? | Notes |
|----|-----------|----------------|-------|
| `&&` (and) | `Bool → Bool → Bool` | **YES** (by pattern match in Lean 4) | Conjunction; `false && x` never evaluates `x` |
| `\|\|` (or) | `Bool → Bool → Bool` | **YES** (by pattern match in Lean 4) | Disjunction; `true \|\| x` never evaluates `x` |
| `!` (not) | `Bool → Bool` | N/A | Logical negation |
| `^^` (xor) | `Bool → Bool → Bool` | **NO** (both sides evaluated always) | Exclusive or; no short-circuit |
| `.decide` | `Decidable p → Bool` | Depends on `p` | Propositional to Bool |

**UB Conditions:**
- None; all operations are total and decidable

**Key Lemmas:**
- `Bool.and_comm : (x && y) = (y && x)`
- `Bool.or_comm : (x || y) = (y || x)`
- `Bool.not_and_self : (!x && x) = false`
- `Bool.and_eq_left_iff_imp : ((a && b) = a) ↔ (a → b)`

**Short-Circuit Behavior in Lean:**
```lean
-- In Lean 4, && and || are syntactic sugar for:
def and : Bool → Bool → Bool
  | false, _ => false
  | true, b => b

def or : Bool → Bool → Bool
  | true, _ => true
  | false, b => b

-- Both provably short-circuit (right operand not evaluated if left matches)
```

---

### 1.5 Float (Double Precision)

| Field | Value |
|-------|-------|
| **Full Path** | `Init.Data.Float` (stdlib) |
| **Type Constructor** | `Float : Type` |
| **Universe** | Sort 0 |
| **Instances** | `BEq`, `Ord` |

**Core Operations:**

| Op | Signature | Notes |
|----|-----------|-------|
| `.add` | `Float → Float → Float` | IEEE 754 addition |
| `.sub` | `Float → Float → Float` | IEEE 754 subtraction |
| `.mul` | `Float → Float → Float` | IEEE 754 multiplication |
| `.div` | `Float → Float → Float` | IEEE 754 division; `/0 = inf` or `nan` |
| `.sqrt` | `Float → Float` | IEEE 754 square root |

**UB Conditions:**
- **Division by zero**: `f / 0.0 = inf` (or `nan` if `f = 0`), not panic
- NaN propagation follows IEEE 754 semantics

**Comparison:**
- `.< , .<= , .> , .>= : Float → Float → Bool`
- `.compare : Float → Float → Ordering` (follows IEEE; NaN < everything)

---

## 2. Container Types

### 2.1 Array α

| Field | Value |
|-------|-------|
| **Full Path** | `Init.Data.Array` (stdlib core) |
| **Type Constructor** | `Array α : Type u → Type u` |
| **Universe** | `α : Type u` ⟹ `Array α : Type u` |
| **Instances** | `Inhabited`, `BEq` (if `BEq α`), `Membership α (Array α)` |
| **Underlying** | Immutable vector backed by persistent data structure |

**Core Operations:**

| Op | Signature | Computable | Panic Condition | Notes |
|----|-----------|-----------|-----------------|-------|
| `.size` | `Array α → Nat` | Yes | — | O(1) |
| `.push` | `Array α → α → Array α` | Yes | — | Append to end; O(1) amortized |
| `.pop` | `Array α → Array α` | Yes | — | Remove last; empty → empty |
| `.get i` | `Array α → Nat → (h : i < size) → α` | Yes | None (requires proof) | Bounded access; O(1) |
| `.get! [Inhabited α]` | `Array α → Nat → α` | Yes | **Out of bounds** | Return `default` on failure (not panic) |
| `.get?` | `Array α → Nat → Option α` | Yes | — | Safe access; `none` if out of bounds |
| `.set i v` | `Array α → Nat → α → (h : i < size) → Array α` | Yes | None (requires proof) | Update i-th element |
| `.set! [Inhabited α]` | `Array α → Nat → α → Array α` | Yes | Out of bounds (silent fail) | No-op if out of bounds |
| `.swap i j` | `Array α → Nat → Nat → (hi hj : bounds) → Array α` | Yes | None (requires proofs) | Swap two elements |
| `.swapIfInBounds` | `Array α → Nat → Nat → Array α` | Yes | — | Safe swap; no-op if either out of bounds |
| `.map f` | `(α → β) → Array α → Array β` | Yes | — | Pointwise map |
| `.filter p` | `(α → Bool) → Array α → Array α` | Yes | — | Keep elements where `p` is true |
| `.foldl f init` | `(β → α → β) → β → Array α → β` | Yes | — | Left fold; `foldl f b [a₁,a₂,...] = f(f(b, a₁), a₂)...` |
| `.foldr f init` | `(α → β → β) → Array α → β → β` | Yes | — | Right fold |
| `.range n` | `Nat → Array Nat` | Yes | — | `Array.range n = #[0, 1, ..., n-1]` |
| `.take n` | `Array α → Nat → Array α` | Yes | — | First `n` elements; extra indices OK |
| `.drop n` | `Array α → Nat → Array α` | Yes | — | Skip first `n` elements |
| `.isEmpty` | `Array α → Bool` | Yes | — | `size = 0` |

**UB Conditions:**
- `.get i h`: Requires `i < size` proof; **no panic** if proof is correct
- `.get! i`: **Does not panic on out of bounds**; returns `default` (or `panic` if no `Inhabited` instance available; see below)
- `.get? i`: Safe; returns `Option`
- `.set! i v`: **Silent no-op if `i ≥ size`**; doesn't panic in Lean, but may differ from C++ behavior
- `.pop` on empty: Returns empty (not panic)
- **Array size limit**: Documented to support up to `USize.size` elements; exceeding is UB in Lean semantics (but practically unchecked)

**Key Lemmas:**
```lean
Array.size_push : (xs.push v).size = xs.size + 1
Array.size_pop : xs.pop.size = xs.size - 1  (xs.empty.size - 1 = 0)
Array.size_set : (xs.set i v h).size = xs.size
Array.getElem_push_lt (h : i < xs.size) : (xs.push v)[i] = xs[i]
Array.getElem_push_eq : (xs.push v)[xs.size] = v
Array.mem_def : a ∈ as ↔ a ∈ as.toList  -- Array membership is via toList
```

**Conversion to List:**
- `.toList : Array α → List α` (O(n), creates list from array)
- `List.toArray : List α → Array α` (O(n))

**Sorting Operations (Batteries extension):**

| Op | Signature | Computable | Notes |
|----|-----------|-----------|-------|
| (No `qsort` in Lean stdlib) | — | — | Use `Batteries` or implement custom |
| (Lean 4.28.0 stdlib: no direct quicksort) | — | — | InsertionSort available in `Init.Data.Array.InsertionSort` |

**Array in Batteries** (`Batteries.Data.Array.Basic`):
- `.minWith`, `.minD`, `.min?`, `.minI`
- `.maxWith`, `.maxD`, `.max?`, `.maxI`
- Additional sorting/filtering lemmas

---

### 2.2 List α

| Field | Value |
|-------|-------|
| **Full Path** | `Init.Data.List.Basic` (stdlib) |
| **Type Constructor** | `List α : Type u → Type u` |
| **Universe** | `α : Type u` ⟹ `List α : Type u` |
| **Instances** | `Inhabited` (empty list), `Append`, `Functor`, `Monad` |
| **Structure** | Singly-linked list (recursive inductive) |

**Core Operations:**

| Op | Signature | Notes |
|----|-----------|-------|
| `.length` | `List α → Nat` | O(n) |
| `.head!` | `List α → [Inhabited α] → α` | First element; default if empty |
| `.tail` | `List α → List α` | Rest of list; empty → empty |
| `.cons a l` | `α → List α → List α` | Prepend; `a :: l` |
| `.append l₁ l₂` | `List α → List α → List α` | Concatenate; `l₁ ++ l₂` |
| `.map f` | `(α → β) → List α → List β` | Pointwise map |
| `.filter p` | `(α → Bool) → List α → List α` | Keep elements where `p` true |
| `.foldl f b` | `(β → α → β) → β → List α → β` | Left fold |
| `.foldr f b` | `(α → β → β) → List α → β → β` | Right fold |
| `.reverse` | `List α → List α` | Reverse; O(n) |
| `.sort cmp` | `(α → α → Ordering) → List α → List α` | Insertion sort |

**Relation to Array:**
- `.toArray : List α → Array α` (cheap in practice)
- `Array.toList : Array α → List α` (cheap)
- `.toList : Array α → List α` is definitionally the list storage

---

### 2.3 Prod α β (Product / Pair)

| Field | Value |
|-------|-------|
| **Full Path** | `Init.Data.Prod` (stdlib) |
| **Type Constructor** | `Prod α β : Type u → Type v → Type (max u v)` |
| **Notation** | `α × β` for `Prod α β` |
| **Instances** | `Prod.Decidable`, `Inhabited` (if both components inhabited) |
| **Structure** | Dependent pair `⟨fst : α, snd : β⟩` |

**Core Operations:**

| Op | Signature | Notes |
|----|-----------|-------|
| `.fst` | `Prod α β → α` | First component |
| `.snd` | `Prod α β → β` | Second component |
| `.mk a b` | `α → β → Prod α β` | Constructor; `⟨a, b⟩` |
| `⟨a, b⟩` | Notation for `.mk` | Tuple literal |
| `(a, b)` | Alternative notation | Same as `Prod.mk a b` |

**Extensionality:**
```lean
Prod.ext : ∀ {x y : Prod α β}, x.fst = y.fst → x.snd = y.snd → x = y
```

**Conversion:**
- `Prod.toLex : Prod α β → lex(α, β)` (lexicographic order)

---

## 3. Algorithmic Operations

### 3.1 Comparison and Ordering

**Full Path**: `Init.Cmp`, `Data.Ordering`

| Type | Signature | Notes |
|------|-----------|-------|
| `Ordering` | `LT \| EQ \| GT` | Three-way comparison result |
| `compare : α → α → Ordering` | Total comparison (requires `Ord` instance) | Reflexive, transitive, antisymmetric |
| `(· < ·) : α → α → Bool` | Less-than (requires `Ord`) | Decidable on concrete types |
| `(· ≤ ·) : α → α → Bool` | Less-or-equal | Decidable |
| `max, min : α → α → α` | Maximum/minimum (requires `Ord`) | Computable |

**Instances:**
- `Nat`, `Int`, `UInt64`, `Float` all have `Ord` instances
- Array/List elements' comparisons depend on element's `Ord`

**Key Lemmas:**
```lean
Nat.lt_irrefl : ¬(a < a)
Nat.lt_trans : a < b → b < c → a < c
```

---

### 3.2 Min/Max Functions

**Full Path**: `Init.Order.Defs` or `Batteries.Data.Array.Basic` (for array variants)

| Function | Signature | Notes |
|----------|-----------|-------|
| `max : Ord α ⟹ α → α → α` | Binary maximum | Computable |
| `min : Ord α ⟹ α → α → α` | Binary minimum | Computable |
| `Array.minI` | `Ord α, Inhabited α ⟹ Array α → α` | Minimum of array; default if empty |
| `Array.maxI` | `Ord α, Inhabited α ⟹ Array α → α` | Maximum of array; default if empty |
| `Array.minD` | `Ord α ⟹ Array α → α → α` | Minimum with explicit default |
| `Array.maxD` | `Ord α ⟹ Array α → α → α` | Maximum with explicit default |

---

### 3.3 Sorting (Batteries Extensions)

**Full Path**: No quicksort in Lean stdlib; Batteries has sorting infrastructure

| Function | Signature | Notes |
|----------|-----------|-------|
| (Lean stdlib has insertion sort) | `InsertionSort α : Ord α` | O(n²) but simple |
| (Recommend custom quicksort) | Custom `def qsortWith : (α → α → Bool) → Array α → Array α` | Must implement for C++ `qsort` translation |

**Recommended approach:**
```lean
-- Define your own quicksort with Lean termination proof
def qsortWith [Inhabited α] (cmp : α → α → Bool) : Array α → Array α :=
  -- Pivot selection + partition + recursive calls
  -- Must provide `termination_by Array.size` proof
```

---

## 4. Type Conversions

### 4.1 Nat ↔ UInt64

| Direction | Function | Signature | Notes |
|-----------|----------|-----------|-------|
| `Nat → UInt64` | `UInt64.ofNat` | `Nat → UInt64` | Wraps `n % 2^64` |
| `UInt64 → Nat` | `UInt64.toNat` or `.toNat` | `UInt64 → Nat` | Lossless; `UInt64.ofNat n = UInt64.ofNat (m)` iff `n ≡ m (mod 2^64)` |

**Key Lemmas:**
```lean
UInt64.ofNat_toNat : (x.toNat).ofNat = x  -- Round-trip OK
UInt64.toNat_ofNat : (UInt64.ofNat n).toNat = n % 2^64
```

### 4.2 Nat ↔ Int

| Direction | Function | Notes |
|-----------|----------|-------|
| `Nat → Int` | `Int.ofNat` | Injection; `ofNat n = n` as integer |
| `Int → Nat` | `Int.natAbs` | Absolute value as Nat; `natAbs (-5) = 5` |

### 4.3 Array ↔ List

| Direction | Function | Complexity | Notes |
|-----------|----------|-----------|-------|
| `Array α → List α` | `Array.toList` | O(n) | Creates list from array |
| `List α → Array α` | `List.toArray` | O(n) | Creates array from list |

---

## 5. Standard Library Functions (Mathlib Scope)

### 5.1 gcd, lcm (in Mathlib or Batteries)

**Full Path**: `Mathlib.Algebra.EuclideanDomain.Basic` or `Batteries.Data.Nat.Gcd`

| Function | Signature | Notes |
|----------|-----------|-------|
| `Nat.gcd : Nat → Nat → Nat` | GCD of two naturals | Euclidean algorithm |
| `Nat.lcm : Nat → Nat → Nat` | LCM of two naturals | `lcm a b = a * b / gcd a b` |

**Key Lemmas:**
```lean
Nat.gcd_comm : gcd a b = gcd b a
Nat.dvd_gcd_iff : (gcd a b ∣ c) ↔ (a ∣ c) ∧ (b ∣ c)
```

### 5.2 Decidable Predicates (Classical)

**Full Path**: `Init.Core`, `Batteries`

| Function | Signature | Notes |
|----------|-----------|-------|
| `Nat.Prime : Nat → Prop` | Primality (Mathlib) | Noncomputable; requires `Fact (Nat.Prime p)` |
| `Decidable.decide : Decidable p → Bool` | Extract Bool from proof | Works for decidable propositions |
| `Classical.decPred : Decidable p → DecidablePred p` | Non-constructive version | For `Nat.find` on non-decidable predicates |

---

## 6. Critical UB & Edge Case Summary

### Division and Modulo

| Operation | Lean Behavior | C++ Behavior | Translation Note |
|-----------|---------------|--------------|------------------|
| `Nat / 0` | Returns `0` | **UB** | Require `b ≠ 0` guard in C++ translation |
| `UInt64 / 0` | Returns `0` | **UB** | Require `b ≠ 0` guard |
| `Int / 0` | Returns `0` | **UB** | Require `b ≠ 0` guard |
| `Nat % 0` | Returns `a` | **UB** | Require `b ≠ 0` guard |
| `UInt64 % 0` | Returns `a` | **UB** | Require `b ≠ 0` guard |
| `UInt64 << n` (n ≥ 64) | Returns `0` | **UB** (often undefined) | Document that shifts are well-defined; equiv to `n % 64` (some archs) |
| `UInt64 >> n` (n ≥ 64) | Returns `0` | **UB** | Document behavior |

### Array Access

| Operation | Lean Behavior | C++ Behavior | Translation Note |
|-----------|---------------|--------------|------------------|
| `.get i` (out of bounds) | Requires proof of `i < size` | **UB** (buffer overflow) | Always use bounded `.get i h` in proofs; `.get! i` (unsafe) returns default |
| `.get! i` (out of bounds) | Returns `default` | **UB** | Marked unsafe; not panic, but returns default |
| `.set! i v` (out of bounds) | No-op (silently ignored) | **UB** (buffer overflow) | Check bounds in C++ translation; Lean's `.set! i v` doesn't panic or fail |
| `.pop` on empty | Returns empty (no-op) | **UB** (undefined behavior) | Always safe in Lean |

### Type Wraparound

| Operation | Lean Behavior | C++ uint64_t | Translation Note |
|-----------|---------------|--------------|------------------|
| `UInt64.max + 1` | Wraps to `0` | **Undefined behavior** | Lean wraps silently; C++ is UB |
| `UInt64.ofNat n` | `n % 2^64` | Implicit truncation | Matches C++ cast semantics |

---

## 7. Computable vs. Noncomputable

| Type/Function | Computable? | Notes |
|---------------|-----------|-------|
| `Nat`, `UInt64`, `Int` (basic ops) | Yes | All primitive operations compile to Lean runtime |
| `Bool && , \|\|` | Yes | Short-circuit evaluation in Lean 4 syntax |
| `Array.map`, `.filter`, `.foldl` | Yes | Pure functions; no IO required |
| `Array.range` | Yes | Generates array 0..n-1 |
| `Nat.Prime` | No | Requires proof; use `Fact (Nat.Prime p)` premise |
| `GCDMonoid.gcd` (Mathlib) | No | Noncomputable; use `Nat.gcd` for Nat |
| Sorting (custom) | Yes | If termination proof provided via `termination_by` |

---

## 8. Usage Frequency in CLPoly Translation

| Type/Op | Frequency | Use Case |
|---------|-----------|----------|
| `Nat` | **High** | Loop counters, array indices, degree |
| `UInt64` | **High** | C++ `uint64_t` coefficients in polynomials |
| `Array` | **Very High** | Polynomial coefficient storage, factor lists |
| `Bool` | **Very High** | Conditional tests, loop guards |
| `Prod` | **Medium** | Returning multiple values (gcd, quotient) |
| `List` | **Low** | Mostly use Array; List for proofs |
| `Int` | **Low** | C++ `int` arguments; prefer `Nat` or `UInt64` |
| Float | **Low** | Not used in integer factorization; skip |
| `max`/`min` | **Medium** | Degree bounds, loop limits |
| `gcd`/`lcm` | **Medium** | Polynomial GCD (Euclidean algorithm) |

---

## 9. Key Imported Modules

**Lean 4.28.0 Stdlib:**
```lean
import Init.Data.Nat
import Init.Data.Int
import Init.Data.Bool
import Init.Data.Array
import Init.Data.List
import Init.Data.Prod
import Init.Cmp
```

**Batteries (optional but recommended):**
```lean
import Batteries.Data.UInt
import Batteries.Data.Array.Basic
import Batteries.Data.Array.Lemmas
```

**Mathlib (for advanced operations):**
```lean
import Mathlib.Data.Nat.Factorization
import Mathlib.Algebra.EuclideanDomain.Basic
```

---

## 10. Decision Tree: Which Type for C++ ↔ Lean

**C++ `uint64_t`** → Lean `UInt64`
- Operations: `+`, `-`, `*`, `/`, `%`, `<<`, `>>`
- All wrap on overflow (match C++ behavior)
- **Note:** `/0` is Lean-safe (returns 0), C++ UB; add guards

**C++ `int`** → Lean `Int` or `Nat`
- Prefer `Nat` if always positive (indices, counts)
- Use `Int` for signed arithmetic (rare in CLPoly)

**C++ `std::vector<T>`** → Lean `Array α`
- Immutable in Lean; mutations via functional updates
- `.push`, `.pop`, `.size` map directly
- `.get!` is **not** panic-safe; use `(i : Nat) → (h : i < size) → T` for proofs

**C++ `std::pair<A, B>`** → Lean `Prod α β` or `α × β`
- Notation: `⟨a, b⟩`
- Access: `.fst`, `.snd`

**C++ `if (a && b)`** → Lean `if a && b then ...`
- Short-circuits in Lean 4 (by pattern match)
- Safe to use; C++ semantics preserved

**C++ `return a / b;` (with UB on b=0)** → Lean with guard
```lean
def divide (a b : UInt64) (hb : b ≠ 0) : UInt64 := a / b
-- Or use Option: def divide? (a b : UInt64) : Option UInt64 := ...
```

---

## 11. Common Pitfalls & Solutions

| Pitfall | Lean Symptom | Solution |
|---------|--------------|----------|
| `.get!` expected to panic | Returns `default` instead | Use bounded `.get i h` in proofs; `.get!` is unsafe |
| Division by zero is valid | Lean computes `a / 0 = 0` | Add explicit guard `(h : b ≠ 0)` to function signature |
| Array out-of-bounds no-op | `.set!` silently ignored on oob | Check bounds before calling; use `.set i v h` with proof |
| Shift by ≥ 64 returns 0 | Expectation: maybe UB? | Lean defines shifts as wrapping; document when translating |
| `Array.range n` includes 0 | Off-by-one errors | `Array.range n = #[0, 1, ..., n-1]`; check C++ `std::iota` semantics |
| `Bool.xor` doesn't short-circuit | Performance surprise | Use custom `if-then-else` for conditional xor if needed |
| Comparing `Nat.Prime` | Requires proof instance | Use `Fact (Nat.Prime p)` as context var; never case-split on primality directly |

---

## 12. References & Further Reading

- **Lean 4 Standard Library**: https://github.com/leanprover/std4
- **Mathlib4**: https://github.com/leanprover-community/mathlib4 (search `Nat.gcd`, `Nat.Prime`)
- **Batteries**: https://github.com/leanprover-community/batteries (Batteries.Data.UInt, Array extensions)
- **Lean 4 Manual**: https://leanprover.github.io/lean4/doc/ (type system, tactics, standard ops)
- **CLPoly CLAUDE.md**: Existing Lean code structure & conventions (see `proof/CLAUDE.md`)

---

**Document Version**: 1.0  
**Last Updated**: 2026-04-19  
**Maintainer**: CLPoly Verification Team  
**Target**: cpp2lean v2 translator reference

