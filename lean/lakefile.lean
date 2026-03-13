import Lake
open Lake DSL

package CLPoly where
  leanOptions := #[⟨`autoImplicit, false⟩]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "stable"

@[default_target]
lean_lib CLPoly where
  srcDir := "."
