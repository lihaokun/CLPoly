import Lake
open Lake DSL

package CLPoly where
  leanOptions := #[⟨`autoImplicit, false⟩]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "stable"

@[default_target]
lean_lib CLPoly where
  srcDir := "."

-- B2B (back-to-back) 测试驱动：stdin/stdout NDJSON，由 proof/b2b/runner/runner.py
-- spawn 后逐条调用。
lean_lib B2B where
  srcDir := "."

-- B2B driver 走解释器（native 编译 346-def mutual block 的 Corpus 在启动时
-- 卡死）。runner.py 用 `lake env lean --run B2B/Driver.lean`。
-- TestTypes 不依赖 Corpus，可走 native。
lean_exe b2b_test_types where
  root := `B2B.TestTypes
  supportInterpreter := true
