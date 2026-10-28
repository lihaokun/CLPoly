# B2B (back-to-back) 语义测试

验证 Lean 翻译模型与 C++ 实现端到端语义一致。

## 目录布局

```
proof/b2b/
  Makefile                 # cpp_driver 构建
  driver/
    cpp_driver.cc          # C++ stdin/stdout NDJSON binary
  types/                   # 类型 parse/serialize 库（Step 1+）
  registry/                # 函数 dispatch 表（Step 2+）
  vectors/                 # 测试向量（Step 4+）
  runner/
    runner.py              # 编排：spawn 两端 driver，diff 输出
    diff.py                # 单条结果对比
  third_party/nlohmann/    # 单 header JSON 库（vendored from miniconda）
```

Lean driver 在主 lake 项目（避免 mathlib 重复 fetch）：
```
proof/lean/B2B/Driver.lean              # Lean 主入口（走解释器）
proof/lean/B2B/Registry.lean            # 函数 dispatch
proof/lean/B2B/Types.lean               # 类型库
proof/lean/B2B/TestTypes.lean           # 类型库 round-trip 测试
```

**Lean 端走解释器**（`lake env lean --run B2B/Driver.lean`）：native 编译 346-def
mutual block 的 Corpus.lean 在启动时卡死（运行时初始化巨型 mutual closure 表
导致 sigsuspend 死循环）。解释器模式启动 ~2s，之后每条请求 ~ms，对 B2B 批量测
试可接受。

## 协议

NDJSON（每行一条 JSON）：

请求：
```
{"id":"<任意串>", "fn":"<函数名>", "args":[<参数列表>]}
```

响应：
```
{"id":"...", "ok":true, "ret":{"type":"...","val":...}}
{"id":"...", "ok":false, "err":"<错误描述>"}
```

参数/返回值的 `{"type","val"}` 表示见 `types/`。

## 构建

```bash
# C++ 端
cd proof/b2b && make cpp_driver
./cpp_driver < some_request.ndjson

# Lean 端
cd proof/lean
lake build b2b_driver
./.lake/build/bin/b2b_driver < some_request.ndjson
```

## 运行（Step 3+）

```bash
python3 proof/b2b/runner/runner.py proof/b2b/vectors/__make_zp.json
```
