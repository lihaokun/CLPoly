# Well-founded uint64 prime enumeration

## C++ correspondence

`next_prime_64` calls GMP `mpz_nextprime` and throws `overflow_error` if its
result exceeds 64 bits.  `prev_prime_64` calls `mpz_prevprime` and throws
`domain_error` when no smaller prime can be requested.  Consequently the
strict L1 interface returns `RawExec UInt64`; a total `UInt64 → UInt64`
boundary would erase observable C++ behavior.

## Natural-language termination and correctness proof

- Upward scanning starts at `p + 1`.  While the candidate is inside the
  uint64 range and not prime, increment it.  The measure
  `UInt64.size - candidate` strictly decreases.  Crossing the boundary returns
  `arithmeticOverflow`.
- Downward scanning starts at `p - 1`.  While the candidate is nonzero and not
  prime, decrement it.  The candidate itself is the decreasing measure.  Zero
  returns `arithmeticDomain`.
- In either successful branch the tested proposition is literally
  `Nat.Prime candidate`; conversion back to `UInt64` is injective because the
  candidate is in range.  Hence every successful result is prime.
- The upward result is strictly greater than its input and the downward result
  strictly smaller.  Therefore the select-prime rank is
  `UInt64.size - p.toNat` in the ascending branch and `p.toNat` in the
  descending branch.
- `PrimeEnumerationTermination.next_decreases` is conditional on an actual
  successful next-prime call.  Error branches return immediately and need no
  recursive decrease.  This exactly matches the generated control flow.

The compatibility functions used by the unverified generated corpus use the
same well-founded scans but fold exceptional exhaustion back to their input.
The verified path never uses that fold; it uses the strict `RawExec` entries.

## 度量

- 耗时：~1.5 小时（C++ 审计、接口修改、证明与编译）
- 迭代：8 轮编译-修复循环
- Lean 新增/修改行数：约 145 行
- 对应 C++ 行数：约 35 行（两个 GMP 包装及 select-prime iterator）
- 放弃的方案：无条件要求有限 `UInt64` 上每次后继都降低 rank；该接口
  无法实例化，因为 C++ 在边界抛异常而不是继续返回机器字
