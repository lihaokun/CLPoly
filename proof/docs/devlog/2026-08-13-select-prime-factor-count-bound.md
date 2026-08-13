# Select-prime factor-count bound

The strict `tryCandidateRaw` refinement now proves that every successful
factor array fits strictly below `UInt64.size`.  This is derived from the
actual execution result, rather than assumed as a callback contract:

- strict DDF and EDF establish that every decoded factor is irreducible and
  monic;
- hence every factor has positive polynomial degree;
- the decoded product is associated to the monic modular input, and monicity
  turns this association into equality;
- additivity of degree over the factor list bounds its length by the input
  degree;
- the existing C++ degree precondition (`< 2^62`) implies the machine-size
  bound (`< 2^64 = UInt64.size`).

This closes the fact needed to prove that the first accepted candidate really
updates C++ `best`, whose initial count is `UINT64_MAX`.  No fuel, partial
definition, semantic oracle, or abstract factor witness was introduced.
