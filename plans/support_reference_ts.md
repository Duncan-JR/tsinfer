# Support reference_ts for split matching

## Objective

Allow the evaluation pipeline to match samples against an ancestor tree sequence
already produced by this tsinfer 1.0 implementation. The result must equal an
uninterrupted ancestor-and-sample match, apart from provenance records.

This is a local pipeline patch. Implement the existing `match` reference argument
and configured reference path for this specific handoff. No C changes are needed.

## Supported input contract

The caller produces `ancestors_ts` using the existing sequence-length workaround:

1. Configure both the ancestor and sample sources in `match.sources`.
2. Compute the jobs with `pipeline.compute_groups_json(cfg)`.
3. Find the first sample group. All ancestor jobs must precede it.
4. Call `pipeline.match(cfg, group_stop=first_sample_group)`.
5. Keep the returned raw tree sequence, or dump it for the sample-matching call.

Because the sample source is configured in step 4, `match` obtains the full contig
length from that source. It also creates the sample individual and population
rows when it builds the initial tree sequence, before matching any samples.

For the MVP:

- Use one ancestor source and one sample source, with contemporary samples and
  strictly older ancestors. Ancestors have `node_flags=0` and
  `create_individuals=False`.
- Use the same sample VCZ, source name, sample selection, order, ploidy, individual
  settings, and population settings for the second call.
- The reference retains its sites, allele meanings, sequence length, metadata,
  individual and population tables, and tsinfer root structure.
- Use biallelic inputs supported by the current matcher. The reference must
  already satisfy the current C matcher's mutation constraints.
- Set `path_compression=False`. Work with raw TSs before post-processing or site
  augmentation.
- Run sample matching without a work directory. If the ancestor stage uses one,
  do not reuse it for the second call: its job manifest describes the full run.

These are documented caller assumptions, not a new validation framework. This
patch establishes the unmodified ancestor-TS handoff. Later experiments that alter
the reference are responsible for preserving these assumptions.

## Existing code to change

- `tsinfer/pipeline.py`: reference resolution, job bounds, and starting TS selection
  in `match`; pass reference positions through the existing job-building calls.
- `tsinfer/config.py`: defer the ancestor-or-reference requirement until `match`
  can see its Python argument.
- `tests/test_pipeline.py`: one new equivalence test family.
- `tests/test_config.py`: remove the obsolete constructor/parser rejection tests
  affected by deferring that requirement. Add no new configuration test family.

Keep `matching.py`, `vcz.py`, the C extension, and `run()` unchanged. Reuse the
existing readers, grouping, `Matcher`, `extend_ts`, and match-file writing.

## Implementation steps

### 1. Resolve the reference before building match jobs

Add a small private resolver in `pipeline.py`, shared by `match` and
`compute_groups_json`:

- Return the explicit Python `reference_ts` when it is not `None`.
- Otherwise load `cfg.match.reference_ts` with `tskit.load` when it is not `None`.
- Otherwise return `None`.

The Python argument remains a `tskit.TreeSequence | None`; paths belong in
`MatchConfig.reference_ts`. The explicit argument takes precedence over the
configured path. Do not require writing an in-memory reference to disk.

`Config.__post_init__` currently rejects `ancestors=[]` unless the configured
reference path is set. It cannot know about `match(cfg, reference_ts=...)`.
Move only this existing requirement into `match`, after reference resolution:
there must be either an ancestor configuration or a resolved reference. Retain
the other configuration checks. Remove the existing config tests whose sole
expectation is rejection of this now-constructible configuration.

Keep `ancestral_state` as an existing required config field. The caller can retain
the original sample annotation configuration; making it optional throughout the
config system is unnecessary for this patch. Reference-mode matching gets its
ancestral states from the TS and does not read that annotation store.

### 2. Give reference mode its positions and ancestral alleles

In `match`, branch on the resolved reference before opening an ancestor store.

For a supplied reference:

- Read its site positions into the existing `int32` positions array.
- Collect `site.ancestral_state` in site order for the reader.
- Use `reference_ts.sequence_length` for the match sequence length and logging.
- Preserve `sequence_intervals` through the reference's existing metadata. There
  is no need to reconstruct it or use it to construct a new root.

For the ordinary run, retain the current ancestor-VCZ metadata loading and
sequence-length selection. Read `sample_time` and compute the artificial root
times only on this ordinary root-building path.

No expression in reference mode should dereference `cfg.ancestors[0]`. The sample
configuration must work with `ancestors=[]` and no ancestor source in `sources`.

### 3. Pass reference positions into job construction

Add an optional `positions` argument to `_build_jobs_with_metadata`. When supplied,
use it for the existing default job bounds. When `None`, keep the current lookup
from the ancestor store at this leaf function.

Pass the positions already obtained by `match` into the job builder. Allow
`_setup_workdir` to forward this argument to its existing fresh-job construction
call, without altering its checkpoint logic. Existing callers can pass `None`.

Have `compute_groups_json` use the resolver and pass reference positions when a
configured reference exists. Its ordinary ancestor-and-sample use remains the
way to calculate the first-stage stopping point.

Generate jobs only for `cfg.match.sources`. For the second call this contains only
the original sample source, so ancestors are not matched again. Sample job indices
and group numbers can restart from zero; they are not reference node IDs. Preserve
the sample ordering so appended node order agrees with the uninterrupted run.

### 4. Start from the reference and reuse its metadata tables

At the existing starting-TS selection, keep checkpoint precedence, then use the
resolved reference if present, otherwise call `matching.make_root_ts` as before.
The supported second-stage call has no workdir, so it takes the reference branch.

Do not add roots, reorder nodes, clear flags, clear mutation times, simplify, or
change sequence length. `extend_ts` already dumps tables before adding new nodes,
so the input TS remains unchanged.

Do not append the individual/population metadata rows returned by the new
job-building call to the reference. They already exist from the first stage.
With the same sample source and selection, and no ancestor individuals, the
sample-only job builder assigns the same individual and population IDs as the
full job builder. Use those IDs directly. Do not offset them by the existing
table lengths: that would duplicate rows or reference nonexistent rows.

This is the reason to preserve the sample configuration and metadata tables in
the input contract. Supporting new cohorts or arbitrary existing individual
tables would require a different policy and is outside this patch.

### 5. Keep the existing allele reader and group loop

Construct `HaplotypeReader` from the sample source, reference positions, and
reference ancestral states. In the supported biallelic case, the same sample VCZ
supplies the same allele vocabulary as in the full run. The ancestral state stays
code zero and the sole derived allele stays code one.

No new allele-merging code or `AlleleMapper` interface is required. References
containing novel alleles or requiring multiallelic reconciliation are outside the
contract. Existing C restrictions on repeated mutations are unchanged.

Run `_process_group` and `extend_ts` normally. Keep the existing JSONL match output;
its parent IDs refer to the reference TS node IDs. Use separate match-file paths
for the ancestor and sample calls because each invocation opens its file in write
mode. Group and haplotype indices in those files are local to the invocation.

Document the supported contract and argument precedence in the `match` docstring.
The existing CLI already calls `pipeline.match`, so loading the configured path
there enables its reference mode without a new CLI command or flag. Do not extend
the full `run()` pipeline to bypass ancestor generation in this patch.

## Only new test family: split matching equals one run

Add `TestReferenceTsEquivalence` in `tests/test_pipeline.py`, using pytest fixtures
and the existing VCZ helpers. Every new test must compare the two execution paths;
do not add arbitrary-reference, invalid-input, or reference-transformation suites.

Use a small deterministic, recombining fixture with several ancestors, positive
ancestor times, and contemporary samples. Ensure its last inference position is
well before its contig length, so a truncated sequence length is observable.
Use compatible biallelic data and include ordinary sample metadata in the fixture.

For each case:

1. Create the sample and ancestor VCZs once, and construct the full config with
   both sources and path compression disabled.
2. Run `pipeline.match(full_cfg)` to obtain the expected raw TS.
3. Compute the full jobs. Find the first sample group and assert directly that
   all ancestor jobs are before it and all sample jobs are at or after it.
4. With a separate first-stage invocation, stop at that group and obtain the raw
   ancestor TS. Assert it has the sample contig length and no matched sample nodes.
   Keep the pre-created individual/population rows.
5. Construct a fresh sample-only config with `ancestors=[]`, unchanged sample
   settings, and no workdir. Exercise the in-memory argument or configured path.
   Do not pass the old ancestor store or the first-stage checkpoint manifest.
6. Match samples and compare the complete tables with the uninterrupted result,
   ignoring provenance only. Do not simplify, reorder, or discard metadata to
   make the comparison pass. Counts alone are insufficient.

Use a small parameterization of this same family for haploid/diploid samples and
an in-memory reference/configured reference path. The path case should dump and
reload the ancestor TS. At least one equivalence case should have configured
individual metadata and populations so incorrect row duplication is caught.
Use fixed parameters and one matching thread to make failures easy to interpret.

The provenance will differ because the configurations and group numbering differ.
Node, edge, site, mutation, individual and population tables, table schemas,
top-level metadata, and sequence length must agree. Also assert the reference TS
was not changed by sample matching.

## Scope boundaries

Do not implement a general reference normalization layer, root insertion, age
filtering, node-ID remapping, metadata merging, input coercion, or fallback repair.
Do not alter the HMM, recurrent-mutation support, or sequence-length semantics of
ordinary ancestor-only runs. Do not implement ancestor restriction experiments in
this change. Do not introduce standalone reader/config validation test families.

Follow `CLAUDE.md`: module imports at file scope, explicit `None` comparisons,
simple operations with named intermediate values, and 89-character lines. Avoid
new abstractions when the existing match flow needs only a branch and positions
passed to its helpers. Do not refactor existing tuple returns as unrelated cleanup;
use a dataclass if a genuinely new multi-value result becomes necessary.

## Verification when implementing

Run the new equivalence family first, then the existing pipeline/config/CLI tests:

```bash
uv run pytest tests/test_pipeline.py::TestReferenceTsEquivalence -v
uv run pytest tests/test_pipeline.py tests/test_config.py tests/test_cli.py -v
uv run ruff check tsinfer/pipeline.py tsinfer/config.py \
    tests/test_pipeline.py tests/test_config.py
uv run ruff format --check tsinfer/pipeline.py tsinfer/config.py \
    tests/test_pipeline.py tests/test_config.py
```

Before committing an implementation, run the test suite with coverage as required
by `CLAUDE.md`. No C rebuild is necessary unless the implementation expands beyond
this plan and changes C code.

Completion means the reference route can run with no ancestor VCZ configured, and
all equivalence cases reproduce the uninterrupted raw TS apart from provenance.
