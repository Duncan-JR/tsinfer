# Efficient stitched matching MVP

## Objective and scope

Implement reduced-HMM and exact-greedy stitched matching without reconstructing
an ancestor tree from the beginning for every focal or gap region. Keep
orchestration in `tsinfer-match-eval`: it owns haplotype tasks, focal/gap
partitions, budget loops, reference selection, and assembly. Tsinfer supplies a
small opt-in mechanism for retaining and advancing its existing native topology.

Base tsinfer work on the local `fix_reference_ts` branch. Its bounded observation
loops, endpoint filtering, regional output copy, and index cleanup already
exist. Preserve those fixes. Follow each repository's
`AGENTS.md`: straightforward functions, explicit intermediate variables,
dataclasses for structured Python results, module imports at the top, explicit
`None` checks, pathlib, Zarr v3, and `uv run` for Python tooling. Put new workflow
parameters in `config.yaml`. Do not add `Co-authored-By` commit trailers.

Include three supporting optimizations with clear, localized benefits:
region-sized native path buffers, run-based stitched paths, and bounded matching
results with one final table assembly. These remove unnecessary allocation,
dense path materialization, and cohort-wide retention around the new procedure.

Keep the existing per-haplotype simplification and scientific method. Shared
reduced references, reference pruning, small-focal-island removal, multiprocessing,
sparse node resets, reader-cache redesign, batched candidate archives, compact
native traceback indexing, and allocator tuning are outside this MVP. Ordinary
full-HMM worker setup is outside the work too. There is no threshold experiment,
new confidence rule, or change to the focal/gap partition.

## Procedure and invariants

A worker processes one whole haplotype at a time. Within each mismatch budget,
process its nonempty regions from left to right, preserving inference-interval
splits. Parallelism remains across haplotypes; regional calls are sequential.

For reduced-HMM matching, keep two independent native matchers:

| Matcher | Work when it is next needed |
| --- | --- |
| Full-reference matcher | Advance its own topology to the next gap start, then run an independent HMM through that gap |
| Haplotype-specific reduced matcher | Advance its own topology to the next focal start, then run an independent HMM through that focal region |

For example, on `gap 1 -> focal 1 -> gap 2 -> focal 2`, the full matcher matches
gap 1, later advances through focal 1, and matches gap 2. The reduced matcher
advances through the leading gap, matches focal 1, later advances through gap 2,
and matches focal 2. Simplification changes node IDs and edges, so these matchers
cannot exchange topology. Map reduced copying parents back to original IDs only
when collecting results.

For greedy matching, construct exact focal copying runs from the existing
intervals. Match gaps with the full matcher in genomic order. Before a later
gap, that matcher advances topology through intervening focal spans without
running an HMM there. No reduced matcher or dummy focal HMM call is needed.

The following invariants define the implementation:

- Parent, child, and sibling arrays retained by a matcher describe the same tree,
  at a documented position on that matcher's own reference.
- Its insertion/removal cursors identify the next unprocessed forward events.
- Each regional HMM starts with fresh likelihoods and traceback decisions.
  Retained topology does not condition one region on another.
- Forward matching stops at `positions[end - 1]`, the final processed site, not
  at the output right endpoint `positions[end]`. Missing calls remain neutral.
- Traceback changes temporary parents, leaving the retained forward topology
  intact. The same regional path, copied alleles, and mutations must result.
- A new haplotype or budget begins a fresh traversal. The original reference
  remains the copying reference throughout sample assembly.

This removes repeated prefix/suffix reconstruction. Edge changes still have to
be applied through skipped spans, and each matched region still needs forward
HMM work and reverse traceback. Per-region dense probability resets and one
node-sized parent copy remain. Do not describe the result as zero tree traversal
or runtime proportional only to the number of gap sites.

## Execution and commit sequence

The repository sections below are independently executable. Complete the
**tsinfer sequence T1 -> T2 -> T3** first, rebuilding and verifying its native
extension. Then execute the **match-eval sequence M1 -> M2 -> M3** against that
checkout. Make one commit per step, including its essential checks. Do not squash
the changes into one commit or introduce additional optimization phases.

| Order | Repository / step | Commit message |
| --- | --- | --- |
| 1 | tsinfer / T1 | `Retain native topology across independent regional matches` |
| 2 | tsinfer / T2 | `Size regional path buffers to the matched interval` |
| 3 | tsinfer / T3 | `Separate match result appending from table finalisation` |
| 4 | tsinfer-match-eval / M1 | `Reuse reference traversal for stitched focal and gap matching` |
| 5 | tsinfer-match-eval / M2 | `Assemble stitched paths directly from copying runs` |
| 6 | tsinfer-match-eval / M3 | `Bound stitched matching results with haplotype batches` |

Before implementing either sequence, save small pre-change stitched results for
both methods, and HMM results for two budgets, outside committed source. Use
these as exact integration references. Do not create a benchmark framework or a
second production matching implementation to preserve the old behavior.

## Changes in tsinfer

### T1. Retain topology and protect it from traceback

**Commit:** `Retain native topology across independent regional matches`

**Files:** `lib/ancestor_matcher.c`, `lib/tsinfer.h`, `_tsinfermodule.c`,
`tsinfer/matching.py`, and the existing C/Python matching tests.

#### Interface and lifetime

Add a keyword-only `reuse_topology=False` option to the Python
`matching.AncestorMatcher` wrapper and expose the corresponding native option
through the extension constructor and a C matcher flag. Match-eval explicitly
uses `True`. Keep the existing `find_path(h, start, end, match_out)` wrapper
signature and low-level return tuple. The ordinary full-HMM path uses the
existing default; do not change `Matcher._match_one` or its worker setup.

Expose `reset_topology()` on the matcher, backed by a small native reset function.
It clears the retained topology and cursor validity without reallocating the
matcher, indexes, probability arrays, or scratch. Per-region probability reset
still happens inside `find_path`; this method is not a general session API.

Add only scalar traversal bookkeeping to `ancestor_matcher_t`: current genomic
position, forward insertion/removal cursors, and a validity flag. Use signed
cursor types compatible with the existing loops. Reuse the existing parent,
child, and sibling arrays. No new tree object, reference index, per-site snapshot,
or state transfer between matchers is needed.

#### Separate topology from probability reset

Split the current reset responsibilities into small functions:

1. Topology reset clears parent/child/sibling arrays, resets cursors, and marks
   the retained position invalid. Initialize on first use and explicit reset.
2. Per-match reset clears allelic and recombination working state, resets the
   traceback allocator and likelihood-node count, and retains the existing
   allocator-size policy. It does not clear a valid retained topology.
3. Initialize likelihood/cache arrays from the actual starting topology with the
   current dense loop, including isolated-node and root-zero behavior. Keep
   recombination-state reset out of the topology-only reset function.

For default matchers, retain the current independent-call initialization and
traceback setup. For reuse-enabled matchers, use the following traversal path.
Keep shared likelihood and traceback loops; do not duplicate the whole engine.

#### Establish the starting topology

If topology is invalid, clear it and use the existing filtered starting-tree
scan once: install only edges active at `positions[start]`, and skip already
expired removal events. Save both cursors and the position. A request earlier
than the retained position uses this same cold initialization; this simple
fallback preserves ordinary arbitrary-order use without a random-access index.

Otherwise, advance from the retained cursors to `positions[start]`. At each
intervening breakpoint, apply removals before insertions, in the existing index
order, using `remove_edge` and `insert_edge`. Process events at the target
position to obtain half-open edge semantics. Do not call site likelihood
updates, decode observations, store traceback, or adjust old likelihoods during
this advance. Fresh likelihood initialization follows it.

Recompute the next breakpoint from the remaining insertion/removal events and
sequence length. Preserve existing exhausted-index handling, including references
with no edges. Keep this topology-only advance private to the native engine;
match-eval triggers it by requesting the next region, rather than calling a
public seek operation for every focal span.

The forward routine then uses its current likelihood logic and edge transitions.
On success, retain the position and cursors at the final processed site. Do not
advance to `positions[end]` just because output edges reach that coordinate.

#### Trace back without rewinding retained parents

The current traceback clears parents and scans edge-array suffixes to rebuild
the endpoint tree, then rewinds those parents to the region start. Child/sibling
arrays remain at the forward endpoint. Replace that setup in reuse mode:

1. Reuse the node-sized `likelihood_nodes_tmp` buffer as parent scratch after
   forward likelihood coalescing has finished. Rename it to `node_scratch` and
   document its two non-overlapping lifetimes: forward coalescing and reverse
   parent links. Verify that no subsequent forward operation needs its contents.
2. Copy the retained endpoint parents into scratch once. Pass an explicit parent
   pointer to traceback rather than swapping `self->parent` or copying it back.
3. Start reverse insertion at `forward_removal - 1` and reverse removal at
   `forward_insertion - 1`, because traceback uses the opposite event arrays.
   These are local reverse cursors; they must not overwrite retained cursors.
4. Remove endpoint parent clearing and suffix reconstruction from this reuse
   path. Preserve the existing recombination scratch reset, reverse event
   ordering, copied-allele reconstruction, switch decisions, and output order.
5. Keep the output right coordinate at `positions[end]`, while starting from the
   tree at `positions[end - 1]`. Handed-off cursors exclude later unobserved
   topology changes. Only scratch parents travel backwards through the region.

After return, all retained topology arrays still describe the forward endpoint.
Scratch contents are disposable; the next coalescing pass overwrites what it
uses. There is no additional node-sized allocation and no traversal to restore
the endpoint after traceback. The parent copy is still O(number of nodes) per
region and must be included in the performance assessment.

#### Boundary cleanup and failure behavior

Keep traceback record reuse within a single region only. Pass the region start
through the existing site-update/store routines and compare with the previous
record only when `site > start`, not merely `site > 0`. This prevents reading
records from another call after allocator reset.

Clear traceback headers and maximum-node entries for the requested slice on
both success and failure, before a later call can reuse allocator memory. On a
native matching error, invalidate retained topology and propagate the existing
error. The next call performs ordinary cold initialization and probability
reset. Do not add retries, checkpoint recovery, or retained partial likelihoods.

Validate `0 <= start < end <= num_sites` at the extension boundary before native
access. This is the necessary bound for the interface, not a reason to add
repeated shape/range checks inside loops. Keep existing haplotype-length checks.

#### Essential tests

Extend existing tests rather than creating a new suite:

- Add one C sequence case with reuse enabled: ascending separated regions,
  a region boundary at a breakpoint, and explicit reset followed by another
  match. Compare paths and matched alleles with default independent matchers;
  check that retained parent/child/sibling arrays describe the endpoint tree and
  cursors advance monotonically within the pass. Exercise earlier-request cold
  initialization and failure recovery in this case or the existing error tests.
- Reuse `test_matching_regions` fixtures for single sites, unobserved endpoint
  trees, exhausted indexes, and empty edge arrays; the zero-edge fixture already
  exists. Keep default-mode assertions where parents finish at the start, and
  reuse-mode assertions where retained parents finish at the endpoint.
- Parameterise Python `TestRegionalMatching` for the new mode, retaining its
  exact expected paths, missing/outside observation checks, and impossible-match
  recovery. Add one small explicit-reset check and parameterised invalid bounds.
  These cover the binding and wrapper without another combinatorial matrix.

Inspect C coverage for new branches. Do not accept changed expected paths merely
because fresh and reused matchers agree with each other; existing fixed expected
paths remain an independent correctness check.

### T2. Allocate native paths for the region

**Commit:** `Size regional path buffers to the matched interval`

**File:** `_tsinfermodule.c`, with an extension-facing test in the existing tests.

After bounds validation, allocate `left`, `right`, and `parent` with capacity
`end - start`. At most one copying segment can start at each processed site.
Preserve returned path length, dtypes, reverse order, and absolute BP coordinates.
Slices returned by the Python wrapper otherwise retain chromosome-sized backing
arrays, so this also reduces retained path memory.

Keep the native matched-allele array full-length and keep the wrapper's
`match_out[start:end] = matched[start:end]`. Do not rebase native site indices,
compact `traceback`/`max_likelihood_node`, or change their accessors. Native output
arrays are uninitialised, so allocation capacity must not be reported as a
whole-array write on each call.

**Essential test:** one parameterised low-level test for single-site and
full-span capacities and valid path entries. Bounds and regional correctness
are already checked in T1. Rerun those checks after rebuilding; no new C-library
algorithm is introduced by this step.

### T3. Provide simple append and finalise helpers

**Commit:** `Separate match result appending from table finalisation`

**File:** `tsinfer/matching.py` and existing assembly tests. Required by M3.

Extract two module-level helpers from `extend_ts`:

- `append_match_results(tables, paired_results, allele_mapper)` sorts the supplied
  batch by `haplotype_index`, appends the existing node metadata, flags,
  individual/population assignments and edges, decodes that batch's mutation
  alleles, and appends mutations. Release its temporary mutation data on return.
- `finalise_match_tables(tables, provenance_record=None)` performs the existing
  sort, index construction, mutation-parent computation, provenance addition,
  and tree-sequence creation, in the same order.

`extend_ts` dumps tables once, calls these helpers for its one group, and returns
its final TS. Preserve its signature and the pipeline's group-by-group reference
updates. Keep metadata and allele decoding logic in one place. Derive positions
from the tables when appending; do not introduce a table writer class or another
cache object to avoid a modest per-batch immutable-array cost.

The helpers retain no match results between calls. Global batch order is a
caller contract: sorting one batch cannot repair batches appended out of order.
The evaluation caller will keep matching against its original immutable
reference, while tables receive samples.

**Essential test:** one `TestExtendTs` case comparing two canonically ordered
batches with a single `extend_ts` call on the same results. Include a mutation
and existing metadata/individual/population assignments. Compare final tables
excluding provenance. Existing assembly tests cover the rest.

## Changes in tsinfer-match-eval

### M1. Use persistent traversal for both stitched methods

**Commit:** `Reuse reference traversal for stitched focal and gap matching`

**Files:** `lib/matching.py` and its execution notes in `README.md`.
**Prerequisites:** T1 and T2.

Keep the existing division into `StitchedInputs`, `StitchedWorker`, and a
whole-haplotype task. Reference indexes and probability arrays already belong to
shared prepared inputs; do not rebuild them per region. Each thread owns its
mutable full matcher. Construct it with `reuse_topology=True`, and add one
reusable full-site `match_out` array to its worker dataclass.

For each haplotype:

1. Read the haplotype once, preserve current called-domain/all-missing behavior,
   and inspect existing mask slices for nonempty focal coverage.
2. For method `hmm`, construct the current haplotype-specific reduced TS, original
   node mapping, indexes, and reuse-enabled matcher once if any requested budget
   has coverage. Keep current simplification flags and retained-node selection.
   For method `greedy`, prepare no reduced TS or matcher.
3. At the start of every budget, reset topology on the full matcher and on the
   reduced matcher when present. This also handles the first budget of a new
   haplotype. Do not reallocate either matcher for these resets.
4. Walk existing mask runs in order, clip them to the called domain, skip empty
   runs, and retain existing inference-interval splits. For a gap, call the full
   matcher. For an HMM focal region, call the reduced matcher. The native call
   itself advances the selected reference from its retained position to start.
5. For greedy focal regions, use the existing exact interval cover and make no
   HMM call. The full matcher will advance through skipped focal spans when the
   next gap is requested. No advance is needed after its final gap.
6. Read only `match_out[start:end]` after each call. Preserve the current reduced
   parent mapping and merge behavior until M2 changes their representation.

Keep one `ThreadPoolExecutor` for whole-haplotype tasks and existing background
reader threads. Never submit focal/gap regions separately, nest a second pool,
share mutable matcher state, or move edge updates into a Python per-edge loop.
The full matcher and its scratch survive successive tasks on the same thread;
the reduced matcher is released after its haplotype. Simplification still holds
the GIL; report that remaining limitation instead of adding a process backend.

Update the README to explain the two independent retained traversals, fresh HMM
starts, and topology-only fast-forwarding. Keep method labels and output formats.
M1 can use the existing dense assembly briefly so that traversal correctness is
reviewable separately from M2.

**Essential verification:** compare small saved pre-change results exactly for
`hmm` and budget-zero `greedy`, using one and two matching workers. For HMM, also
check two budgets and consecutive haplotypes on one worker to exercise resets.
Include leading/trailing gaps, alternating focal/gap runs, and empty focal
coverage. Inspect regional call counts against masks after inference-boundary
splitting. Use manual checks; do not add a match-eval test framework.

### M2. Keep stitched copying paths as runs

**Commit:** `Assemble stitched paths directly from copying runs`

**File:** `lib/matching.py`. **Prerequisite:** M1.

The current wrapper expands HMM paths into per-site parents, copies regional
alleles, fills full-site parent/copied arrays, and compresses them again. Remove
that round trip while preserving the matching intervals and scientific outputs.

Add a small `CopyingRuns` dataclass with numeric `left_site`, `right_site`, and
`parent` arrays. Runs use half-open bounds on the existing inference-site axis.
Use a region-result dataclass for runs and owned mutation results; do not add a
class hierarchy or support both dense and run representations indefinitely.

In `_match_region`:

1. Run the existing regional call using worker scratch. Reverse native path
   arrays into genomic order and convert their BP boundaries with vectorised
   `searchsorted`. Clip site bounds to the exact region requested.
2. Map reduced parents through `original_nodes` once here when a mapping is
   supplied; full parents already use original IDs. Keep the operation explicit.
3. Compare called `haplotype[start:end]` against `match_out[start:end]` immediately.
   Add `start` to local mismatch indices and construct owned native `Mutation`
   objects. No copied-allele array survives into the stitched result.
4. Return compact numeric runs and mutations. No result may reference reusable
   matched-allele scratch. Do not materialize one parent per site.

Change the Numba greedy sweep to emit one numeric run per chosen interval instead
of filling dense parents. Preserve the stable left-bound sorting, maximum-reach
choice, and equal-reach tie behavior. Allocate output capacity from the number
of eligible input intervals in its Python wrapper, pass arrays into the compiled
sweep, and return the number written. Wrap the valid output in `CopyingRuns`
outside Numba. Each selected interval consumes a new input choice, so this
capacity suffices; no chromosome-sized allocation or resizing machinery is
needed. Return compact owned arrays rather than retaining large unused capacity.

Compute greedy focal runs once for the haplotype's budget-zero mask. Use a running
cursor through those ordered runs while walking focal/gap components; append
focal results and HMM gap results in genomic order. Exact greedy focal runs add
no mismatches under the existing validated interval contract.

Replace `merge_paths` with collection of ordered `CopyingRuns` blocks and a final
conversion step. Concatenate each numeric array once after the regions finish,
then coalesce in one linear pass; do not concatenate growing arrays at every
region. Merge adjacent runs only when site bounds touch and parent IDs agree.
Such output merging must never remove an independent matching boundary. Keep
inference-boundary splits for HMM calls, while preserving existing output seam
merging. Convert final runs to ordinary native `PathSegment` objects once:
site-zero left endpoint maps to BP zero; the final site-axis endpoint maps to
sequence length; other endpoints use inference positions. Collect mutations in
site order. Remove obsolete dense parent/copied arrays and their helpers.

**Essential verification:** rerun the saved comparisons from M1. Add a small
manual seam example with two touching equal-parent runs and one changed-parent
run, including outer BP endpoints and a single-site region. Native tests already
cover missing observations and traceback; do not duplicate their matrix here.

### M3. Bound submitted work and assemble results by batch

**Commit:** `Bound stitched matching results with haplotype batches`

**Files:** `lib/matching.py`, `Snakefile`, `config.yaml`, the example configuration
if present, and `README.md`. **Prerequisites:** T3 and M2.

Add one positive integer `haplotype_batch_size` in configuration and pass it from
the stitched rule to `match_stitched_samples`. Validate it once at the library
boundary. Set its configured value to keep the allocated workers occupied when
there are enough haplotypes; do not derive hidden defaults at multiple layers.
This parameter controls matching tasks/results, not candidate generation or
reference construction.

1. Keep the existing canonical job list and reader schedule for the entire job.
   Load reference indexes, masks/intervals, identities, and reader once. Reuse the
   existing cache and read workers; do not replace caches at batch boundaries.
2. Create one result table collection per budget by dumping the original
   reference. Open each budget's existing JSONL output once. Use ordinary helper
   functions and context-managed file handles, not a streaming writer framework.
3. Keep one matching executor and thread-local worker state alive across batches.
   Submit only the next contiguous batch of canonical jobs. With one worker,
   use the same batch loop sequentially. There is no future for every cohort
   haplotype and no cohort-wide `results` list.
4. Collect the completed batch and sort by `haplotype_index`. For each budget,
   append its results through `tsinfer.matching.append_match_results` and write
   JSONL records using the existing schema. Refactor `_write_stitched_matches`
   to accept an open file rather than reopening/overwriting a path per batch.
5. Release batch futures and result objects before admitting the next batch.
   Continue matching against the original reference, never growing result tables.
6. After the final batch, finalise each budget's tables once, dump its tree
   sequence, and release those final tables/TS before finalising the next budget.
   Do not invoke `extend_ts` on a growing TS after each batch.

Keep sample-node order, source/sample/ploidy identities, metadata, allele states,
and original parent IDs identical. Expected working result memory scales with
batch size. Final tables necessarily scale with cohort output and budget count;
input masks/associations and the existing reader's retained chunks are still
cohort/reference-sized. State those limits explicitly rather than claiming
complete cohort-bounded preparation memory. Batched candidate archives and
site-streaming are not part of this change.

**Essential verification:** compare one batch with several small batches on the
same saved workload, using one and two workers. Compare final scientific tables
and canonical JSONL records exactly, ignoring provenance timing only. Confirm
that a second batch reuses workers and full matchers and that memory no longer
retains every haplotype's Python results. Existing worker failures propagate;
do not add retries or resumable batch state.

## Checks, measurements, and completion

### Required repository checks

For each tsinfer commit, rebuild with `uv sync` after native/extension edits.
Run focused existing matching/assembly tests during development. Before committing,
run `uv run ruff check` and `uv run ruff format --check` on changed Python files,
then `uv run pytest tests/ --cov=tsinfer --cov-report=term-missing`, as required by
`CLAUDE.md`. Inspect coverage for changed Python code. Added tests are limited
to the essential cases listed above.

For T1 native coverage, configure a separate build once:

```sh
meson setup lib/build_stitched_mvp -Db_coverage=true
ninja -C lib/build_stitched_mvp test
cd lib
gcovr build_stitched_mvp -r . --exclude 'subprojects/*' --exclude 'tests/*'
```

Reuse this build and rerun its tests after T2. Inspect coverage of new C-library
logic and use the build's existing strict warning flags. Python coverage does
not establish C coverage. Do not add allocator fault injection or a new memory
checking framework for this proof of concept.

In match-eval, use `uv run ruff check lib` and `uv run ruff format --check lib`
and the listed manual comparisons. Follow its `AGENTS.md`; there is no new
automated test suite. Keep diagnostic scripts/results outside production code.
Force `match_stitched_samples` after library edits: library files are not
Snakemake inputs, so otherwise existing products can conceal the changed code.

### Minimal performance assessment

Use one fixed ancestor reference and haplotype set with many short focal/gap
regions. Preserve masks, AC cutoff, inference intervals, HMM parameters, and
per-haplotype reduced-node selection. Compare complete runs of both stitched
methods with their saved implementations. Check HMM budget reuse separately;
greedy remains restricted to budget zero.

Use the same matching worker count, reader-cache budget, and inputs in baseline
and stitched comparisons. The current workflow requests `matching_threads` for
the ordinary full HMM but all workflow cores for stitched matching; override
these allocations consistently for the experiment. Do not optimize ordinary
full-HMM setup to make this comparison. Use the same tsinfer revision with its
default reuse-disabled path for that baseline.

Record at least three runs and report median elapsed time and peak process RSS.
Include setup, haplotype-specific simplification, matching, and assembly in total
cost; a faster isolated HMM call is not sufficient. Use a short diagnostic to
separate reduction preparation, regional matching, and output work. Do not time
every tiny interval in production or build an instrumentation subsystem.

Confirm from the native control flow, cursor checks, and a short profile that
later ordered calls do not replay earlier edge prefixes or scan reference
suffixes to initiate traceback. There is one cold start per used reference per
haplotype/budget pass, followed by forward cursor advancement. Report the
remaining dense node reset/copy and simplification costs. `total_memory` currently
omits some native arrays, so use measured RSS for memory claims.

Check performance with one worker and a small multi-worker allocation. Keep
parallelism across whole haplotypes; enough cohort tasks normally provide
concurrency without subdividing regions. Do not add another scheduling backend
in response to poor scaling during this MVP.

### Completion criteria

The MVP is complete when both stitched methods preserve their existing paths,
mutations, and called-allele reconstruction; retained topology remains correct
through skipped spans and resets; repeated endpoint reconstruction is removed;
dense stitched paths and cohort-wide Python result retention are removed; and
complete timing/RSS measurements document the effect.

Investigate any exact-output difference before accepting it, including tie
ordering and BP seams. Do not silently update expected results. If per-haplotype
simplification, dense probability resets, or input archives become the remaining
bottleneck, report them as limits of this MVP rather than extending its scope.

## Implementation references

Use the local source at execution time; these are the relevant entry points:

- tsinfer `lib/ancestor_matcher.c`: `ancestor_matcher_reset`, forward matching,
  traceback storage, likelihood coalescing, and `ancestor_matcher_run_traceback`.
- tsinfer `lib/tsinfer.h`: `ancestor_matcher_t` and native matcher flags.
- tsinfer `_tsinfermodule.c`: matcher construction and output allocation.
- tsinfer `tsinfer/matching.py`: wrapper and direct table assembly.
- tsinfer `lib/tests/tests.c`: `test_matching_regions` and matching error cases;
  `tests/test_matching.py`: regional matching and assembly fixtures.
- tsinfer-match-eval `lib/matching.py`: `StitchedInputs`, `StitchedWorker`,
  `_match_region`, greedy cover, whole-haplotype execution, and JSONL output.
- tsinfer-match-eval `Snakefile`, `config.yaml`, and `README.md`: allocation,
  method/budget dispatch, and execution documentation.
