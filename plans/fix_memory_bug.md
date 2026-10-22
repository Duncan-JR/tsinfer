# Bug report: concurrent cache misses load the same genotype chunk repeatedly

## Summary

`ScheduledCache.get()` treats a one-second wait timeout as permission to load a
missing chunk synchronously, even when another thread is already loading it.
Concurrent matching threads can therefore allocate and decode many copies of the
same genotype chunk. Large Full HMM runs exhibit transient private-memory peaks
above 130 GiB, and the cache summary understates the number of actual loads.

Affected revision: `6beae28f54acd6e657062b32a294c889496066c5`, which is both
the revision identified in the review and the local HEAD inspected for this plan.
The affected code is `ScheduledCache` in `tsinfer/vcz.py`; its caller,
`VCZHaplotypeReader._load_sample_chunk()`, explicitly promises deduplication.

## Reproduction

1. Register a loader for one scheduled `(source_name, chunk_idx)` key. Give the
   cache enough memory for that chunk and enough references for all readers.
2. Start the cache worker and hold its loader open for more than one second.
3. Have several threads call `get()` for the same key while the loader is blocked.
4. Observe that each timed-out caller can invoke the loader again. Release the
   loader and count invocations, rather than just entries in `_chunks`.

A tiny array and an event-controlled loader reproduce the race; multi-gigabyte
VCZ data are not required. A second form occurs when memory pressure makes a
caller initiate the first synchronous load: other callers and a newly unblocked
background worker can race that load as well.

## Expected behavior

For a given key, one thread owns an outstanding load. Other readers wait for its
result or receive its error. A slow load does not trigger another load or an
unnecessary eviction. A genuinely missing, unclaimed chunk can still load
synchronously to break a memory-pressure deadlock. Reloading after a legitimate
eviction remains supported.

## Observed evidence and impact

The review's `cache_audit.csv` records:

| Full HMM workload | Populate records | Stall warnings | Peak private memory |
| --- | ---: | ---: | ---: |
| Out-of-Africa, n=10,000, 10 Mbp, error | 33 | 32 | 143.81 GiB |
| 1KGP chr20p, n=3,202 | 33 | 32 | 131.69 GiB |
| 1KGP chr20p, n=1,012 | 33 | 32 | 79.23 GiB |
| 1KGP chr20p, n=320 | 33 | 32 | 53.25 GiB |

These workloads use one physical sample chunk. The large simulation peaks around
25 seconds; its peak after 60 seconds is 9.44 GiB and its overall median is
8.74 GiB. The review reports that the cache summary nevertheless counts one load.
Duplicate allocations explain the transient peaks, but these measurements do not
establish the corrected peak or wall time. Logged load durations overlap and
cannot be subtracted from elapsed time. Comparisons against the stitching
implementation need a corrected Full HMM baseline.

## Root cause

- `_worker_loop()` advances the schedule and reserves `_pending_bytes`, but does
  not register which key is being loaded.
- `get()` checks only `_chunks` after its timeout, then releases the condition
  lock around the loader without reserving the key or pending bytes.
- Checking `_chunks` after loading only prevents some duplicate insertions; the
  expensive allocations and decoding have already happened. The background
  insertion is unconditional, which can also corrupt retained-byte accounting
  when it overwrites a result inserted by a synchronous caller.
- `_total_loads` and `_total_load_time` are updated only by background workers.
- A background loader exception currently escapes its worker. Simply adding an
  in-flight marker without publishing failures would leave readers waiting for
  a result that can never arrive.

Evidence: section 6 of
[the performance report](/Users/duncan/work/dphil-analysis/data/ch5/mutation_age_performance_review_2026_10_22/mutation_age_performance_report.html)
and its
[cache audit](/Users/duncan/work/dphil-analysis/data/ch5/mutation_age_performance_review_2026_10_22/cache_audit.csv).

---

# Implementation plan

## Scope and design

This document is the only change in the planning task. The subsequent fix should
modify only `tsinfer/vcz.py` and `tests/test_vcz.py`. Keep the public API, chunk
layout, allele mapping, scheduling policy, and C extension unchanged.

Use the existing `threading.Condition` to protect ownership and publication.
Add an in-flight key set and a per-key error dictionary to `ScheduledCache`.
Successful results continue to live only in `_chunks`, so no second container
retains arrays after eviction. A completed error remains available for the
lifetime of this cache; this fix does not introduce automatic retry semantics.
Separate containers are sufficient: no new state class, future, executor, or
general cache abstraction is necessary.

Maintain these invariants under the condition lock:

1. A key has at most one active loader, across workers and synchronous callers.
2. Claiming a key and reserving its estimated bytes happen before releasing the
   lock. `_pending_bytes` includes both kinds of loads and is released once.
3. Completion publishes either the array or the exception, clears ownership,
   and notifies all waiters in the same critical section.
4. `_current_bytes` accounts for arrays retained by `_chunks`. It is not a
   process-memory limit: callers may retain evicted arrays and decoding uses
   temporary buffers.
5. An evicted successful key can be claimed again if still needed. A failed key
   is not silently retried by the next waiter or worker.

## 1. Make both paths claim the key

In `_worker_loop()`, retain the existing budget check and cursor advancement.
When a chunk fits, mark its key in flight alongside the pending-byte reservation
under the lock. Extend `_advance_cursor()` to skip in-flight and failed keys as
well as cached/unneeded keys. This handles a synchronous caller claiming a key
before a worker reaches it, including the wakeup caused by `_force_evict()`.

In `get()`, retain the fast cached-hit path and the existing timeout for deciding
when to attempt fallback. On every loop iteration, check for a published error.
After a timeout, recheck cached/error/in-flight state under the lock:

- If cached, return the existing result.
- If failed, raise the stored loader exception.
- If in flight, continue waiting without evicting or invoking the loader.
- Otherwise, for a registered source, claim the key, reserve its estimate, and
  perform the existing fallback eviction under the lock. Only then release the
  lock and execute the synchronous load.

Keep the existing one-chunk forced-eviction policy for this narrow fix. Do not
turn fallback into another strict budget wait: that would reintroduce the
deadlock covered by `test_memory_pressure_deadlock`. Reserving fallback bytes
ensures background workers see those allocations, but fallback can still exceed
the nominal budget to make progress. Correct the `get()` docstring's claim that
one eviction necessarily keeps loading within budget.

Do not replace the timeout with a longer timeout, serialize I/O under the
condition lock, or remove synchronous recovery. Those changes do not supply
the required ownership invariant.

## 2. Share load completion and failure handling

Extract one small private load helper used by both owners. Its precondition is
that the key and estimated bytes have already been reserved. It runs the loader
outside the condition lock, measures elapsed time, then reacquires the lock to
publish completion and release the reservation. Use explicit intermediate
variables and a short docstring describing the lock/ownership contract.

On success, insert the result once and add `data.nbytes` to `_current_bytes`.
Remove the existing duplicate-after-load check: ownership prevents duplicates
before allocation. On a normal loader exception (`Exception`), publish the
original exception for that key without inserting a result. In both cases,
release pending bytes, remove the in-flight marker, and `notify_all()`.

Have the synchronous caller re-enter its result/error loop after completion.
The background worker may continue to the next scheduled key after publishing
an error; consumers of the failed key raise it. Keep failure handling confined
to loader execution, so bookkeeping defects are not disguised as I/O errors.
Do not add a cache-wide cancellation or shutdown redesign to this patch.

Count each actual loader invocation in `_total_loads`, including synchronous
loads and failed attempts. Accumulate its elapsed time once on completion,
including failures. Do this in the common path under the lock; remove the old
worker-only increments. Retain existing successful load and stall diagnostics,
emitting a stall warning only when fallback actually claims a key. Document
that the summary counts attempts, not just retained results. No additional
metrics framework is needed.

Update the class docstring, which currently says `get()` never triggers a load,
to describe deduplicated background loading and synchronous recovery accurately.

## 3. Minimum regression tests

Use `TestScheduledCache` in `tests/test_vcz.py`, tiny `_arr()` values, and the
existing threading imports. Add only two parameterized test methods, each with
background-owned and synchronous-owned cases. Use a local pytest fixture if
shared setup/cleanup is needed; avoid a reusable concurrency-test framework.

**Concurrent success, two ownership cases.** Hold the first loader with an event
and start three readers for its key. In the background case, wait for the worker
to enter the loader before starting readers. In the synchronous case, first fill
a one-chunk cache with another key, then let one reader reach the existing
fallback and enter the blocked loader before starting the other readers.

Ensure waiting readers cross the one-second timeout while the owner is blocked.
Prefer a test-local wrapper around this cache's condition `wait` method that
signals an event after a real timed wait returns false; preserve its actual
release/reacquire behavior. This gives a concrete synchronization point rather
than assuming a short sleep schedules every thread. Do not change the
production timeout or add a production test hook. Release the loader after the
required timeout observations, then join every reader with bounded waits.

Assert one loader invocation for the target key, identical returned array
objects, correct contents, zero pending bytes, cleared ownership, and retained
bytes equal to the array size. Verify load counts include the synchronous case
(and the initial filler load). In that case, after fallback has claimed the
target, allow the real background worker to wake on the freed space; it must
skip that key. Join the worker before the final invocation-count assertion so
late background duplication cannot escape the test.

**Concurrent failure, two ownership cases.** Reuse the same ownership scenarios
with an event-controlled loader raising a specific `ValueError`. All readers
must finish and report that error; a subsequent `get()` must raise it without
another invocation. Assert no target result is cached, ownership is cleared,
pending bytes return to zero, retained-byte accounting is correct, and the
failed attempt is counted. Capture exceptions explicitly from reader threads
so pytest does not merely emit thread warnings while missing a failure.

Use `finally`/fixture teardown to release the loader event, join reader threads,
and shut down the cache even when an assertion fails. Bound event waits and
assert every joined thread has stopped. Cleanup must not conceal assertions.

Retain the existing memory-pressure deadlock, concurrent waiter, read-ahead,
refcount eviction, byte-tracking, and haplotype-reader tests. They already cover
the surrounding behavior. Do not add RSS thresholds, large simulations, new C
tests, a benchmark harness, or another end-to-end test for this Python race.

## 4. Validation and acceptance

During implementation, demonstrate that the new success regression fails on
the original code due to duplicate calls, then passes with the fix. Run:

```bash
uv run pytest tests/test_vcz.py::TestScheduledCache -v
uv run ruff check tsinfer/vcz.py tests/test_vcz.py
uv run ruff format --check tsinfer/vcz.py tests/test_vcz.py
uv run pytest tests/ --cov=tsinfer --cov-report=term-missing
```

The full-suite coverage run is required by `AGENTS.md` before committing.
Inspect coverage for the ownership, synchronous completion, and exception paths;
adjust the two regression methods if a new branch is unexercised. No C rebuild
is needed for these Python-only changes. Follow repository style: imports at
module scope, explicit `None` checks, simple control flow, and 89-column Python
formatting. Avoid unrelated formatting or test cleanup.

Acceptance requires one invocation per outstanding key, successful waiters
sharing the result, failed waiters terminating with the loader error, correct
byte/load accounting, and continued progress in the existing deadlock test.
One invocation for the entire cache lifetime is not required when legitimate
eviction and later demand cause a reload.

After the code fix, separately rerun the affected Full HMM benchmarks with the
same inputs, worker counts, and comparable cold/warm cache conditions. Compare
populate records with cache summaries, private-memory traces, and wall times;
check inferred outputs for equivalence allowing provenance differences. Record
the corrected revision and cache conditions. This is follow-up measurement,
not an additional repository test or work to perform in this planning task.
Do not promise a particular memory reduction or speedup before measuring it.

The plan deliberately excludes mutation-age accuracy changes, stitching
optimizations, multiprocessing changes, storage rechunking, and a general
memory-budget redesign. The concurrency here coordinates shared arrays and
condition state, so preserving the existing threads is the parsimonious fix.
