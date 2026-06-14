# Status — gmulti federation implementation

Companion to [DESIGN_gmulti_federation.md](DESIGN_gmulti_federation.md). This file tracks **only the four implementation phases currently in scope** — what's landed, what's planned, what's deferred.

**Branch**: `feature/gmulti-federation-design` (off `feature/giotto-view`)
**Worktree**: `/Users/george/Documents/GitHub/GiottoClass-federation-design`
**Test count baseline**: 174 PASS on parent branch (before this work); 242 PASS as of phase 3 (+ 1 pre-existing snapshotSave failure unrelated to this work).

---

## Landed

### Phase 1 — `@mapping` slot + auto-discovery + accessor

Commits `057f5135` (slot + auto-discovery + full-list setter) and `fc0b631a` (axis-scoped + entry-scoped setters).

- New slot `gmulti@mapping`: two-list structure (`spat_unit` / `feat_type`), each holding per-sample-named char vectors mapping gmulti-level handle → child-level slot name. Declares which children participate in each (spat_unit / feat_type) and reconciles per-child name variation.
- `.gm_discover_mapping(children)`: auto-populates the symmetric trivial mapping from children's `@cell_ID` / `@feat_ID` slot keys at construction. User-edited mappings survive bare re-init; assigning `NULL` triggers fresh re-discovery.
- `gmultiMapping()` / `gmultiMapping<-` accessor (exported generics). Three setter forms:
  - `gmultiMapping(mg) <- full_list` — full replacement
  - `gmultiMapping(mg, "spat_unit") <- axis_list` — replace one axis
  - `gmultiMapping(mg, "spat_unit", "cell") <- c(B191 = "cell", B215 = "poly")` — replace one entry
- Setter validates per-sample names + child slot existence; rejects unknown entries with clear errors.
- Joint-slot invalidation on mapping mutation is **per-universe** — only the affected `(spat_unit, feat_type)` universe's joint state drops; unrelated universes survive byte-identically.
- `show()` displays a one-line summary per axis.

### Phase 2 — Federation consults `@mapping`

Commit `f1d7beb8`.

- New `.gm_resolve_axis(g, axis, handle)`: returns participating samples + per-sample child-level slot names. Primary path is `@mapping` lookup; legacy fallback scans child slots when the handle isn't declared.
- New `.gm_resolve_participation(g, spat_unit, feat_type)`: composes two axis resolutions into the per-sample `(su, ft)` federation plan.
- Three federation helpers refactored to delegate via the new helpers:
  - `.gm_assemble_expression`
  - `.gm_assemble_cell_metadata`
  - `.gm_assemble_feat_metadata`
- `@mapping` is now load-bearing — a user-edited mapping unifying e.g. B's `"transcripts"` feat_type under gmulti-level `"rna"` makes `getCellMetadata(mg, feat_type = "rna")` actually pull B's transcripts slot.
- Legacy gmulti objects without `@mapping` (mapping cleared) keep working via per-child default fallback.

### Phase 3 — Access layer: `sample =` arg + `"sample::name"` parser

Commit `71a8e6a5`.

- `.parse_sample_qualified_name(name)`: splits on the first `::` for optional sample prefix.
- `.gm_slice_to_sample(x, sample, gobject)`: per-subobject-class slicer for joint cmeta / expr / dimreduc / nnnet via the `sample::cell_id` namespacing convention.
- `getExpression(gmulti)`: adds `sample =` AND parses `values =` for `"sample::name"` prefix. Conflicting `sample =` + prefix errors clearly.
- `getCellMetadata(gmulti)` and `getFeatureMetadata(gmulti)`: add `sample =` arg. (featmeta sample is no-op slice + validation only since featIDs are passthrough.)
- Five spatial-domain getters (`getSpatialLocations`, `getSpatialNetwork`, `getPolygonInfo`, `getFeatureInfo`, `getGiottoImage`): `sample =` as canonical alias for the legacy `object =` arg. Conflicts error.

After phase 3, the access pattern looks like:

```r
# canonical
getCellMetadata(mg, sample = "B191")
getExpression(mg, sample = "B191", values = "raw")

# prefix shortcut (on getters with a name/values arg)
getExpression(mg, values = "B191::raw")

# spatial domain (sample is alias of legacy object=)
getSpatialLocations(mg, sample = "B191")
```

---

## Planned — phase 4 (next up)

**Dispatcher cleanup: `samples =` arg + drop `:::` reach into GiottoClass internals.**

Touches: `GiottoVisuals/R/gmulti.R` (the `.gg_multi_dispatch_spatial` helper), and the top-level spatial plot fns that accept `space =` for sample selection.

Concrete changes:
1. Rename `space` → `samples` in `.gg_multi_dispatch_spatial`. Today `space =` is typed as a character vector of sample names — semantic misnomer. After phase 4, `space =` only takes real defined-space names from `@spaces`.
2. Implement `.resolve_samples(gobject, samples, space)`: when `space = "atlas"` is passed and `samples = NULL`, derive `samples` from `names(atlas@samples)` (the auto-injection convention already documented in the view/space design memory).
3. Drop the scratch-child injection. Per-panel slicing now goes through the access-layer `sample =` arg from phase 3:
   ```r
   plots <- lapply(samples, function(s) {
       a <- named
       a$gobject <- gobject  # gmulti, not a child
       a$sample <- s
       do.call(plot_fn, c(a, dots))
   })
   ```
4. Delete `.gm_inject_joint_metadata` from `R/gmulti.R` (in GiottoClass). The only caller was the GiottoVisuals dispatcher; with phase 3 in place that caller doesn't need it anymore.
5. Drop the `GiottoClass:::` reference in `GiottoVisuals/R/gmulti.R:99`. R CMD check NOTE on GiottoVisuals goes away.
6. Update each spatial plot fn (`spatPlot2D`, `spatInSituPlotPoints`, `dimPlot2D`, etc.) to accept `samples =` and pass it through to the dispatcher. Update plot fn signatures to also accept `sample =` (forwarded to getters internally) for the per-panel call from the dispatcher.

Tests:
- Dispatcher receives `samples =` and iterates correctly.
- `space = "atlas"` derives samples via auto-injection.
- Conflicting `samples =` + `space@samples` keys (sample not in atlas) error.
- Plot fn output unchanged from current behavior on existing test cases.
- GiottoVisuals R CMD check no longer flags the `:::` use.

Estimated scope: ~50 lines of edits across two repos (GiottoClass + GiottoVisuals), bulk in plot fn signatures. No new classes, no new slots.

---

## Deferred — phases 5–6 (NOT in scope right now)

These are real follow-ons but not blocking. Land them when a concrete workflow needs them.

### Phase 5 — Federated-read wrapper class

When getters like `getExpression(mg, sample = NULL)` return cross-sample federations, today they're either an eager combined object (current behavior post-assembly) or a list-of-substores when assembly hasn't happened.

The wrapper class (`federatedReadHandle` or similar) would let the materialization decision defer to the consumer — duckdb/sedonadb queries lower into a single SQL plan over the list-of-substores; arrow consumers concat at use time. Mirrors `unionParquetGeomStore` in GiottoDisk.

Why deferred: list-of-substores works today via the existing assembly path. The wrapper class is a perf/laziness optimization for specific consumers, not a correctness fix.

### Phase 6 — Pointer-class for `@spatial_info`

A `gmultiSpatialAlias` class for ad-hoc cross-sample groupings within `@spatial_info` that don't follow the federation pattern (e.g. `"tumor_focus_polys" = B191's tumor_roi + B215's epithelium_roi`).

Why deferred: `@mapping` covers the standard federation case for spat_unit / feat_type axes. Pointer-class only matters for content-level aliases that don't decompose along those axes — long-tail use case, not blocking.

---

## Out of scope for THIS implementation effort (separate design pass)

Items the design doc explicitly identifies as needing their own design — do not pull into this work:

- **Composable views** (`view1 + view2`). Memory line 125 in `project_giottoview_design_shape.md`.
- **Cross-sample aggregation infrastructure** (the "fan-out + reduce" dispatch shape — needed for atlas-frame polygon aggregation across samples).
- **GiottoLens consumption of `@mapping`** — downstream package, separate branch.
- **Migration helper** for existing saved gmulti objects without `@mapping`. Auto-init on load might cover it; defer the decision until a concrete migration case appears.
- **`@h5_file` slot deprecation cleanup** — adjacent but independent.

---

## Files touched so far

```
R/gmulti.R                         (+ ~510 lines net across phases 1-3)
NAMESPACE                          (+ 4 exports: gmultiMapping, gmultiMapping<-)
man/giottoMulti-class.Rd           (regenerated)
man/gmultiMapping.Rd               (new)
tests/testthat/test-gmulti.R       (+ 26 new test_that blocks)
```

No changes outside GiottoClass yet. Phase 4 will touch GiottoVisuals as well.

---

## Pointers

- Full design rationale + sdata comparison: [DESIGN_gmulti_federation.md](DESIGN_gmulti_federation.md)
- Memory entry that indexes both files: `project_gmulti_federation_design.md`
- Foundation memory (view/space split — load-bearing for §2 of the design doc): `project_giottoview_design_shape.md`

*Last updated after phase 3 (2026-06-14). Update on phase 4 landing.*
