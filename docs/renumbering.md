# Renumbering a mesh for memory locality

A response to `docs/directions/fvcom-mesh-renumbering_ja.md`, which asks for
four things and calls the first the priority. This implements the first and
explains, with the FVCOM source, why the other three do not reach FVCOM.

## 1. What FVCOM actually does with a mesh's numbering

Two facts decide everything here, and both are in the source rather than in
anybody's expectation.

**FVCOM partitions the domain itself, at run time.** `setup_domain.F:120`:

```
!  DECOMPOSE DOMAIN BY ELEMENTS USING METIS
     CALL DOMDEC(NGL,NVG,NPROCS,EL_PID,MSR)
```

`DOMDEC` calls METIS with `ncommon = 2` -- the element dual -- on the mesh as
it was read. So nothing written beside the mesh changes the partition: a
`.dgraph` file, vertex weights and halo statistics would all be ignored,
because FVCOM never looks for them. That is items 2 to 4 of the request, and
they are not implemented. They would matter for a different consumer, or for
a modified FVCOM that took its partition from outside.

**The file's numbering is the memory layout inside each rank.** `genmap.F:76`
builds every rank's local numbering by walking the global ids in order and
keeping the ones it owns:

```
      DO I=1,NGL
         IF(EL_PID(I) == MYID) THEN
            N = N + 1
            NTEMP(N) = I
```

and the node map at line 94 does the same over `1..MGL`. Each rank's arrays
are therefore the file's order, filtered -- so two neighbours far apart in
the file are far apart in that rank's memory. That is item 1, and it is real.

There is no ordering constraint to respect: FVCOM reads the open-boundary
list as given and maps it through `NLID`, so any permutation is legal.

## 2. Why this project needs it

Not as a polish step. A local refinement keeps its retained elements first
and in base order -- that is what makes the frozen zone verifiable, and it is
worth keeping -- and appends the patch's new nodes at the end. Every edge
between a new node and a retained one then spans thousands of indices.

| mesh | nodes | bandwidth | mean \|i-j\| | p99 |
|---|---:|---:|---:|---:|
| goto2023 base | 3,210 | **80** | 32.6 | 74 |
| the same base, refined at Banzu | 5,409 | **3,596** | 369.4 | 2,014 |

The base is *well* numbered. The refinement took its bandwidth up by a factor
of 45. Renumbering puts back what the patch took out.

## 3. Where the sweep starts matters more than the algorithm

The first implementation used scipy's `reverse_cuthill_mckee`, which chooses
its own pseudo-peripheral starting node. The owner then pointed out that SMS
does something else: **select the open-boundary nodestring, then renumber.**

That is a different design, and measurement says it is the better one.

| starting front | goto2023 base | Banzu refined |
|---|---|---|
| the mesh as it is | 80 / 32.6 / 74 | 3,596 / 369.4 / 2,014 |
| scipy's automatic choice | 128 / 43.2 / 125 | 195 / 79.2 / 190 |
| **the open boundary (SMS)** | **78 / 33.9 / 73** | **195 / 80.3 / 190** |
| the whole boundary, coast included | 781 / 294.2 / 761 | 802 / 239.1 / 781 |

(bandwidth / mean |i-j| / p99)

Three things fall out of that table.

**The production base was numbered this way.** Seeded at the open boundary,
RCM reproduces the base's own bandwidth -- 78 against 80 -- where an
automatic start gives 128. The base's node 1 sits on the coastline at the
head of the bay, 59 % of its first hundred ids are boundary nodes against
24 % overall, and its thirteen OBC nodes are ids 3015-3150 of 3210 -- a
sweep that ENDED at the mouth, which is what the reversal in RCM does to an
OBC-seeded one.

**Seeding the whole boundary would be a mistake.** A front 781 nodes wide
gives a bandwidth of 781, ten times the mesh's own, and ordering the front
along the coastline walk instead of by degree does not help (784). A
selection marks where to start, not what to number first -- so "select the
open boundary and the coastline" cannot be what SMS does with the selection.

**The convention has a second payoff, which may be its real motive.** The
owner's reading is that SMS's workflow exists to make the open boundary
condition easier to set up, not to reduce bandwidth. The measurement supports
that: an OBC-seeded RCM leaves the open-boundary nodes as ids 3198-3210 --
**one contiguous block at the end of the mesh**. The production base is
close to it (3015-3150, clustered but not contiguous). So the convention
gives both, and there is no trade to make.

`--seed obc` is therefore the default, `--seed auto` keeps scipy's choice,
and `--seed boundary` exists so the claim above can be checked rather than
believed.

## 4. Morton and Hilbert

Both are implemented because both were asked for. Neither won anything.

| grid | Morton | Hilbert | RCM |
|---|---|---|---|
| 17 x 17 | 145 / 11.8 | 227 / 12.6 | **17 / 8.1** |
| 33 x 33 | 545 / 22.4 | 879 / 25.3 | **33 / 15.3** |
| 65 x 65 | 2,113 / 43.7 | 3,463 / 50.9 | **65 / 29.5** |

(bandwidth / mean |i-j|)

The textbook argument is that Hilbert beats Morton because the Z curve jumps
the width of the domain at every power-of-two boundary. On these meshes it
does not -- Morton's mean edge span is the better of the two -- and neither
comes close to RCM, which optimises the quantity being measured rather than
a proxy for it. The docstrings say so, and the tests assert the measurement
rather than the expectation.

## 5. The guard

`renumber_mesh` refuses to renumber a mesh it would make worse, and the
criterion is the mean **and** the tail: Morton on the goto2023 base improves
the mean edge span, 32.6 to 27.8, while taking the 99th percentile from 74 to
610. A mesh whose typical neighbour is nearer and whose worst neighbours are
eight times further away has not been improved. `--force` overrides; the base
is left alone by default, which is the correct answer for it.

## 6. Does it make FVCOM faster?

Open, and honestly so. The first attempt did not answer it: the two cases ran
concurrently on one shared node, once each, and the renumbered mesh came back
**14 % slower** (545 s against 477 s). Two things were wrong with that as an
experiment, one of which should have been obvious beforehand:

* they shared the node, so each one's time depended on the other's;
* **renumbering changes the partition.** METIS is not invariant under a
  permutation of its input, so a 64-rank comparison measures the memory
  layout *and* whatever decomposition METIS happened to choose. It is not a
  controlled experiment for the claim in the request.

`jobs/octopus/416_renumber_benchmark.sh` is the controlled one: the cases run
sequentially inside one job, alternating order, repeated, at **1 rank** --
where `setup_domain.F` skips `EL_PID` entirely and there is no partition, so
only the layout differs -- and at 64, which is what production would see.

Physics agreed in the uncontrolled run, for what that is worth: 961 records
each, and the five gauges matched to 3 micrometres in amplitude.
