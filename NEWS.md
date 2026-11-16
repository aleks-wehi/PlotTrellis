# PlotTrellis 1.0.0

## Breaking changes

- `pt_figure()` takes plain track names plus `samples =`, `comparison =` and `order =`. The
  `"atac@LP"` form and the `parallel` argument are gone. `pt_figure(proj, "Elf5", "atac@LP")`
  stops with a message showing the new call to write, and so does `parallel =` (named, or a
  TRUE or FALSE in its old fourth place). Rewrite old scripts. Tracks and samples can also be
  given by number: `pt_tracks()` now has one row per track with an `id`, and the new
  `pt_samples()` lists the samples. A samples.txt can set the opening view with
  `# start_samples:`, `# start_comparison:` and `# start_order:` lines.
- Splice-junction files keep their TPM-scaled comparison track, but loop files (paired-end
  contacts, such as ChIA-PET or HiCCUPS loops) no longer get a `_diff_tpm` track, and their labels
  no longer say "reads". A saved script that names a loop file's `_diff_tpm` track now stops
  with "Unknown track". A `# start_tracks:` line that names one (or any track the project does
  not have, a misspelling included) leaves that track out with one warning naming it and keeps
  the rest of the start view (if none are left, there is no opening view).
- Figures are drawn in a second R process (it needs the `callr` package), started when the app
  launches and kept running. That is what lets the new Cancel button stop a render, and the
  app stays responsive while a figure draws. It uses about as much memory again as the app
  itself (about 1.7 GB peak on the example), a figure asked for within about 20 seconds of
  starting waits for the new process to load its packages, and the first figure after a Cancel
  takes about 15 to 20 seconds longer for the same reason. Options set in your own `.Rprofile`
  (a proxy, timeouts) do not reach that drawing process. After a cancelled zoom the region boxes
  show the new region while the figure is still the old one.
- PlotTrellis needs R 4.5.0 or later and plotgardener 1.16.0 or later. With older versions the
  transcript tracks failed with "unused argument (transcriptFilter = ...)".
- Genome files built from a GTF (the annotation database, transcript lookups, chromosome
  lengths) now go in the project's own cache folder (`cache/genome/<hash>/`, one subfolder per GTF
  file, a new one if the file is replaced or touched), never beside the GTF. Each existing
  project builds them once more after the upgrade (1 to 3 minutes for a full annotation) and
  fetches chromosome lengths from UCSC, so it needs the internet that once. The old files next
  to your GTF are no longer used and can be deleted. "Clear all caches" now only touches the
  project's cache folder.
- A blank genome label is no longer accepted: the wizard's Save and a samples.txt import refuse
  it, and so does opening a project whose project_settings.json has one (`pt_project()` and
  `plottrellis_composer()` stop with the same message, naming the file). `run_plottrellis()`
  opens the setup wizard on such a project instead, saying why in the console. A blank label
  used to be fine while the genome's chromosome lengths were already cached (or shipped, like
  the examples) and only stopped the first launch of a new GTF; now it is refused every time.
  A samples.txt with no genome label line still gets the default (mm39); only an explicitly
  blank one is refused. The message says to enter the UCSC name of the build, for example mm39,
  mm10 or hg38.

## New and improved

- The setup wizard is a guided first-time setup: a landing screen, numbered Genome, Data folder
  and Samples sections, a "Which GTF? Where do I get it?" guide, and live checks under the GTF,
  OrgDb and genome label fields as you type. Help icons open on hover or click.
- Importing a samples.txt that has one column per sample (instead of one column per kind of
  data) now warns, names those columns and shows the layout to use.
- A transcript has the same colour in every transcript track of a figure.
- The app has a new look drawn from the logo: dark-green title bar, cream page, green action
  buttons. Nothing behaves differently.
- Slow or one-off buttons (Save, Add, Remove, Import, Export, Clear caches and others) run their
  action once when clicked twice quickly, and are greyed out until it finishes.
- "Clear all caches" and "Reset project configuration" are in a closed "Advanced" part of the
  Edit project settings dialog, which also shows the project folder.
- The vignette has a troubleshooting section (what to delete to rebuild or start fresh), a
  "Your genome" section (which GTF, where to get it, which release, which
  OrgDb package and genome label), an explanation of the samples.txt columns, and a prompt that
  has an AI assistant write a samples.txt with you (also a second vignette,
  `vignette("samples-txt-ai-prompt")`). Its code blocks have a Copy button that stays put when a
  long line scrolls.
- A short README and a CONTRIBUTING.md (how to report a bug).
- Figures can be cancelled: the loading screen has a Cancel button that stops the render within
  a second, keeps the figure you had and says it was cancelled (Generate, zoom and Apply to
  figure). See the breaking changes for what this costs.
- `pt_figure(comparison =)` or `pt_figure(samples =)` on a figure with no track that needs a
  sample now warns "has no effect" instead of being ignored silently. The figure is still made.
- Before a sample is picked, the empty Available tracks list says what to pick first.
- After a failed Save, the wizard's "Not saved yet" list follows the form: a problem you fix
  leaves the list straight away, and the checks run afresh on the next Save. Save into a
  project folder you cannot write to shows a message with the reason instead of ending the
  session, and settings (also "Save as project default" in the app) are written through a
  temporary file so a failed write keeps the old settings.
- A 7-column loop file whose last column is a count (ENCODE ChIA-PET) draws each loop at that
  height; the vignette says which column sets the arc height.
- Saving the setup wizard without changes keeps the project's opening view (the gene, samples
  and tracks the app starts on); it is dropped only if its genome changed or one of its samples
  or tracks was removed, and a samples.txt imported in the wizard replaces it with the view its
  `# start_` lines give (none if it has none), as `pt_project()` does. project_samples.txt now
  carries that view as start lines after its other `#` lines, so the restore offered after
  "Reset project configuration" gives it back. Start lines that do not fit the samples do not
  stop the Save: there is no opening view, and the "Settings saved" dialog says why.
- A Hi-C colour cap of 0 or below, or a junction opacity cap of 0 or below, says so in that
  track ("the Hi-C colour cap must be above 0 (it is 0)") instead of a red drawing error or
  fully opaque arcs.
- What a figure prints while it draws (a failed track's reason, notes such as "named
  transcripts leave out ...", warnings) shows in the R console again when the figure is done,
  fails or is cancelled, as it did before figures moved to the second R process. The lines
  arrive all at once at the end, not one by one.
- The composer's track filter shows buttons only for the track types the project has (the
  mammary example no longer offers BED, HiC or Junctions).
- Several R processes opening the same new Ensembl project at once no longer fail with
  "database is locked": each waits up to 10 seconds for another's chromosome rename and checks
  again before renaming. If one waits longer, or the file is read-only, it warns and says which.
- The app starts without bslib's "Navigation containers expect a collection of nav_panel()"
  warning, which used to print twice.

## Fixes

- A page width of 0 or below in project_settings.json takes the default (9 inches) with the
  "Ignoring unusable value ... for page_w" warning, instead of stopping the figure. A narrow positive width is kept.
- "Reset project configuration" in a project folder you cannot write to says "Nothing was
  reset" with the reason and the app stays as it was, instead of saying the file was deleted
  and restarting into the unchanged project.
- After a failed Save in the wizard ("Couldn't write ..."), the message goes at the next Save
  or Back, and the live "Not saved yet" list works again.
- Switching the gene box to Coordinates and back keeps its error ("not found" or an invalid
  flank) instead of clearing it.
- After "Clear all caches", a figure drawn in R from a project opened before the clear (a `pt_project()` object) says to run
  `pt_project()` on it again instead of "DB file ... not found". The dialog says "All caches
  have been deleted" only when they were, lists what it could not delete, and no longer leaves
  a .nfs file behind on a network drive.
- A .cool Hi-C file whose last chromosome has no contacts shows "this file has no contacts on
  ..." in the track instead of "index exceeds HDF5-array dimension". (A .cool file that names
  its chromosomes "1" instead of "chr1" is still not matched.)
- Rebuilding a project from its own project_samples.txt (`pt_project()` or `run_plottrellis()`
  with it) keeps its opening view, and the copy keeps its `# start_` lines. A copy written with `# start_samples: NA, NA` (an opening view of tracks that
  need no sample) imports without the "samples has no effect" warning, and such a view is
  written without a `# start_samples:` line from now on.
- Collapsing or opening a section inside another (Figure settings or R script) no longer
  flips the arrow of the section around it.
- With an Ensembl GTF, Gene info shows CCDS ids and the protein domain track finds the canonical
  transcript by CCDS id, as with GENCODE. Protein domain results already cached for an Ensembl
  GTF keep their protein-length-matched transcript until you use "Clear all caches".
- A region whose end is not after its start (or a gene flank that makes one) stops with a plain
  message instead of "argument is of length zero".
- A coverage track whose signal in view is tiny (top below 0.00005) draws on its own axis
  instead of failing.
- Clearing the gene box no longer leaves the old gene's region behind for Generate to draw.
- Saving or re-importing a project no longer rounds settings to 4 decimals.
- A transcript row opened in Gene / protein info stays open when the UniProt lookup finishes.
- The wizard's Remove and Edit buttons act on the sample on the row you clicked, however fast
  you click.
- The Gene info UniProt-CCDS cross-reference also matches CCDS ids that carry no version
  (Ensembl's), so it works with Ensembl GTFs.
- A sample whose quant.sf is missing or unreadable no longer fails every transcript track: the
  named transcripts come from the samples that can be read, only that sample's own expression
  row says which file is missing, and the R console names the sample left out.
- The per-gene caches of transcripts and protein domains remember which GTF they came from, so
  a GTF edited or replaced in place gets fresh results on the next start. Genes with no
  reviewed UniProt entry, no domain features or no transcript match are cached like any other
  result, so they stop asking UniProt on every render.
- A display setting in project_settings.json with the wrong type or length (`page_w` "wide",
  `fontsize` `[]`) falls back to its default with one warning naming it, instead of stopping
  the figure. A box the app itself saved empty takes the default without a warning.
- A Hi-C view with no values under its normalisation says so in the track, with what to
  change, instead of drawing an empty panel.
- A non-numeric Gene flank % (for example in a hand-edited project file) falls back to the
  last saved value or the default 15% instead of ending the session.
- A Hi-C view with a very small maximum draws instead of failing, and a Hi-C track that fails
  keeps its normal height.
- With an Ensembl GTF, chromosomes are named `chr1`, `chrX`, `chrM` as with GENCODE: the axis of
  a gene search reads `chr2`, not `2`, and a region picked by chromosome list or
  `region = "chr2:..."` draws its gene annotation tracks. An index built before this is
  renamed the next time the project is opened.
