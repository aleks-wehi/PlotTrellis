# Gene / protein info panel on the composer's Output page: the UniProt link-out,
# the per-transcript annotation table (transcript_annotation.R, plus the UniProt
# canonical-isoform check in uniprot_canonical.R) and the canonical protein
# annotation table (protein_domains.R). Not part of the figure or the downloaded
# R script.
#
# It always describes the gene of the RENDERED figure (state$rendered_gene, set by
# do_render() in app.R; zoom and Apply keep it), never the gene currently typed on
# Compose, and is hidden for a coordinates-only figure.
#
# Loading: the transcript table comes from the local TxDb/lookup files, so it is
# built as soon as a figure is rendered, before the panel is opened. The UniProt
# parts can need a network call. When the gene is already in their disk cache
# they are filled in straight away; otherwise the panel shows a loading note and
# the fetch runs in session$onFlushed(), i.e. after that note has reached the
# browser, and only once the panel (or tab) is actually open, so a figure whose
# panel nobody opens never triggers a UniProt call.
#
# Nothing here reads composer globals: everything comes in as arguments.

# The Output-page panel. `panel` is the composer's collapsible_panel() shell.
.pt_gene_info_ui <- function(panel) {
  conditionalPanel("output.gene_info_gene",
    tags$div(style = "margin-top:12px;",
      panel("gene_protein_info", "Gene / protein info", collapsed = TRUE,
        tags$p(class = "text-muted", "For ", tags$strong(textOutput("gene_info_gene", inline = TRUE)),
               ", the gene of the figure above."),
        tabsetPanel(type = "pills",
          tabPanel("Gene info",
                   # Shown until the first value arrives.
                   uiOutput("gr_extras_ui", container = function(...) tags$div(...,
                     tags$div(class = "text-muted", style = "font-style:italic;", "Loading transcript annotations...")))),
          tabPanel("Canonical protein annotation", uiOutput("gr_protein_ui"))
        )
      )
    )
  )
}

.PT_GI_UNIPROT_LOADING <- "Fetching UniProt annotations (first time for this gene can take a few seconds)..."

# Whether the UniProt canonical-isoform cache already holds this gene (same key as
# uniprot_canonical_match()).
.pt_gi_uniprot_cached <- function(gene, orgdb) {
  !is.null(.uniprot_canonical_cache_load()[[paste0(orgdb, "|", gene)]])
}

# Whether gene_protein_domains() has this gene on disk (same path as its cache_file).
.pt_gi_pd_cached <- function(gene, genome_label) {
  file.exists(file.path(PROTEIN_DOMAINS_CACHE_DIR, genome_label, paste0(gene, ".rds")))
}

# genome: list(orgdb_pkg, tx_name_lookup_path, tx_meta_lookup_path, gtf_path,
# genome_label). reviewed_only: a reactive giving the live uniprot_reviewed_only.
.pt_gene_info_server <- function(input, output, session, state, txdb, cds_by_tx, genome, reviewed_only) {
  gene <- reactive({
    g <- state$rendered_gene
    if (is.null(g) || length(g) != 1 || is.na(g) || !nzchar(g)) NULL else g
  })

  output$gene_info_gene <- renderText(gene() %||% "")
  outputOptions(output, "gene_info_gene", suspendWhenHidden = FALSE)

  # Built when the figure is rendered, not when the panel opens (see the header).
  ann_r <- reactive({
    g <- gene()
    if (is.null(g)) return(NULL)
    tryCatch(list(ann = transcript_annotations(g, txdb, genome$orgdb_pkg, genome$tx_name_lookup_path,
                                               genome$tx_meta_lookup_path)),
             error = function(e) list(error = conditionMessage(e)))
  })

  # Session results of deferred fetches, by gene: list(value, failed, stamp). A failed
  # (uncached) lookup is retried after the next render of the figure.
  uc_store <- reactiveVal(list())
  pd_store <- reactiveVal(list())
  pending  <- new.env()
  stored <- function(store, g) {
    hit <- store()[[g]]
    if (!is.null(hit) && hit$failed && !identical(hit$stamp, state$svg_path)) NULL else hit
  }
  is_open <- function(id) !isTRUE(session$clientData[[paste0("output_", id, "_hidden")]])
  # Runs fetch() once this flush has sent the loading note, then stores the result.
  # later() rather than fetching inside onFlushed itself, so the server loop
  # re-renders straight after (it flushes after running later callbacks).
  defer <- function(store, tag, g, fetch, cached) {
    key <- paste(tag, g)
    if (isTRUE(pending[[key]])) return(invisible())
    pending[[key]] <- TRUE
    session$onFlushed(function() later::later(function() {
      on.exit(pending[[key]] <- FALSE, add = TRUE)
      if (session$isClosed()) return()
      value <- fetch()
      s <- isolate(store())
      s[[g]] <- list(value = value, failed = !cached(), stamp = isolate(state$svg_path))
      store(s)
    }), once = TRUE)
    invisible()
  }

  uc_fetch <- function(g, ann) tryCatch(
    uniprot_canonical_match(g, ann$tx_id, cds_by_tx, ann$ccds_id, orgdb = genome$orgdb_pkg),
    error = function(e) NULL)
  pd_fetch <- function(g) tryCatch(
    gene_protein_domains(g, txdb, genome$gtf_path, genome$orgdb_pkg, genome$genome_label),
    error = function(e) list(status = "error", message = conditionMessage(e)))

  output$gr_extras_ui <- renderUI({
    gene <- gene()
    if (is.null(gene)) return(NULL)
    res <- ann_r()
    ann <- res$ann
    uc_loading <- FALSE

    if (!is.null(ann) && nrow(ann) > 0) {
      # UniProt check is a live network call (see uniprot_canonical.R) --
      # failure degrades to NA/"unavailable" per row, never blanks out the
      # GTF-derived columns above, which don't depend on it.
      hit <- stored(uc_store, gene)
      uc <- if (!is.null(hit)) hit$value
            else if (.pt_gi_uniprot_cached(gene, genome$orgdb_pkg)) uc_fetch(gene, ann)
            else {
              uc_loading <- TRUE
              if (is_open("gr_extras_ui"))
                defer(uc_store, "uc", gene, function() uc_fetch(gene, ann),
                      function() .pt_gi_uniprot_cached(gene, genome$orgdb_pkg))
              NULL
            }
      if (is.null(uc) || nrow(uc) != nrow(ann)) {
        uc <- data.frame(uniprot_accession = rep(NA_character_, nrow(ann)),
                          uniprot_match = rep(NA, nrow(ann)),
                          uniprot_unique = rep(NA, nrow(ann)),
                          uniprot_match_count = rep(NA_integer_, nrow(ann)),
                          uniprot_source = rep(NA_character_, nrow(ann)),
                          uniprot_canonical_isoform_id = rep(NA_character_, nrow(ann)),
                          uniprot_canonical_isoform_name = rep(NA_character_, nrow(ann)),
                          uniprot_row_isoform_id = rep(NA_character_, nrow(ann)),
                          uniprot_row_isoform_name = rep(NA_character_, nrow(ann)))
      }
      ann[names(uc)] <- uc
    }

    tag_pill <- function(text) tags$span(
      style = paste("display:inline-block; background:#eef2f7; color:#345; border-radius:3px;",
                     "padding:1px 6px; margin-right:4px; font-size:0.9em;"),
      text
    )

    uniprot_source_label <- c(
      ensembl_xref = "direct UniProt-Ensembl cross-reference",
      ccds_xref    = "UniProt-CCDS cross-reference",
      length_match = "protein length match only, not a direct cross-reference"
    )

    tagList(
      tags$a(href = uniprot_search_url(gene, genome$orgdb_pkg, isTRUE(reviewed_only())), target = "_blank",
             class = "btn btn-default btn-sm",
             paste0("UniProt: ", gene, " (opens in new tab)")),
      if (!is.null(res$error))
        tags$div(class = "text-danger", style = "margin-top:8px;",
                 "Couldn't load the transcript annotations: ", res$error),
      if (uc_loading)
        tags$div(class = "text-muted gene-info-uniprot-loading", style = "margin-top:8px; font-style:italic;",
                 .PT_GI_UNIPROT_LOADING),
      if (!is.null(ann) && nrow(ann) > 0) {
        tags$div(style = "margin-top:8px; font-size: var(--bio-font-code, 12px);",
          tags$div(style = "display:flex; gap:10px; font-weight:600; color:#666; border-bottom:1px solid #ccc; padding-bottom:3px; margin-bottom:2px;",
            tags$span(style = "flex:0 0 90px;", "Transcript"),
            tags$span(style = "flex:0 0 60px;", "Type"),
            tags$span(style = "flex:1;", "Tags")
          ),
          lapply(seq_len(nrow(ann)), function(i) {
            src <- ann$uniprot_source[i]

            # The specific UniProt sequence this whole check compares
            # against -- an entry-level fact (isoform id + name, or just the
            # bare accession for a single-isoform entry), true regardless of
            # whether THIS transcript carries a direct tag for it.
            canon_id   <- ann$uniprot_canonical_isoform_id[i]
            canon_name <- ann$uniprot_canonical_isoform_name[i]
            canon_desc <- if (is.na(canon_id)) NA_character_ else
              if (is.na(canon_name)) canon_id else sprintf("isoform %s (%s)", canon_name, canon_id)

            uni_text <- if (uc_loading) {
              "loading..."
            } else if (is.na(ann$uniprot_match[i])) {
              "unknown (UniProt lookup unavailable)"
            } else if (ann$uniprot_match[i] && ann$uniprot_unique[i]) {
              sprintf("yes, this is UniProt's canonical %s (via %s)", canon_desc, uniprot_source_label[[src]])
            } else if (ann$uniprot_match[i]) {
              sprintf("matches UniProt's canonical %s (via %s), but %d transcripts of this gene share that exact sequence, so it doesn't uniquely identify this one",
                      canon_desc, uniprot_source_label[[src]], ann$uniprot_match_count[i])
            } else if (!is.na(src) && src != "length_match") {
              row_id   <- ann$uniprot_row_isoform_id[i]
              row_name <- ann$uniprot_row_isoform_name[i]
              row_desc <- if (is.na(row_id)) ann$uniprot_accession[i] else
                if (is.na(row_name)) row_id else sprintf("isoform %s (%s)", row_name, row_id)
              sprintf("no, this is %s, confirmed by %s; UniProt's canonical is %s",
                      row_desc, uniprot_source_label[[src]], canon_desc)
            } else if (!is.na(ann$uniprot_accession[i])) {
              sprintf("no; UniProt's canonical is %s", canon_desc)
            } else {
              "no reviewed UniProt entry found"
            }

            row_tags <- character(0)
            if (isTRUE(ann$ensembl_canonical[i])) row_tags <- c(row_tags, "Ensembl canonical")
            if (!is.na(ann$appris_tier[i])) row_tags <- c(row_tags, paste("APPRIS", ann$appris_tier[i]))
            if (isTRUE(ann$uniprot_match[i])) {
              row_tags <- c(row_tags, if (isTRUE(ann$uniprot_unique[i])) "UniProt canonical" else "UniProt (shared match)")
            }

            tags$details(
              style = "border-bottom: 1px solid #eee; padding: 3px 0;",
              tags$summary(
                style = "display:flex; gap:10px;",
                tags$span(style = "flex:0 0 90px; font-weight:600;", ann$transcript_name[i]),
                tags$span(style = "flex:0 0 60px;", ann$coding_status[i]),
                tags$span(style = "flex:1;", lapply(row_tags, tag_pill))
              ),
              tags$div(style = "margin: 4px 0 4px 14px; color:#444;",
                tags$div(tags$strong("Ensembl canonical: "), if (ann$ensembl_canonical[i]) "yes" else "no"),
                tags$div(tags$strong("APPRIS: "), if (is.na(ann$appris_tier[i])) "none" else ann$appris_tier[i]),
                tags$div(tags$strong("CCDS: "), if (is.na(ann$ccds_id[i])) "none" else ann$ccds_id[i]),
                tags$div(tags$strong("UniProt canonical: "), uni_text),
                tags$div(style = "margin-top:3px; font-style:italic; color:#888;",
                  "Source: GENCODE annotation tags, UniProt cross-reference (see UniProt canonical line above for which kind).")
              )
            )
          })
        )
      }
    )
  })
  outputOptions(output, "gr_extras_ui", suspendWhenHidden = FALSE)

  # ── Canonical protein annotation: UniProt domain/region/motif deep dive ─────
  # Same progressive-disclosure style as gr_extras_ui above (one <details> per
  # feature). Uses gene_protein_domains() (protein_domains.R) directly -- the
  # same resolver the "Canonical protein annotation" TRACK calls, read from its
  # per-gene RDS cache after the first lookup of a gene. One row per domain_id
  # (a multi-exon feature collapsed to one row, not one row per exon segment) --
  # matches how the track itself groups them.
  output$gr_protein_ui <- renderUI({
    gene <- gene()
    if (is.null(gene)) return(NULL)

    hit <- stored(pd_store, gene)
    pd <- if (!is.null(hit)) hit$value
          else if (.pt_gi_pd_cached(gene, genome$genome_label)) pd_fetch(gene)
          else NULL
    if (is.null(pd)) {
      if (is_open("gr_protein_ui"))
        defer(pd_store, "pd", gene, function() pd_fetch(gene),
              function() .pt_gi_pd_cached(gene, genome$genome_label))
      return(tags$div(class = "text-muted gene-info-uniprot-loading", style = "font-style:italic;",
                      .PT_GI_UNIPROT_LOADING))
    }

    pd_msg <- switch(pd$status,
      no_uniprot_entry     = "No reviewed UniProt entry found for this gene",
      no_domains           = "No annotated domain/region/motif features for this gene's canonical isoform",
      no_canonical_transcript = "UniProt's canonical isoform could not be matched to any annotated transcript for this gene",
      ok = NULL,
      "Protein domain lookup failed (network issue?) -- will retry on next render"
    )

    method_label <- if (identical(pd$status, "ok")) switch(pd$method,
      ccds_id_join                    = "via CCDS ID",
      length_match_fallback           = "via protein-length match (CCDS not tagged in GTF yet)",
      length_match_fallback_ambiguous = "via protein-length match (ambiguous, CDS shared by multiple transcripts)",
      pd$method
    ) else NULL

    tagList(
      if (!is.null(pd_msg))
        div(style = "color:#888; font-style:italic;", pd_msg)
      else {
        domains <- pd$domains
        dids <- unique(domains$domain_id)
        dids <- dids[order(vapply(dids, function(d) min(domains$start[domains$domain_id == d]), numeric(1)))]

        tagList(
          div(style = "font-size:0.9em; color:#666; margin-bottom:6px;",
              sprintf("canonical: %s (%s)", pd$canonical_tx, method_label)),
          tags$div(style = "font-size: var(--bio-font-code, 12px);",
            tags$div(style = "display:flex; gap:10px; font-weight:600; color:#666; border-bottom:1px solid #ccc; padding-bottom:3px; margin-bottom:2px;",
              tags$span(style = "flex:0 0 14px;", ""),
              tags$span(style = "flex:1;", "Feature"),
              tags$span(style = "flex:0 0 60px;", "Exons")
            ),
            lapply(dids, function(did) {
              seg  <- domains[domains$domain_id == did, ]
              seg  <- seg[order(seg$start), ]
              ty   <- seg$type[1]
              col  <- PD_TYPE_COLOURS[[ty]]
              nseg <- nrow(seg)

              tags$details(
                style = "border-bottom: 1px solid #eee; padding: 3px 0;",
                tags$summary(
                  style = "display:flex; gap:10px; align-items:center;",
                  tags$span(style = sprintf(
                    "flex:0 0 12px; width:12px; height:12px; border-radius:2px; background:%s; display:inline-block;",
                    col
                  )),
                  tags$span(style = "flex:1;", seg$note[1]),
                  tags$span(style = "flex:0 0 60px; color:#888;", if (nseg > 1) sprintf("%d exons", nseg) else "")
                ),
                tags$div(style = "margin: 4px 0 4px 22px; color:#444;",
                  tags$div(tags$strong("Type: "), ty),
                  tags$div(tags$strong("Position: "), sprintf("aa %d-%d", seg$aa_start[1], seg$aa_end[1])),
                  tags$div(tags$strong("Exon segments: "), nseg)
                )
              )
            })
          )
        )
      }
    )
  })
  outputOptions(output, "gr_protein_ui", suspendWhenHidden = FALSE)
}
