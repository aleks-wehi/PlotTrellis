

    // Track-filter buttons + bucket are scoped per [data-available-id] block.
    function applyTrackFilter(blockEl) {
      var $block  = $(blockEl);
      var availId = $block.attr('data-available-id');
      var f = $block.find('.filter-btn.active').data('filter') || 'all';
      $block.find('#' + availId + ' .rank-list-item').each(function() {
        var key  = $(this).attr('data-rank-id') || '';
        var type = TRACK_TYPE_MAP[key] || '';
        var matches = f === 'all' || TRACK_FILTER_MAP[type] === f;
        if (matches) {
          $(this).show();
        } else {
          $(this).hide();
        }
      });
    }

    function applyAllTrackFilters() {
      $('[data-available-id]').each(function() { applyTrackFilter(this); });
    }

    // Flags, rather than hides, gene-required tracks while the region mode is
    // Coordinates: one may already be in Figure tracks. A track that is also
    // broken keeps applyBrokenFlag()'s tooltip, whichever of the two runs first.
    function applyGeneRequiredFlag() {
      var isCoords = $('input[name="region_mode"]:checked').val() === 'coords';
      $('#tracks-block .rank-list-item').each(function() {
        var key      = $(this).attr('data-rank-id') || '';
        var type     = TRACK_TYPE_MAP[key] || '';
        var flagIt   = isCoords && GENE_REQUIRED_TYPES.indexOf(type) !== -1;
        $(this).toggleClass('track-gene-required', flagIt);
        if (TRACK_BROKEN_MAP[key]) return;
        $(this).attr('title', flagIt ? "Needs 'Gene name' region mode." : '');
      });
    }

    // Every .rank-list-item in the block, Figure tracks included (a broken track
    // can already be placed there), not just the Available list. Never cleared:
    // TRACK_BROKEN_MAP is fixed at page load.
    function applyBrokenFlag() {
      $('[data-available-id]').each(function() {
        $(this).find('.rank-list-item').each(function() {
          var key    = $(this).attr('data-rank-id') || '';
          var broken = TRACK_BROKEN_MAP[key] || false;
          if (broken) {
            $(this).addClass('track-broken');
            $(this).attr('title', "File path is missing or unreadable.");
          }
        });
      });
    }

    // Locks the region_mode radios and the gene box until a sample-count mode and
    // its condition(s) are chosen. Driven by the server's condition_lock message,
    // since whether the conditions are ready depends on inputs this script has
    // no simple DOM check for. The Coordinates panel needs no lock of its own: it
    // can't be reached until region_mode can be clicked.
    function setConditionGateLocked(locked) {
      var $regionModeInputs = $('input[name="region_mode"]');
      $regionModeInputs.prop('disabled', locked);
      $regionModeInputs.closest('.radio-inline').toggleClass('disabled-radio-label', locked);
      var geneSel = document.getElementById('gene_name');
      if (geneSel && geneSel.selectize) {
        if (locked) geneSel.selectize.disable(); else geneSel.selectize.enable();
      }
    }

    // The track picker locks (both Available and Figure tracks, not just
    // visually dimmed) whenever the server-side region_lock message says the
    // mode+conditions+region gate isn't satisfied yet -- see the
    // `region_lock` custom message handler below for how this actually
    // disables dragging via each rank-list's own Sortable instance.
    function setPickerLocked(locked) {
      setupRepeatableTrack();
      addListButtons();
      $('#add_all_tracks, #clear_all_tracks').prop('disabled', locked);
      ['available_tracks', 'figure_tracks'].forEach(function(id) {
        var el = document.getElementById(id);
        if (!el) return;
        var inst = window.Sortable && window.Sortable.get(el);
        if (inst) inst.option('disabled', locked);
        $(el).toggleClass('rank-list-locked', locked);
      });
    }

    // Gene annotation ('genes') can sit in Figure tracks any number of times, every
    // one the same track, while Available always holds exactly one: dragging it
    // out of Available leaves one behind (Sortable's clone pull), and one dropped
    // on Available that already has one is removed (a drop's 'add' event runs
    // before the list reports itself in 'sort'). bucket_list() gives each list
    // its group as a plain name, so the clone pull is set here on the live
    // Sortable instance, once (update_rank_list() keeps the instances).
    function isRepeatableTrack(el) { return el.getAttribute('data-rank-id') === 'genes'; }
    function holdsRepeatableTrack(list, except) {
      return Array.prototype.some.call(list.children, function(c) { return c !== except && isRepeatableTrack(c); });
    }
    function setupRepeatableTrack() {
      var avail = document.getElementById('available_tracks');
      var inst = avail && window.Sortable && Sortable.get(avail);
      if (!inst || inst.ptRepeatable) return;
      inst.ptRepeatable = true;
      inst.option('group', { name: inst.options.group.name,
        pull: function(to, from, dragEl) { return isRepeatableTrack(dragEl) ? 'clone' : true; }, put: true });
      var onAdd = inst.options.onAdd;
      inst.option('onAdd', function(evt) {
        if (isRepeatableTrack(evt.item) && holdsRepeatableTrack(avail, evt.item)) evt.item.parentNode.removeChild(evt.item);
        if (onAdd) onAdd.call(this, evt);
      });
    }

    // Moves tracks of one list, in the order given, to the bottom of the other
    // (Available -> Figure, Figure -> Available), then lets both lists report
    // themselves through sortable's own onSort handlers, as after a drag (and
    // set sortable's empty-list class, which a drag sets on move). A locked
    // picker ignores it, like drag. Gene annotation as a drag does: a copy goes
    // to Figure tracks, and one sent back to an Available that has one is removed.
    // Returns the last track placed, or null.
    function moveTracks(fromId, items) {
      var from = document.getElementById(fromId);
      var to = document.getElementById(fromId === 'available_tracks' ? 'figure_tracks' : 'available_tracks');
      var src = from && window.Sortable && Sortable.get(from), dst = to && window.Sortable && Sortable.get(to);
      if (!src || !dst || src.option('disabled') || dst.option('disabled') || items.length === 0) return null;
      var last = null;
      items.forEach(function(el) {
        var item = el;
        if (isRepeatableTrack(el) && fromId === 'available_tracks') item = el.cloneNode(true);
        else if (isRepeatableTrack(el) && holdsRepeatableTrack(to, el)) item = null;
        if (item) { to.appendChild(item); last = item; } else from.removeChild(el);
      });
      [from, to].forEach(function(el) { Sortable.utils.toggleClass(el, 'rank-list-empty', el.children.length === 0); });
      [src, dst].forEach(function(s) { if (s.options.onSort) s.options.onSort.call(s, {}); });
      return last;
    }

    // Double-click moves one track.
    $(document).on('dblclick', '#available_tracks .rank-list-item, #figure_tracks .rank-list-item', function() {
      if (window.getSelection) window.getSelection().removeAllRanges();
      var placed = moveTracks(this.parentNode.id, [this]);
      if (placed) placed.scrollIntoView({ block: 'nearest' });
    });

    // 'Add all' at the bottom of Available tracks (every track its filter shows,
    // in its order) and 'Clear all' at the bottom of Figure tracks (every track),
    // each centred inside its list's box, below the list: update_rank_list()
    // replaces only the list, so the buttons stay. Same moves as double-click.
    // Disabled until the server's first region_lock message (setPickerLocked()).
    function addListButtons() {
      [['available_tracks', 'add_all_tracks', 'Add all'], ['figure_tracks', 'clear_all_tracks', 'Clear all']].forEach(function(b) {
        var box = document.getElementById('rank-list-' + b[0]);
        if (!box || document.getElementById(b[1])) return;
        var btn = document.createElement('button');
        btn.type = 'button';
        btn.id = b[1];
        btn.className = 'btn btn-default track-list-action';
        btn.textContent = b[2];
        btn.disabled = true;
        box.appendChild(btn);
      });
    }
    $(document).on('click', '#add_all_tracks', function() {
      var shown = Array.prototype.filter.call(document.getElementById('available_tracks').children,
                                              function(el) { return el.style.display !== 'none'; });
      moveTracks('available_tracks', shown);
    });
    $(document).on('click', '#clear_all_tracks', function() {
      moveTracks('figure_tracks', Array.prototype.slice.call(document.getElementById('figure_tracks').children));
    });

    $(document).ready(function() {
      addListButtons();
      setTimeout(function() { setupRepeatableTrack(); applyAllTrackFilters(); applyGeneRequiredFlag(); applyBrokenFlag(); }, 400);
    });
    $(document).on('shiny:inputchanged', function(e) {
      if (e.name === 'available_tracks' || e.name === 'figure_tracks') {
        setTimeout(applyAllTrackFilters, 80);
      }
      if (e.name === 'available_tracks' || e.name === 'figure_tracks' || e.name === 'region_mode') {
        setTimeout(applyGeneRequiredFlag, 80);
        setTimeout(applyBrokenFlag, 80);
      }
    });
    $(document).on('click', '.filter-btn', function() {
      var $block = $(this).closest('[data-available-id]');
      $block.find('.filter-btn').removeClass('active');
      $(this).addClass('active');
      applyTrackFilter($block[0]);
    });

    // A button of comparisons only (Comparisons) exists in two-sample mode only
    // (see the conditionalPanel around it, .pt_filter_buttons()) -- if it was active
    // when the user switches to single/compose, it's now hidden with no way to click
    // back to 'All', so reset the filter for them instead of leaving it stuck.
    $(document).on('change', 'input[name="view_mode"]', function() {
      if ($(this).val() === 'two') return;
      var $block = $('#tracks-block');
      var $cmp   = $block.find('.shiny-panel-conditional .filter-btn');
      if ($cmp.hasClass('active')) {
        $cmp.removeClass('active');
        $block.find('.filter-btn[data-filter="all"]').addClass('active');
        applyTrackFilter($block[0]);
      }
    });

    // Tooltips: init on load and after Shiny re-renders
    $(function() { $('[data-toggle="tooltip"]').tooltip(); });
    $(document).on('shiny:value', function() {
      setTimeout(function() { $('[data-toggle="tooltip"]').tooltip(); }, 250);
    });

    // -- Loading overlay --
    // Injected into <body> directly so position:fixed works as a true viewport overlay
    // (children of the navbar would sit inside the page layout, which breaks fixed positioning)
    $(document).ready(function() {
      $('body').append(
        '<div id="loading-overlay" style="display:none;">' +
        '<div class="plottrellis-spinner"></div>' +
        '<p style="margin-top:12px;font-size:14px;font-weight:600;color:var(--pt-text);">Rendering\u2026</p>' +
        '<button type="button" id="cancel_render_btn" class="btn btn-default">Cancel</button>' +
        '</div>'
      );
      // Cancel stops the render in the background R process; the server then hides the overlay.
      $('#cancel_render_btn').on('click', function() {
        if (!plottrellisRendering) return;
        this.disabled = true;
        this.textContent = 'Cancelling\u2026';
        Shiny.setInputValue('cancel_render', Date.now(), {priority: 'event'});
      });
    });

    var plottrellisRendering = false;

    var ZOOM_BTN_IDS = ['zoom_in_15x','zoom_in_3x','zoom_in_10x',
                        'zoom_out_15x','zoom_out_3x','zoom_out_10x'];

    function showSpinner() {
      plottrellisRendering = true;
      var cancel = document.getElementById('cancel_render_btn');
      if (cancel) { cancel.disabled = false; cancel.textContent = 'Cancel'; }
      document.getElementById('loading-overlay').style.display = 'block';
      ZOOM_BTN_IDS.forEach(function(id) {
        var el = document.getElementById(id);
        if (el) { el.disabled = true; el.classList.add('disabled'); }
      });
    }

    function hideSpinner() {
      if (!plottrellisRendering) return;
      plottrellisRendering = false;
      document.getElementById('loading-overlay').style.display = 'none';
      ZOOM_BTN_IDS.forEach(function(id) {
        var el = document.getElementById(id);
        if (el) { el.disabled = false; el.classList.remove('disabled'); }
      });
    }

    // Capture-phase gate: runs before Shiny's own click binding (which is
    // bound directly on the button), so a click that arrives while a render
    // is already in flight never reaches Shiny at all, it can't queue up
    // and fire a second render after the first one finishes.
    document.addEventListener('click', function(e) {
      // .generate-btn, not a hardcoded #generate id -- Generate and Apply to
      // figure both carry this class, so neither render skips the
      // spinner/double-click gate.
      var gateSelector = '.generate-btn, #' + ZOOM_BTN_IDS.join(', #');
      var el = e.target.closest(gateSelector);
      if (!el) return;
      if (plottrellisRendering) {
        e.stopImmediatePropagation();
        e.preventDefault();
        return;
      }
      showSpinner();
    }, true);

    // Primary signal: Shiny has delivered new SVG or error content to the DOM.
    // Double rAF gives the browser one paint cycle so SVG is visible before we uncover it.
    $(document).on('shiny:value', function(e) {
      if (!plottrellisRendering) return;
      if (e.name === 'svg_display_ui' || e.name === 'render_error_ui') {
        requestAnimationFrame(function() { requestAnimationFrame(hideSpinner); });
      }
    });

    // Fallback: explicit signal from R (fires after on.exit in do_render)
    Shiny.addCustomMessageHandler('hide_spinner', function(x) { hideSpinner(); });

    // Server-side region_lock drives the track
    // picker's lock state (see setPickerLocked() above and the
    // has_mode()/has_conditions_ready()/resolved_region()-driven observe()
    // in server()).
    Shiny.addCustomMessageHandler('region_lock', function(x) { setPickerLocked(x.locked); });

    // Server-side condition_lock drives the
    // region_mode/gene_name lock state (see setConditionGateLocked() above
    // and the has_mode()/has_conditions_ready()-driven observe() in
    // server()).
    Shiny.addCustomMessageHandler('condition_lock', function(x) { setConditionGateLocked(x.locked); });

    // Safety net: an output's error also dismisses the spinner. Not a silent one (a req()
    // that is not met yet sends an empty message): outputs re-run while a figure draws.
    $(document).on('shiny:error', function(e) { if (e.error && e.error.message) hideSpinner(); });

    // Show output section, expand it, scroll to it after first render.
    Shiny.addCustomMessageHandler('show_output', function(x) {
      var sec = document.getElementById('output-section');
      if (sec) sec.style.display = 'block';
      $('#output-body').collapse('show');
      setTimeout(function() {
        if (sec) sec.scrollIntoView({behavior: 'smooth'});
      }, 350);
    });

    // Keep collapse arrow in sync with open/close state
    // A section nested in another (Figure settings in Compose) fires these on the outer
    // section's body too as the event bubbles: act only on the section that toggled, and
    // only on its own arrow.
    function setSectionCollapsed(body, collapsed) {
      var $section = $(body).closest('.plottrellis-section');
      $section.children('.panel-heading').find('.collapse-arrow').text(collapsed ? '\u25B8' : '\u25BE');
      $section.toggleClass('collapsed-panel', collapsed);
    }
    $(document).on('hide.bs.collapse', '.plottrellis-section .collapse', function(e) {
      if (e.target === this) setSectionCollapsed(this, true);
    });
    $(document).on('show.bs.collapse', '.plottrellis-section .collapse', function(e) {
      if (e.target === this) setSectionCollapsed(this, false);
    });

    // Coord text inputs: only fire Shiny update on blur or Enter (not every keystroke)
    $(document).on('shiny:bound', function(e) {
      var id = e.target.id;
      if (id !== 'region_start_text' && id !== 'region_end_text') return;
      var el = e.target;
      $(el).off('.textInputBinding');  // remove default keyup/input reactive listeners
      $(el).on('blur.coordInput', function() {
        Shiny.setInputValue(id, el.value, {priority: 'event'});
      });
      $(el).on('keydown.coordInput', function(evt) {
        if (evt.which === 13) {
          evt.preventDefault();
          Shiny.setInputValue(id, el.value, {priority: 'event'});
        }
      });
    });
  
