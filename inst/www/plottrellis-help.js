// Help popovers (R/app-help-popover.R): the circled-i button beside a field label, in the setup
// wizard and the composer's Figure settings alike. Bootstrap 3's own popover (bundled with
// shiny), trigger "manual": hovering shows it, leaving hides it, a click pins it open (so it
// works on touch and the text can be read at leisure), a second click, Escape or a click
// anywhere else unpins it. Handlers are delegated to the document, so an icon Shiny renders
// later needs no re-initialisation; an icon whose output is re-rendered takes its popover down.
(function ($) {
  var PINNED = 'pt-help-pinned';
  var TEMPLATE = '<div class="popover pt-help-popover" role="tooltip"><div class="arrow"></div>' +
                 '<h3 class="popover-title"></h3><div class="popover-content"></div></div>';

  function popover($el) {
    if (!$el.data('bs.popover')) {
      $el.popover({
        trigger: 'manual', placement: 'auto top', container: 'body', html: false,
        template: TEMPLATE,
        content: function () { return $(this).attr('data-pt-help') || ''; }
      });
    }
    return $el;
  }

  function unpin($keep) {
    $('.pt-help.' + PINNED).not($keep || null).removeClass(PINNED).popover('hide');
  }

  $(document)
    .on('mouseenter', '.pt-help', function () { popover($(this)).popover('show'); })
    .on('mouseleave', '.pt-help', function () {
      var $el = $(this);
      if (!$el.hasClass(PINNED)) $el.popover('hide');
    })
    .on('click', '.pt-help', function (e) {
      // A button inside a <label>: stop the click reaching the label's input
      // (a checkbox would toggle) and the document handler below.
      e.preventDefault();
      e.stopPropagation();
      var $el = popover($(this));
      if ($el.hasClass(PINNED)) {
        $el.removeClass(PINNED).popover('hide');
      } else {
        unpin($el);
        $el.addClass(PINNED).popover('show');
      }
    })
    .on('keydown', function (e) { if (e.key === 'Escape') unpin(); })
    // A modal closing takes its icons out of the page: unpin first, or a pinned popover
    // would stay behind with no icon left to close it.
    .on('hide.bs.modal', function () {
      $('.pt-help').each(function () { if ($(this).data('bs.popover')) $(this).removeClass(PINNED).popover('hide'); });
    })
    .on('click', function (e) {
      if (!$(e.target).closest('.pt-help-popover').length) unpin();
    })
    // Shiny is about to replace this output's contents: drop the popovers of the icons in it,
    // or a pinned one would stay on the page with no icon left to close it.
    .on('shiny:value', function (e) {
      $(e.target).find('.pt-help').each(function () {
        var $el = $(this);
        if ($el.data('bs.popover')) $el.removeClass(PINNED).popover('destroy');
      });
    });
})(jQuery);
