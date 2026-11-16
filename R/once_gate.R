# Repeat-click gate: a button with class "pt-once" runs its action once per click, however
# fast it is clicked again. Used by the composer app and the setup wizard; each puts
# .pt_once_js() in its page head and calls .pt_once_server() in its server.
#
# Client: a click on a busy button is dropped before Shiny sees it, and so is every later click
# of a double-click that started on a gated button, wherever it lands (on the button again, or
# on the dialog the first click opened). Any other click marks the button busy (greyed out
# with the "disabled" class; not the disabled property, which would drop keyboard focus) and
# sends `pt_once_click` with its id, after the button's own value (in the same input message,
# or for a button with its own onclick, the next one; the server handles messages in order).
# Server: after the flush that ran the button's action, `pt_once_ready` releases the button,
# so a deliberate later click works again. A download link has no server round trip (the
# file comes over HTTP), so it is held for 2 seconds instead.
.pt_once_js <- function() tags$script(HTML("
  (function() {
    function gated(e) { return e.target && e.target.closest ? e.target.closest('.pt-once') : null; }
    function release(el) {
      delete el.dataset.ptBusy;
      el.classList.remove('disabled');
      el.removeAttribute('aria-disabled');
    }
    var seqGated = false;  // the first click of the current click sequence was on a gated button
    document.addEventListener('click', function(e) {
      var el = gated(e);
      if (e.detail <= 1) seqGated = !!el;
      if ((e.detail > 1 && seqGated) || (el && el.dataset.ptBusy)) {
        e.stopImmediatePropagation();
        e.preventDefault();
      }
    }, true);
    document.addEventListener('click', function(e) {
      var el = gated(e);
      if (!el || el.dataset.ptBusy) return;
      if (el.tagName === 'A') {
        el.dataset.ptBusy = '1';
        el.classList.add('disabled');
        setTimeout(function() { release(el); }, 2000);
        return;
      }
      // Before Shiny is connected the click goes nowhere, so nothing would release the button.
      if (!(window.Shiny && Shiny.shinyapp && Shiny.shinyapp.isConnected())) return;
      el.dataset.ptBusy = '1';
      el.classList.add('disabled');
      el.setAttribute('aria-disabled', 'true');
      Shiny.setInputValue('pt_once_click', el.id, {priority: 'event'});
    });
    Shiny.addCustomMessageHandler('pt_once_ready', function(id) {
      var el = document.getElementById(id);
      if (el && el.dataset.ptBusy) release(el);
    });
  })();
"))

.pt_once_server <- function(input, session) {
  observeEvent(input$pt_once_click, {
    id <- input$pt_once_click
    session$onFlushed(function() session$sendCustomMessage("pt_once_ready", id), once = TRUE)
  })
}
