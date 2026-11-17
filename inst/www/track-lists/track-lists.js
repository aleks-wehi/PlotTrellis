// The composer's two track lists (Available, Figure) on SortableJS: the Shiny side of
// bucket_list(), add_rank_list() and update_rank_list() in R/track-lists.R. Adapted from the R
// package sortable 0.6.0 (MIT; see LICENSE), so the lists behave as they did with it: the same
// Sortable options and callbacks, the same inputs ('<list>:sortablejs.rank_list' and
// '<group>:sortablejs.bucket_list') and the same update message.

// The update binding: update_rank_list() replaces a list's items (and title) and reports the
// list's new ids as input '<list>' (sortable's rank-list binding, unchanged).
var ptRankListBinding = new Shiny.InputBinding();
$.extend(ptRankListBinding, {
  find: function(scope) { return $(scope).find(".rank-list-container"); },
  initialize: function(el) { },
  getValue: function(el) { },
  setValue: function(el, data) {
    if (data.text) {
      $(el).find(".rank-list-title").text(data.text);
    }
    if (data.labels) {
      const short_id = data.id.replace(/^rank-list-/, "");
      $('#' + short_id).html(data.labels);
      const label_ids = $('#' + short_id).children().map((idx, child) => {
        return $(child).attr("data-rank-id") || $.trim(child.innerHTML);
      }).toArray();
      Shiny.setInputValue(short_id, label_ids, {priority: 'event'});
    }
  },
  subscribe: function(el, callback) { },
  receiveMessage: function(el, data) { this.setValue(el, data); }
});
Shiny.inputBindings.register(ptRankListBinding);

// The bucket's own binding (header text only; bucket_list() is given none).
var ptBucketListBinding = new Shiny.InputBinding();
$.extend(ptBucketListBinding, {
  find: function(scope) { return $(scope).find(".bucket-list-container"); },
  initialize: function(el) { },
  getValue: function(el) { },
  setValue: function(el, data) {
    if (data.header) {
      $(el).find(".bucket-list-header").text(data.header);
    }
  },
  subscribe: function(el, callback) { },
  receiveMessage: function(el, data) { this.setValue(el, data); }
});
Shiny.inputBindings.register(ptBucketListBinding);

window.ptTrackLists = (function() {
  // An item's id: its data-rank-id, else its text.
  function childIdOrText(child) {
    return $(child).attr('data-rank-id') || $.trim(child.innerText);
  }
  // Calls each function with this call's `this` and arguments; one that throws is logged and
  // the rest still run.
  function chain(fns) {
    return function() {
      var self = this, args = arguments;
      fns.forEach(function(fn) {
        try {
          fn.apply(self, args);
        } catch (e) {
          if (window.console && window.console.error) window.console.error(e);
        }
      });
    };
  }
  // Every list of the group, as input '<group>:sortablejs.bucket_list' ({<input id>: ids}), once
  // Shiny is up (sortable sent it at once, and on page load that could come before Shiny's
  // set-up and fail).
  function captureBucket(group, cssIds, inputIds) {
    return function(evt) {
      if (typeof Shiny == "undefined") {
        return;
      }
      var ret = {};
      $.map(cssIds, function(cssId, i) {
        var item = $('#' + cssId).get(0);
        ret[inputIds[i]] = item && item.children ? $.map(item.children, childIdOrText) : undefined;
      });
      Shiny.initializedPromise.then(() => {
        Shiny.setInputValue(group + ":sortablejs.bucket_list", ret);
      });
    };
  }
  // This list's ids, as input '<input id>:sortablejs.rank_list', once Shiny is up.
  function captureList(inputId) {
    return function(evt) {
      if (typeof Shiny !== "undefined") {
        Shiny.initializedPromise.then(() => {
          Shiny.setInputValue(inputId + ":sortablejs.rank_list", $.map(this.el.children, childIdOrText));
        });
      }
    };
  }
  // 'rank-list-empty' on each list with no items, after the current event.
  function setEmptyClass(cssIds) {
    return function(evt) {
      setTimeout(function() {
        cssIds.map(function(id) {
          var el = window.document.getElementById(id);
          if (el) {
            Sortable.utils.toggleClass(el, 'rank-list-empty', el.children.length == 0);
          }
        });
      }, 0);
    };
  }
  // One Sortable per list, all in one drag group; each reports itself and the bucket on load
  // and after every sort.
  function createNow(group, cssIds, inputIds) {
    var bucket = captureBucket(group, cssIds, inputIds);
    cssIds.forEach(function(cssId, i) {
      var el = document.getElementById(cssId);
      if (!el) return;
      var inst = Sortable.create(el, {
        animation: 150,
        emptyInsertThreshold: 50 / 4,
        group: group,
        onLoad: chain([bucket, captureList(inputIds[i]), setEmptyClass([cssId])]),
        onSort: chain([bucket, captureList(inputIds[i])]),
        onMove: setEmptyClass(cssIds)
      });
      setTimeout(function() {
        var evt = document.createEvent('Event');
        evt.initEvent("onLoad", true, true);
        inst.options.onLoad.apply(inst, [evt]);
      }, 0);
    });
  }
  // Created when the page's DOM is ready (as soon as called, if it already is).
  function create(group, cssIds, inputIds) {
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", function() { createNow(group, cssIds, inputIds); });
    } else {
      createNow(group, cssIds, inputIds);
    }
  }
  return { create: create };
})();
