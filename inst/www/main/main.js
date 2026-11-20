// PlotTrellis main screen (the page of R/app-screen-ui.R). Its markup is main.html, its styles main.css.
// The script runs once, after the markup; it keeps the whole state of the page in S.
/* ------------------------------------------------------------------ the server
   Every question the page asks goes through `api`. In the app the answers come from the R session
   (one Shiny input, "pt_rpc", carries the questions; the "pt_rpc_reply" message brings each answer back).
   A test page may replace `api` with its own object before PT.start() (see solutions/review/W78/harness). */
const rpcWait = {}; let rpcN = 0, rpcOut = [], rpcT = 0;
function rpc(op, args) { return new Promise(done => { if (rpcGone) return done({ error: GONE }); const id = ++rpcN; rpcWait[id] = done; rpcOut.push({ id, op, args: args || {} });
  // questions asked in the same moment travel together: one input value per send, so none can overwrite another
  if (!rpcT) rpcT = setTimeout(() => { const batch = rpcOut; rpcOut = []; rpcT = 0; Shiny.setInputValue('pt_rpc', JSON.stringify(batch), { priority: 'event' }); }, 0); }); }
let api = {
  wait: 120,   // ms to wait after a keystroke before asking for suggestions
  project() { return rpc('project'); },
  search(text) { return rpc('search', { q: text }); },
  locus(text) { return rpc('locus', { q: text }); },
  nearby(chr, s, e) { return rpc('nearby', { chr, s: Math.round(s), e: Math.round(e) }); },
  render(req) { return rpc('render', req); },
  raster(job, format) { return rpc('raster', { job, format }); },
  cancel() { return rpc('cancel'); },
  geneinfo(gene) { return rpc('geneinfo', { gene }); }
};

/* ------------------------------------------------------------------ state */
const GONE = 'The connection to R was lost. Reload the page.'; let rpcGone = false;
const byId = id => document.getElementById(id);   // not named $: that is jQuery's, which Shiny's own scripts look up
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const fmt = n => Math.round(n).toLocaleString('en-US');
const fmtSpan = n => n >= 1e6 ? (n / 1e6).toFixed(n >= 1e7 ? 0 : 1) + ' Mb' : n >= 1e3 ? Math.round(n / 1e3) + ' kb' : Math.round(n) + ' bp';
const q = new URLSearchParams(location.search);
let P, S, TYPE = {}, REG = [], addable = [];
let SETKEY = '', setSeq = 0;   // the Figure settings as the server last described them: part of what makes a figure the same figure
const SAVED = {};   // each project keeps its own state, so switching project and back loses nothing
const SPAN_MIN = 1000, SPAN_MAX = 5000000;

async function reset(which) {
  if (P && S) SAVED[P.id] = S;
  document.body.classList.add('loading'); byId('genome').hidden = false; byId('genome').textContent = 'Loading project\u2026'; byId('gene').disabled = true;
  P = await api.project(which);
  if (!P || P.error) { byId('genome').textContent = 'The project could not be loaded: ' + (P && P.error || 'no answer'); return; }
  document.body.classList.remove('loading'); byId('gene').disabled = false;
  byId('chr').innerHTML = Object.keys(P.chrLen || {}).map(c => `<option>${esc(c)}</option>`).join(''); byId('ptVer').textContent = P.version ? 'v' + P.version : '';
  { const fx = P.types.filter(T => T.kind === 'fixed').sort((a, b) => annAt(a.id) - annAt(b.id)); let i = 0; P.types = P.types.map(T => T.kind === 'fixed' ? fx[i++] : T); }
  for (const T of P.types) T.name = T.name.replace(/\s*\(stacked\)/, '');
  for (const T of P.types) T.name = T.name.replace(/\s*mirrored(?= diff)/, '');   // "Junctions diff", "Junctions diff (TPM-scaled)"
  TYPE = Object.fromEntries(P.types.map(t => [t.id, t]));
  S = SAVED[P.id] || { region: false, locus: null, view: null, err: null, nearby: [], picked: [], gap: null, fig: [], drawn: null, busy: false, fail: null, tab: 'figure', info: null,
    exp: {}, hist: [], recent: [], sugg: [], suggAt: 0, lastGene: '', set: { key: SETKEY } };
  S.busy = false; S.err = null; S.drawn = S.hist.find(h => h.id === sig()) || null;
  byId('gene').value = S.locus && S.locus.gene || ''; byId('range').value = ''; byId('gene').placeholder = P.suggested[0] || 'Gene'; mirror();
  byId('genome').hidden = true;
  document.body.classList.toggle('many', P.samples.length > 8);
  ['roles', 'samples', 'tracks', 'right', 'hist'].forEach(id => { byId(id).hidden = true; }); byId('stage').dataset.key = '';
  setRegion(S.region); render(); if (S.locus) loadNearby();
  byId('gene').focus();   // only now: its focus handler reads the new project's state
  if (document.activeElement === byId('gene')) updateSugg();   // first open: the suggested genes are already listed
}
const sampleOf = name => P.samples.find(s => s.name === name);
const pair = () => S.picked.length === 2;   // picked[0] is the comparison, picked[1] the reference

// A track is one line of the figure. Key: "type|sample" for a per-sample track, "type" otherwise.
function tk(key) {
  const base = key.replace(/~\d+$/, ''), i = base.indexOf('|'), T = TYPE[i < 0 ? base : base.slice(0, i)], sample = i < 0 ? null : base.slice(i + 1);
  const geneOk = !T.gene || !!(S.locus && S.locus.gene);
  if (T.kind === 'sample') return { key, T, cat: T.cat, name: T.name, sub: sample, col: sampleOf(sample).colour, ok: geneOk && S.picked.includes(sample) && T.samples.includes(sample) };
  if (T.kind === 'comparison') return { key, T, cat: 'cmp', name: T.name, sub: pair() ? pairShort() : '', full: pair() ? S.picked[0] + ' vs ' + S.picked[1] : '', ok: geneOk && pair() && S.picked.every(n => T.samples.includes(n)) };
  return { key, T, cat: 'ann', name: T.name, ok: geneOk };
}
// The two picked names without the words they share at the start and at the end ("naive vs memory"); the sentence at the top has them in full.
function pairShort() { const a = S.picked[0].split(' '), b = S.picked[1].split(' '); let i = 0, j = 0;
  while (i < a.length - 1 && i < b.length - 1 && a[i] === b[i]) i++;
  while (j < a.length - 1 - i && j < b.length - 1 - i && a[a.length - 1 - j] === b[b.length - 1 - j]) j++;
  return a.slice(i, a.length - j).join(' ') + ' vs ' + b.slice(i, b.length - j).join(' '); }
// A picked sample's name without the words all picked names start with and end with ("naive" and "memory" for
// "CD8 T cells naive resting" and "CD8 T cells memory resting"). Never an invented code; a name alone stays whole.
function sn(name) { const p = S.picked; if (p.length < 2 || !p.includes(name)) return name;
  const ws = p.map(x => x.split(' ')), min = Math.min(...ws.map(w => w.length)); let i = 0, j = 0;
  while (i < min - 1 && ws.every(w => w[i] === ws[0][i])) i++;
  while (j < min - 1 - i && ws.every(w => w[w.length - 1 - j] === ws[0][ws[0].length - 1 - j])) j++;
  const w = name.split(' '); return w.slice(i, w.length - j).join(' '); }
// A name too long for its place is cut in the middle: similar samples differ at the end, so the end stays.
function midCut(t, n) { if (t.length <= n) return t; const c = midCut1(t, n);
  return S.picked.some(o => { const x = sn(o); return x !== t && x.length > n && midCut1(x, n) === c; }) ? t : c; }   // a cut that would make two names read alike is not made
function midCut1(t, n) { const sp = t.indexOf(' '), head = sp > 0 && sp <= n / 2 ? t.slice(0, sp) : t.slice(0, Math.floor(n / 3));
  let tail = t.slice(-(n - head.length - 1)); const cut = tail.indexOf(' '); if (cut >= 0 && cut < tail.length / 2) tail = tail.slice(cut + 1);
  return head + ' \u2026 ' + tail; }
// Per-sample tracks of one kind run in the order of the sample list; for a pair, the comparison comes first.
const rank = n => pair() && S.picked.includes(n) ? S.picked.indexOf(n) - 2 : P.samples.findIndex(x => x.name === n);
const fullName = t => t.name + (t.full || t.sub ? ' \u00B7 ' + (t.full || t.sub) : '');
// Every track that can be drawn with the current choice, in the app's usual order.
function pool() {
  const out = [];
  for (const T of P.types) { if (T.kind === 'sample') for (const n of S.picked) out.push(T.id + '|' + n); else out.push(T.id); }
  return out.filter(k => tk(k).ok);
}
const CATS = { cmp: 'Comparisons', expr: 'Expression', ann: 'Annotations' };
// Annotation tracks, most used first, each under its id in the checks' made-up data and its id in the app.
// The first three tracks (five ids) carry a star.
const ANN_ORDER = ['genes', 'tx_ann', 'transcripts_annotated', 'protein', 'protein_domains', 'tx_all', 'transcripts_all', 'tx', 'transcripts'];
const STARRED = ANN_ORDER.slice(0, 5), annAt = id => (ANN_ORDER.indexOf(id) + 1) || 99;
const THINT = '<span class="chunk empty"><svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"><path d="M14 8H3M7 4 3 8l4 4"/></svg>select tracks</span>';
const cats = pl => { const cs = [...new Set(pl.map(k => tk(k).cat))], rank = c => c === 'cmp' ? 1 : c === 'expr' ? 2 : c === 'ann' ? 3 : 0; return cs.sort((a, b) => rank(a) - rank(b)); };
const shown = () => S.fig.filter(k => tk(k).ok);
const inFig = k => S.fig.includes(k);
const reg = keys => REG.push(keys) - 1;
// A track added by a click goes next to its own kind: after the last track of the same type (a difference track: after its coverage tracks), else at the end.
function place(k) { const t = tk(k), T = t.T; let at = -1, first = -1; S.fig.forEach((x, i) => { const o = tk(x); if (o.T !== T) return; if (first < 0) first = i;
    if (T.kind !== 'sample' || rank(o.sub) <= rank(t.sub)) at = i; });
  if (at < 0 && first >= 0) return void S.fig.splice(first, 0, k);
  if (at < 0 && T.kind === 'sample') { const d = S.fig.findIndex(x => { const o = tk(x).T; return o.kind === 'comparison' && (o.data || o.id.replace(/_diff.*$/, '')) === T.id; }); if (d >= 0) return void S.fig.splice(d, 0, k); }
  if (at < 0 && T.kind === 'comparison') { const cov = T.data || T.id.replace(/_diff.*$/, ''); S.fig.forEach((x, i) => { if (x.split('|')[0] === cov) at = i; }); }
  if (at >= 0) return void S.fig.splice(at + 1, 0, k);
  let end = S.fig.length; if (T.kind !== 'fixed') while (end > 0 && tk(S.fig[end - 1]).T.kind === 'fixed') end--;
  S.fig.splice(end, 0, k); }
const byHand = k => { const T = tk(k).T; if (T.kind === 'sample' && !(S.hand || []).includes(T.id)) S.hand = [...(S.hand || []), T.id]; };
function reranked(mutate) { const by = (a, b) => rank(tk(a).sub) - rank(tk(b).sub);
  const runs = P.types.filter(T => T.kind === 'sample').map(T => { const at = S.fig.map((x, i) => tk(x).T === T ? i : -1).filter(i => i >= 0), ks = at.map(i => S.fig[i]); return { at, ks, tidy: !(S.hand || []).includes(T.id) && ks.join('\n') === [...ks].sort(by).join('\n') }; });
  mutate(); for (const r of runs) if (r.tidy) { const ks = [...r.ks].sort(by); r.at.forEach((i, j) => { S.fig[i] = ks[j]; }); } }
function toggle(keys) { changed(() => { if (keys.every(inFig)) { S.fig = S.fig.filter(k => !keys.includes(k)); } else for (const k of keys) if (!inFig(k)) place(k); }); }
const pressed = keys => keys.every(inFig) ? 'true' : keys.some(inFig) ? 'mixed' : 'false';

function reveal(id) { const el = byId(id); if (el.hidden) { el.hidden = false; el.classList.add('pop'); } }
// What would be drawn now. The same choice drawn earlier this session comes straight back from the cache.
const sig = () => S.locus && S.view ? JSON.stringify([P.id, S.locus.gene, S.locus.chr, S.view.s, S.view.e, S.picked, shown().map(k => k.replace(/~\d+$/, '')), S.set]) : '';
function changed(mutate) { const b4 = S.fig.slice(); if (mutate) mutate();
  const rm = new Set(S.rm || []); for (const k of b4) if (!S.fig.includes(k)) rm.add(k); for (const k of S.fig) rm.delete(k); S.rm = S.fig.length ? [...rm] : [];
  const id = sig(); S.drawn = S.hist.find(h => h.id === id) || null; S.fail = null; render(); }
function toast(msg) { const t = byId('toast'); t.textContent = msg; t.hidden = false; clearTimeout(toast.t); toast.t = setTimeout(() => { t.hidden = true; }, 3200); }

const dot = c => `<span class="dot" style="background:${/^[#\w(),.%\s-]+$/.test(c || '') ? c : 'transparent'}"></span>`;   // a colour, never markup
// The In figure rows name the samples in full; a comparison row gives both, each with its swatch.
const label = t => t.T.kind === 'comparison' && t.sub
  ? `<b>${esc(t.name)}</b><span class="sub"><span class="l2 cmp">${S.picked.slice(0, 2).map(x => dot(sampleOf(x).colour) + esc(x)).join(' vs ')}</span></span>`
  : `<b${t.sub ? '' : ' class="solo"'}>${esc(t.name)}</b>${t.sub ? `<span class="sub">${t.col ? dot(t.col) : ''}<span class="l2">${esc(t.sub)}</span></span>` : ''}`;
const tchip = (keys, html, cls = '') => `<button class="tchip ${cls}" aria-pressed="${pressed(keys)}" data-g="${reg(keys)}">${html}</button>`;

/* ------------------------------------------------------------------ the Add half
   renderAdd() fills #addList and returns the tracks "Add all" should add.
   Controls carry data-g="<index into REG>" (a list of track keys to switch on or off). */
// Comparisons and annotations as on/off buttons.
function singles(pl, cat) {
  const ks = pl.filter(k => tk(k).cat === cat); if (!ks.length) return '';
  if (cat === 'ann') { const star = k => STARRED.includes(k), chip = k => tchip([k], (star(k) ? '<span class="star" title="Most used" aria-label="Most used">\u2605</span>' : '') + esc(tk(k).name));
    return `<div class="cap">${CATS[cat]}</div><div class="wrap">${ks.filter(star).map(chip).join('')}</div><div class="wrap">${ks.filter(k => !star(k)).map(chip).join('')}</div>`; }
  return `<div class="cap">${CATS[cat]}</div><div class="wrap">${ks.map(k => tchip([k], esc(tk(k).name))).join('')}</div>`;
}
// One button per kind of data puts it in for every picked sample; the arrow opens the samples.
function renderAdd(pl) {
  const kind = cat => { const ks = pl.filter(k => tk(k).cat === cat); if (!ks.length) return '';
    return `<div class="kind">${tchip(ks, `<span class="trunc">${esc(cat === 'expr' ? tk(ks[0]).name : cat)}</span><span class="n">${ks.filter(inFig).length}/${ks.length}</span>`)}${ks.length > 1 ? `<button class="ib" data-x="${esc(cat)}" aria-label="Samples" aria-expanded="${!!S.exp[cat]}">${S.exp[cat] ? '&#9662;' : '&#9656;'}</button>` : ''}</div>`
      + (S.exp[cat] && ks.length > 1 ? `<div class="wrap per">${ks.map(k => tchip([k], dot(tk(k).col) + `<span class="l2" title="${esc(tk(k).sub)}">${esc(sn(tk(k).sub))}</span>`)).join('')}</div>` : ''); };
  const data = cats(pl).filter(c => !CATS[c]).map(kind).join('');
  // comparisons kept in the figure that cannot be drawn with this choice of samples: still here, marked, a click takes one out
  const kept = S.fig.filter(k => tk(k).T.kind === 'comparison' && !tk(k).ok), two = S.picked.length === 2;
  const keptHtml = kept.length ? `<div class="cap">${pl.some(k => tk(k).cat === 'cmp') ? '' : CATS.cmp}<span class="only">${two ? 'not for this pair' : 'only for 2 samples'}</span></div><div class="wrap kept">${kept.map(k => tchip([k], esc(tk(k).name))).join('')}</div>` : '';
  byId('addList').innerHTML = data + kind('expr') + singles(pl, 'cmp') + keptHtml + singles(pl, 'ann');
  return pl.filter(k => !inFig(k));
}
/* ------------------------------------------------------------------ drawing the page */
// After a re-render the control that had the focus is found again (or whatever now sits in its place), so the keyboard is never thrown back to the top of the page.
function focusKey() { const el = document.activeElement; if (!el || el === document.body || !el.closest) return null;
  const cont = el.closest('#chips, #inList, #addList, #where, #hist, #roles, #stage'); if (!cont) return null;
  const d = el.dataset, q1 = v => '"' + String(v).replace(/["\\]/g, '\\$&') + '"';
  return { cont: cont.id, idx: [...cont.querySelectorAll('button, [tabindex]')].indexOf(el),
    sel: el.id ? '#' + el.id : ['s', 'rm', 'dup', 'k', 'x', 'open', 'zoom', 'pan', 'fit', 'gene'].filter(a => d[a] != null).map(a => `[data-${a}=${q1(d[a])}]`)[0] || null,
    g: d.g != null && REG[+d.g] ? REG[+d.g].join('\n') + '|' + el.className : null }; }
function restoreFocus(fk) { const now = document.activeElement; if (!fk || (now && now !== document.body && document.contains(now))) return;
  const cont = byId(fk.cont); if (!cont || cont.closest('[hidden]')) return; const all = [...cont.querySelectorAll('button, [tabindex]')];
  const el = (fk.sel && cont.querySelector(fk.sel)) || (fk.g && all.find(b => b.dataset.g != null && REG[+b.dataset.g] && REG[+b.dataset.g].join('\n') + '|' + b.className === fk.g)) || all[Math.min(fk.idx, all.length - 1)];
  if (el && !el.disabled) el.focus({ preventScroll: true }); }
/* ------------------------------------------------------------------ saved figures
   S.hist holds the figures drawn this session, newest first; S.drawn is the one on screen (or null).
   An entry: { id, no, label, sub, at, title, svg, html, script, stem, thumb, state }.
   reopen(i) puts figure i back on screen with its choice; forget(i) deletes it; download(name, text, type) saves a file.
   #hist is the row under the canvas: the filmstrip fills it; its flick overlay is added to the body. */
function forget(i) { const gone = S.hist.splice(i, 1)[0]; if (gone === S.drawn) S.drawn = null; render(); }
// The filmstrip. Click a card: its figure opens (pointing at a card shows nothing). Hold Space: the canvas flicks to the figure shown
// before this one (dashed edge). A card's one button deletes it; deleting the one on screen leaves Generate on the canvas.
// While a figure is being drawn no saved figure can be opened, stepped to or flicked to: the cards are disabled until it ends.
// Cards of one gene stand together under the gene's name. The fold button shrinks the strip to a bar of dots (remembered per project).
let gsOn = false, gsFlick = false, gsOther = null, gsLast = null, gsFoldMem = {}, gsCur = null;
const gsI = {
  x: '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"><path d="M4 4l8 8M12 4l-8 8"/></svg>',
  down: '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M3.5 6 8 10.5 12.5 6"/></svg>',
  up: '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M3.5 10 8 5.5l4.5 4.5"/></svg>',
  eye: '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.6"><path d="M1.5 8s2.4-4.5 6.5-4.5S14.5 8 14.5 8 12.1 12.5 8 12.5 1.5 8 1.5 8z"/><circle cx="8" cy="8" r="2"/></svg>' };
const gsById = id => S.hist.find(h => h.id === id) || null;
const gsSvg = h => h.html.replace(/class="fig pop"/, 'class="gs-fig"');
// A card's picture is the figure drawn once onto a small canvas, kept as a PNG. Not the SVG itself: as a picture each SVG is a
// whole document of its own, and its text went through the strip's markup on every redraw, so 40 big figures (Hi-C) held
// millions of nodes and made every click slow. The card's picture is set in place when it is ready; no redraw.
function gsThumb(h) { const u = URL.createObjectURL(new Blob([h.svg], { type: 'image/svg+xml' })), img = new Image();
  const done = src => { URL.revokeObjectURL(u); h.thumb = src; const el = document.querySelector(`#hist .gs-c[data-gsid="${CSS.escape(h.id)}"] img`); if (el) el.src = src; };
  img.onload = () => { try { const c = document.createElement('canvas'); c.width = 276; c.height = Math.max(1, Math.round(276 * (img.naturalWidth ? img.naturalHeight / img.naturalWidth : 0.6)));
      const g = c.getContext('2d'); g.fillStyle = '#FFFFFF'; g.fillRect(0, 0, c.width, c.height); g.drawImage(img, 0, 0, c.width, c.height); done(c.toDataURL('image/png')); }
    catch (e) { done('data:image/svg+xml;charset=utf-8,' + encodeURIComponent(h.svg)); } };   // a canvas that cannot be read: the SVG as before
  img.onerror = () => done(''); img.src = u; }
function gsFolded() { if (!(P.id in gsFoldMem)) { try { gsFoldMem[P.id] = localStorage.getItem('gsFold:' + P.id) === '1'; } catch (e) { gsFoldMem[P.id] = false; } } return gsFoldMem[P.id]; }
function gsSetFold(v) { gsFoldMem[P.id] = v; try { localStorage.setItem('gsFold:' + P.id, v ? '1' : '0'); } catch (e) {} }
// Oldest on the left; one group per gene (or region), groups in the order they first appeared.
function gsGroups() { const out = [], at = {}; for (const h of [...S.hist].reverse()) { if (!(h.label in at)) { at[h.label] = out.length; out.push({ label: h.label, hs: [] }); } out[at[h.label]].hs.push(h); } return out; }
function gsOv(id) { let el = byId(id); if (!el) { el = document.createElement('div'); el.id = id; el.className = 'gs-ov'; el.hidden = true; document.body.appendChild(el); } return el; }
function gsPlace(el) { const b = byId('stage').getBoundingClientRect(); Object.assign(el.style, { left: b.left + 'px', top: b.top + 'px', width: b.width + 'px', height: b.height + 'px' }); }
// What the canvas shows for a moment: a flick (Space held).
function gsShowTemp() { const el = gsOv('gs-peek'), id = gsFlick ? gsOther || (gsLast && gsLast.id) : null, h = id && gsById(id);
  if (!gsOn || !h || h === S.drawn || S.busy || byId('right').hidden) { el.hidden = true; el.dataset.id = ''; return; }
  if (el.dataset.id !== h.id) { el.dataset.id = h.id; el.innerHTML = `<span class="gs-tag">${gsI.eye}${esc(h.label)} \u00B7 ${esc(h.sub)}</span>` + gsSvg(h); }
  const fitH = FITH && S.tab === 'figure', f = el.querySelector('svg.gs-fig'); el.classList.toggle('fit-h', fitH);
  if (f) { if (fitH) f.setAttribute('preserveAspectRatio', 'xMidYMin meet'); else f.removeAttribute('preserveAspectRatio'); }
  gsPlace(el); el.hidden = false; }
function gsCard(h) { const cur = h === S.drawn, i = S.hist.indexOf(h);
  return `<div class="gs-c${cur ? ' cur' : ''}" data-gsid="${esc(h.id)}"><button class="gs-th" data-open="${i}" title="${esc(h.title)}" aria-label="Open ${esc(h.title)}"${S.busy ? ' disabled' : ''}><img alt="" src="${h.thumb}"></button>`
    + `<button class="gs-x" data-hx="${i}" aria-label="Delete ${esc(h.title)}" title="Delete">${gsI.x}</button>`
    + `<span class="gs-m"><span class="trunc">${esc((h.at ? h.at + ' · ' : '') + h.sub)}</span></span></div>`; }
function gsDraw() { const hb = byId('hist'), fold = gsFolded(), n = S.hist.length;
  document.body.classList.toggle('gs-fold', fold);
  if (fold) { hb.innerHTML = `<div class="gs-bar"><button class="gs-ib" data-gsfold aria-label="Show the saved figures" title="Show">${gsI.up}</button><span class="gs-n" title="${n} saved figure${n === 1 ? '' : 's'}">${n}</span></div>`; return; }
  if (!n) { hb.innerHTML = ''; return; }
  hb.innerHTML = `<div class="gs-h"><button class="gs-ib" data-gsfold aria-label="Fold the saved figures" title="Fold">${gsI.down}</button><span class="gs-n" title="Figures this session. Click one to open, hold Space to flick.">${n}</span></div>`
    + gsGroups().map(g => `<div class="gs-g"><span class="gs-gl trunc" title="${esc(g.label)}">${esc(g.label)}</span><div class="gs-row">${g.hs.map(gsCard).join('')}</div></div>`).join('');
  // keep the newest card, or the one just opened, in view
  const want = S.hist[0] && hb.dataset.gsTop !== S.hist[0].id ? S.hist[0] : S.drawn !== gsCur ? S.drawn : null; if (S.hist[0]) hb.dataset.gsTop = S.hist[0].id; gsCur = S.drawn;
  const c = want && hb.querySelector(`.gs-c[data-gsid="${CSS.escape(want.id)}"]`);
  if (c) { const r = hb.getBoundingClientRect(), b = c.getBoundingClientRect(), left = r.left + 36; if (b.left < left) hb.scrollLeft -= left - b.left; else if (b.right > r.right - 4) hb.scrollLeft += b.right - r.right + 4; } }
byId('hist').addEventListener('click', e => { if (gsOn && e.target.closest('[data-gsfold]')) { gsSetFold(!gsFolded()); render(); } });
// Space held: flick to the other figure; released: back. Only where Space would not press a button or type.
const gsSpaceOk = e => { const a = document.activeElement, t = e.target; return gsOn && S.hist.length && !S.busy && !document.querySelector('dialog[open], .modal.in') && !(t.matches && t.matches('input, select, textarea'))
  && !(t.closest && t.closest('#ptSettings, .modal, summary'))
  && (!a || a === document.body || !!a.closest('#hist, #stage') || !a.matches('button, [role="switch"], [tabindex]')); };
window.addEventListener('keydown', e => { if (e.key === ' ' && gsSpaceOk(e)) { e.preventDefault(); e.stopPropagation(); if (!e.repeat && !gsFlick) { gsFlick = true; gsShowTemp(); } } }, true);
window.addEventListener('keyup', e => { if (e.key === ' ' && gsFlick) { e.preventDefault(); gsFlick = false; gsShowTemp(); } }, true);
window.addEventListener('blur', () => { if (gsFlick) { gsFlick = false; gsShowTemp(); } });
window.addEventListener('resize', () => { if (gsOn) gsShowTemp(); });
// Runs on every redraw.
function gsRender() {
  gsOn = true; document.body.classList.add('gs-on');
  if (S.drawn && S.drawn !== gsLast) { if (gsLast && gsLast.id !== S.drawn.id && gsById(gsLast.id)) gsOther = gsLast.id; gsLast = S.drawn; }
  gsDraw(); gsShowTemp();
}
// "Enter draws": shown in the box that has the keyboard, when its content is what is on screen and there is a figure to draw.
function enterHint() { if (!S) return; const box = S.region ? byId('range') : byId('gene'), same = S.locus && (S.region || box.value.trim().toLowerCase() === String(S.locus.gene || '').toLowerCase());
  byId('ehint').hidden = !(document.activeElement === box && same && S.picked.length && shown().length && !S.drawn && !S.busy && byId('sugg').hidden); }
let drag = null;
function render() {
  const fk = focusKey(), had = document.activeElement;
  REG = [];
  // samples, and the sentence that says what is picked
  const n = P.samples.length, p = S.picked;
  byId('chips').className = 'chips ' + (n <= 4 ? 'few' : n <= 8 ? 'some' : 'lots');
  byId('chips').innerHTML = P.samples.map((s, i) => `<button class="chip" aria-pressed="${p.includes(s.name)}" data-s="${esc(s.name)}"${n <= 9 ? ` aria-keyshortcuts="${i + 1}"` : ''}>${dot(s.colour)}${esc(s.name)}</button>`).join('');
  const allOn = p.length === n; byId('allBtn').setAttribute('aria-pressed', String(allOn)); byId('allBtn').title = allOn ? 'Clear the samples' : 'Pick all samples'; byId('allBtn').setAttribute('aria-label', byId('allBtn').title); byId('allKey').textContent = allOn ? '0' : 'A';
  const chunk = (name, role) => `<span class="chunk" title="${esc(name)}${role ? ' (' + role + ')' : ''}">${dot(sampleOf(name).colour)}<span class="trunc">${esc(name)}</span></span>`;
  byId('roles').innerHTML = '<span class="word">in</span>' + (!p.length ? '<span class="chunk empty">select samples<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"><path d="M2 8h11M9 4l4 4-4 4"/></svg></span>' : p.length === 1 ? chunk(p[0])
    : p.length === 2 ? `${chunk(p[0], 'comparison')}<span class="word">relative to</span>${chunk(p[1], 'reference')}<button class="swap" id="swap" aria-label="Switch which sample is the reference" title="Switch which sample is the reference"><svg width="16" height="16" viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"><path d="M2 5h11M10 2l3 3-3 3M14 11H3M6 8l-3 3 3 3"/></svg></button>`
    : `<span class="chunk multi" title="${esc(p.join(', '))}">${(p.length > 5 ? [] : p).map(x => dot(sampleOf(x).colour)).join('')}<span>${p.length} samples</span></span>`)
    + (p.length ? '<span class="unp"><button class="unpick" id="unpick" aria-label="Clear the samples" title="Clear the samples">&#10005;</button></span>' : '');
  const box = S.region ? byId('range') : byId('gene'); box.setAttribute('aria-invalid', String(!!S.err)); box.title = S.err || '';
  byId('clear').hidden = !box.value;
  if (S.locus) { reveal('roles'); reveal('samples'); }
  if (S.locus && p.length) reveal('tracks');
  // tracks in the figure, one row per track, top first
  const inF = shown(), pl = pool();
  // the figure box comes in with the first track; until then "select tracks" stands where it will be
  if (S.locus && p.length && inF.length) reveal('right');
  byId('splitV').hidden = byId('right').hidden; byId('pickT').hidden = byId('tracks').hidden || !byId('right').hidden; if (!byId('pickT').hidden && !byId('pickT').dataset.in) { byId('pickT').dataset.in = 1; byId('pickT').classList.add('pop'); }
  byId('inList').innerHTML = inF.map(k => `<div class="trow in" role="listitem" tabindex="0" data-k="${esc(k)}" title="${esc(fullName(tk(k)))}" aria-label="${esc(fullName(tk(k)))}"><span class="grip"></span>${label(tk(k))}${tk(k).T.viz === 'genes' ? `<button class="ib" data-dup="${esc(k)}" aria-label="Another gene annotation below" title="Another one below">&#10697;</button>` : ''}<button class="ib" data-rm="${esc(k)}" aria-label="Remove ${esc(fullName(tk(k)))}">&#10005;</button></div>`).join('');
  // tracks kept for later: in the figure's list, but not drawable with the current choice
  const hid = S.fig.filter(k => { const t = tk(k); return !t.ok && t.T.kind !== 'comparison' && !(t.T.kind === 'sample' && !S.picked.includes(t.sub)); }), why = t => t.T.gene && !(S.locus && S.locus.gene) ? 'needs a gene'
    : t.T.kind === 'comparison' ? (S.picked.length !== 2 ? 'needs exactly two samples' : 'one of the two samples has no such data')
    : t.T.kind === 'sample' && S.picked.includes(t.sub) ? 'this sample has no such data' : 'its sample is not picked';
  byId('hidN').hidden = !hid.length; byId('hidN').textContent = hid.length + ' hidden';
  byId('hidN').title = hid.map(k => { const t = tk(k); return t.name + (t.T.kind === 'sample' ? ' \u00B7 ' + t.sub : '') + ': ' + why(t); }).join('\n') + '\nThey come back when they can be drawn.';
  byId('clearAll').disabled = !inF.length;

  addable = renderAdd(pl);
  byId('addAll').disabled = !addable.length;
  // figure or gene info
  const gene = S.locus && S.locus.gene;
  if (!gene && S.tab === 'info') S.tab = 'figure';
  byId('viewFlip').hidden = !gene; setFlip(byId('viewFlip'), S.tab === 'info');
  const st = byId('stage'), fig = S.tab === 'figure'; st.classList.toggle('info', !fig);
  byId('where').hidden = !fig || !S.locus; if (fig && S.locus) renderWhere();
  if (!fig && (!S.info || S.info.gene !== gene) && S.infoAsk !== gene) { S.infoAsk = gene; api.geneinfo(gene).then(j => { j.gene = gene; S.info = j; S.infoAsk = null; if (S.tab === 'info') render(); }); }
  // the canvas is only rewritten when what it shows changes (no flicker, no replayed animation)
  const key = [S.tab, S.busy, S.drawn && S.drawn.id, S.fail, inF.length > 0, !fig && S.info && S.info.gene === gene ? gene + S.picked.join() : '', S.hist.length > 0].join('|');
  if (st.dataset.key !== key) { st.dataset.key = key;
    if (!fig) st.innerHTML = S.info && S.info.gene === gene ? infoHTML(S.info) : '<div class="centre"><div class="spin"></div></div>';
    else if (S.busy) st.innerHTML = '<div class="centre"><div class="spin"></div>' + (api.cancel ? '<button class="sbtn" id="cancelBtn">Cancel</button>' : '') + '</div>';
    else if (S.drawn) st.innerHTML = S.drawn.html;
    else if (!inF.length) st.innerHTML = '<div class="centre">' + THINT + '</div>';
    else st.innerHTML = `<div class="centre">${S.fail ? `<div class="err">${esc(S.fail)}</div>` : ''}<button class="gen pop" id="gen">Generate figure</button></div>`; }
  // Download: while a PNG or JPEG is being drawn the button is its Cancel (where the server can stop a draw)
  { const b = byId('dlBtn'); b.disabled = dlBusy ? !api.cancel : !S.drawn; b.classList.toggle('busy', !!dlBusy);
    byId('dlLbl').textContent = !dlBusy ? 'Download' : api.cancel ? 'Cancel ' + dlWhat : 'Drawing ' + dlWhat + '\u2026'; if (!S.drawn || dlBusy) dlMenu(false); }
  st.classList.toggle('blank', fig && !S.drawn); applyFit();
  // figures drawn this session
  gsRender();
  reveal('hist');
  restoreFocus(fk);
  if (had === byId('addAll') && byId('addAll').disabled) byId('clearAll').focus(); else if (had === byId('clearAll') && byId('clearAll').disabled) byId('addAll').focus();
  if (drag && drag.y != null && !drag.row.isConnected && dragRow()) dragTo(drag.y);   // a redraw under a drag: the row goes back under the pointer
  fitNames();
  enterHint();
}
function fitNames() { for (const el of document.querySelectorAll('#roles .chunk .trunc')) { el.style.width = '';
  const r = document.createRange(); r.selectNodeContents(el); const rs = [...r.getClientRects()];
  if (new Set(rs.map(x => Math.round(x.top))).size > 1) el.style.width = Math.ceil(Math.max(...rs.map(x => x.width)) + 1) + 'px'; } }

// Zoom and move buttons, the chromosome, and the genes around the window (click one to go there).
const IC = { zout: '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"><circle cx="7" cy="7" r="4.5"/><path d="M10.5 10.5 14 14M5 7h4"/></svg>',
  zin: '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"><circle cx="7" cy="7" r="4.5"/><path d="M10.5 10.5 14 14M5 7h4M7 5v4"/></svg>',
  left: '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M10 3 5 8l5 5"/></svg>',
  right: '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="m6 3 5 5-5 5"/></svg>' };
function stripRange() { const v = S.view, span = Math.max((v.e - v.s) * 7, 200000), c = (v.s + v.e) / 2; return { a: c - span / 2, b: c + span / 2, span }; }
function renderWhere() {
  const L = S.locus, v = S.view, span = v.e - v.s, { a, span: wide } = stripRange(), pct = x => (x - a) / wide * 100;
  const atHome = v.s === L.view.s && v.e === L.view.e;
  const genes = S.nearby.filter(g => g.e > a && g.s < a + wide).map(g => { const l = Math.max(0, pct(g.s)), r = Math.min(100, pct(g.e));
    return `<button class="gn${g.sym === L.gene ? ' cur' : ''}" data-gene="${esc(g.sym)}" style="left:${l}%;width:${Math.max(0.5, r - l)}%" title="${esc(g.sym)}, ${esc(L.chr)}:${fmt(g.s)}-${fmt(g.e)}" aria-label="Go to ${esc(g.sym)}"><span class="gb"></span><span class="gl">${esc(g.sym)}</span></button>`; }).join('');
  byId('where').innerHTML = `<button class="tbtn sq" data-zoom="3" aria-label="Zoom out 3 times" title="Zoom out 3x" ${span < Math.min(SPAN_MAX, L.chrLen) ? '' : 'disabled'}>${IC.zout}</button>
    <span class="span" title="Width of the window">${fmtSpan(span)}</span>
    <button class="tbtn sq" data-zoom="0.3333333" aria-label="Zoom in 3 times" title="Zoom in 3x" ${span > SPAN_MIN ? '' : 'disabled'}>${IC.zin}</button>
    <button class="tbtn sq" data-pan="-1" aria-label="Move left by half a window" title="Move left" ${v.s > 1 ? '' : 'disabled'}>${IC.left}</button>
    <button class="tbtn sq" data-pan="1" aria-label="Move right by half a window" title="Move right" ${v.e < L.chrLen ? '' : 'disabled'}>${IC.right}</button>
    <button class="tbtn" data-fit="1" ${atHome ? 'disabled style="visibility:hidden"' : ''} title="Back to ${esc(L.label)}"><span class="trunc" style="max-width:86px">${L.gene ? esc(L.gene) : 'Fit'}</span></button>
    <div class="chrom" title="${esc(L.chr)}:${fmt(v.s)}-${fmt(v.e)}"><span>${esc(L.chr)}</span><span class="cbar"><i style="left:${(v.s + v.e) / 2 / L.chrLen * 100}%"></i></span></div>
    <div class="nb" id="nb" title="${fmtSpan(wide)} around the window"><div class="axis"></div><div class="win" style="left:${pct(v.s)}%;width:${Math.max(0.6, pct(v.e) - pct(v.s))}%"></div>${genes}</div>`;
  // gene names that would collide are hidden; the current gene and the widest genes are placed first
  const nb = byId('nb'), box = nb.getBoundingClientRect(), placed = [];
  [...nb.querySelectorAll('.gn')].sort((x, y) => (y.classList.contains('cur') - x.classList.contains('cur')) || (y.offsetWidth - x.offsetWidth)).forEach(g => {
    const l = g.querySelector('.gl'); let r = l.getBoundingClientRect();
    if (r.left < box.left) { l.style.left = '0'; l.style.transform = 'none'; r = l.getBoundingClientRect(); }
    if (r.right > box.right) { l.style.left = 'auto'; l.style.right = '0'; l.style.transform = 'none'; r = l.getBoundingClientRect(); }
    if (placed.some(o => r.left < o.right + 6 && r.right > o.left - 6)) g.classList.add('hide-l'); else placed.push(r);
  });
}
let nearbySeq = 0;
async function loadNearby() { if (!S.locus) return; const { a, span } = stripRange(), seq = ++nearbySeq, chr = S.locus.chr;
  try { const g = await api.nearby(chr, a - span, a + 2 * span); if (Array.isArray(g) && seq === nearbySeq && S.locus && S.locus.chr === chr) { S.nearby = g; if (S.tab === 'figure') renderWhere(); } } catch (e) { /* the strip just stays as it is */ } }
function setView(s, e) { const L = S.locus, w = Math.min(Math.max(e - s, SPAN_MIN), Math.min(SPAN_MAX, L.chrLen)); s = Math.round(Math.min(Math.max(1, s), L.chrLen - w)); changed(() => { S.view = { s, e: Math.round(s + w) }; }); mirror(); loadNearby(); }
// The Region boxes always say the window on screen.
function mirror() { if (S && S.locus && S.view) { byId('chr').value = S.locus.chr; byId('range').value = `${fmt(S.view.s)}-${fmt(S.view.e)}`; } }

// Gene info: transcripts, the canonical protein, a link to UniProt.
function infoHTML(j) {
  if (typeof j.error === 'string') return `<div class="gi"><h3>${esc(j.gene)}</h3><div class="err">${esc(j.error)}</div></div>`;
  const num = v => v == null || v === '' ? '' : (typeof v === 'number' ? fmt(v) : esc(v)), p = j.protein || {}, tx = j.transcripts || [], ft = p.features || [];
  const tpm = S.picked.filter(n => tx.some(t => t.tpm && t.tpm[n] != null));   // TPM columns for the picked samples, when the answer carries them
  return `<div class="gi"><h3>${esc(j.gene)} ${j.mock ? '<span class="mock">made-up values</span>' : ''}</h3>
    ${j.chr ? `<div>${esc(j.chr)}:${fmt(j.s)}-${fmt(j.e)}${j.strand ? ', ' + esc(j.strand) + ' strand' : ''}</div>` : ''}
    <h4>Transcripts</h4><table><thead><tr><th>Transcript</th><th>Biotype</th><th class="num">Exons</th><th class="num">Length (bp)</th><th class="num">Protein (aa)</th>${tpm.map(n => `<th class="num" title="${esc(n)}">TPM, ${esc(n)}</th>`).join('')}<th></th></tr></thead><tbody>${tx.map(t =>
      `<tr><td><b>${esc(t.name || t.id || '')}</b>${t.name && t.id ? `<small>${esc(t.id)}</small>` : ''}</td><td>${esc(String(t.biotype || '').replace(/_/g, ' '))}</td><td class="num">${num(t.exons)}</td><td class="num">${num(t.length)}</td><td class="num">${num(t.aa)}</td>${tpm.map(n => `<td class="num">${t.tpm[n] == null ? '' : t.tpm[n].toFixed(1)}</td>`).join('')}<td>${(t.tags || []).map(x => `<span class="tag">${esc(x)}</span>`).join('')}</td></tr>`).join('')}</tbody></table>
    <h4>Canonical protein${j.uniprot_url ? ` &nbsp;<a href="${esc(j.uniprot_url)}" target="_blank" rel="noopener">Open in UniProt &#8599;</a>` : ''}</h4>
    ${p.length || p.accession ? `<table><thead><tr><th>Accession</th><th>Name</th><th class="num">Length (aa)</th></tr></thead><tbody><tr><td>${p.accession && j.uniprot_url ? `<a href="${esc(j.uniprot_url)}" target="_blank" rel="noopener">${esc(p.accession)}</a>` : esc(p.accession || '')}</td><td>${esc(p.name || '')}</td><td class="num">${num(p.length)}</td></tr></tbody></table>` : j.protein_status ? `<div>${esc(j.protein_status)}</div>` : ''}
    ${ft.length ? `<h4>Protein features</h4><table><thead><tr><th>Feature</th><th class="num">Start</th><th class="num">End</th><th>Source</th></tr></thead><tbody>${ft.map(f => `<tr><td>${esc(f.name)}</td><td class="num">${num(f.start)}</td><td class="num">${num(f.end)}</td><td>${esc(f.source || '')}</td></tr>`).join('')}</tbody></table>` : ''}</div>`;
}

/* ------------------------------------------------------------------ generate, cache, downloads */

// The picture as drawn in the page: scaled to the width, whatever size the SVG declares.
function fitSvg(svg) { const d = new DOMParser().parseFromString(svg, 'image/svg+xml').documentElement; if (d.nodeName !== 'svg') return '<div class="centre"><div class="err">The server did not return a picture.</div></div>';
  if (!d.getAttribute('viewBox')) d.setAttribute('viewBox', `0 0 ${parseFloat(d.getAttribute('width')) || 1000} ${parseFloat(d.getAttribute('height')) || 600}`);
  d.removeAttribute('width'); d.removeAttribute('height'); d.setAttribute('class', 'fig pop'); return new XMLSerializer().serializeToString(d); }
async function generate() {
  if (!S.locus || !shown().length || (S.drawn && S.drawn.id === sig()) || S.busy) return;
  { const hit = S.hist.find(h => h.id === sig()); if (hit) { S.drawn = hit; return render(); } }   // drawn before this session: it comes straight back
  if (dlBusy) return toast('The ' + dlWhat + ' is still being drawn.');
  // Everything the picture, its label and its file name need is fixed here, before the wait: a change made while it draws must not leak in.
  const S0 = S, seq0 = setSeq, L = S.locus, v = { ...S.view }, picked = [...S.picked], names = picked.map(sn).join(', '), who = picked.length > 2 ? (names.length <= 30 ? names : picked.length + ' samples') : picked.length === 2 ? pairShort() : picked[0] || 'annotation only';   // a card says what tells it apart, in few words
  const state = JSON.parse(JSON.stringify({ region: S.region, locus: L, view: v, picked, fig: S.fig, set: S.set, hand: S.hand || [], rm: S.rm || [], exp: S.exp }));
  const snap = { locus: L, view: v, picked, set: { ...S.set }, items: shown().map(k => tk(k)), stem: [L.gene || L.chr, v.s, v.e, ...(picked.length > 3 ? [picked.length + '-samples'] : picked)].join('_').replace(/[^\w.+-]+/g, '-') };
  const req = { id: P.id, gene: L.gene, chr: L.chr, s: v.s, e: v.e, picked, tracks: snap.items.map(t => t.T.kind === 'sample' ? { type: t.T.id, sample: t.sub } : { type: t.T.id }) };
  S.tab = 'figure'; S.busy = true; S.fail = null; render();
  let id = sig(), out; try { out = await api.render(req, snap); } catch (e) { out = { error: String(e.message || e) }; }
  S0.busy = false;
  if (out.cancelled) return render();
  if (out.error || !out.svg) { if (S === S0 && sig() === id) S.fail = typeof out.error === 'string' ? out.error : 'Nothing came back.'; return render(); }
  // The server says which settings it drew with. A setting changed just before Generate may not have been announced to the page
  // yet: the figure is always filed under the settings it has. Those become the page's current settings only when the server has
  // said nothing newer since the request went out (a "pt_settings" that arrived meanwhile is later news than this answer).
  let took = false;
  if (out.key && out.key !== state.set.key) { state.set = { key: out.key }; const a = JSON.parse(id); a[a.length - 1] = state.set; id = JSON.stringify(a);
    if (setSeq === seq0 && SETKEY !== out.key) { SETKEY = out.key; if (S === S0) S.set = { key: out.key }; took = true; } }
  S0.count = (S0.count || 0) + 1;
  const moved = L.gene && L.view && (v.s !== L.view.s || v.e !== L.view.e);
  const entry = { id, no: S0.count, at: moved ? (v.s / 1e6).toFixed(2) + ' Mb' : '', label: L.gene || `${L.chr}:${(v.s / 1e6).toFixed(2)} Mb`, sub: `${who} \u00B7 ${fmtSpan(v.e - v.s)}`,
    title: `${L.gene ? L.gene + ', ' : ''}${L.chr}:${fmt(v.s)}-${fmt(v.e)}, ${picked.join(picked.length === 2 ? ' vs ' : ', ')}, ${snap.items.length} tracks`,
    svg: out.svg, html: fitSvg(out.svg), script: out.script || '', job: out.job, stem: snap.stem, state, thumb: '' }; gsThumb(entry);
  S0.hist = [entry, ...S0.hist.filter(h => h.id !== id)];
  // Over 40, the oldest go, but never the one on screen.
  for (let i = S0.hist.length - 1; S0.hist.length > 40 && i > 0; i--) if (S0.hist[i] !== S0.drawn) S0.hist.splice(i, 1);
  if (L.gene && !S0.recent.includes(L.gene)) S0.recent.unshift(L.gene);
  if (S === S0 && sig() === id) S.drawn = entry;   // the choice may have moved on while it was drawing: then it only joins the strip
  else if (S === S0 && took) { const now = sig(); S.drawn = S.hist.find(h => h.id === now) || null; }   // the settings moved under another choice: what is on screen is looked up again
  render();
}
function reopen(i) { if (S.busy) return; const h = S.hist[i]; Object.assign(S, JSON.parse(JSON.stringify(h.state))); S.set = { key: SETKEY };   // the settings are the panel's, not the old figure's
  setRegion(S.region); byId('gene').value = S.locus.gene || ''; if (S.locus.gene) S.lastGene = S.locus.gene; mirror();
  S.tab = 'figure'; S.err = null; S.gap = null; S.was = null; S.drawn = h; S.fail = null; render(); loadNearby(); }
/* ------------------------------------------------------------------ downloads */
// Where the browser has a save dialog (Chrome, Edge) the user picks the folder and the name. Elsewhere (Firefox, Safari), or where
// the dialog is refused (the page inside another site's frame), the browser downloads the file as it always did. The dialog is asked
// for at the click, before anything slow: a browser opens it only straight after a click or a key. Closing it saves nothing and says
// nothing. make() gives the content (text or a Blob) and may take a while.
// A figure's SVG as a file to take away: without svglite's text width lock (textLength), which the page keeps so text has
// R's widths in any font, and which would squeeze an edited label to its old width (R/viz-svg.R).
const svgFile = svg => svg.replace(/<text\b[^>]*>/g, t => t.replace(/ textLength='[^']*' lengthAdjust='[^']*'/g, ''));
const KIND = { svg: ['SVG image', 'image/svg+xml', ['.svg']], png: ['PNG image', 'image/png', ['.png']], jpg: ['JPEG image', 'image/jpeg', ['.jpg', '.jpeg']], R: ['R script', 'text/plain', ['.R']] };
async function saveFile(name, kind, make) {
  const [what, type, exts] = KIND[kind]; let file = null;
  if (window.showSaveFilePicker) { try { file = await window.showSaveFilePicker({ suggestedName: name, types: [{ description: what, accept: { [type]: exts } }] }); }
    catch (e) { if (e && e.name === 'AbortError') return; file = null; } }
  let data = null, why = '';
  try { data = await make(); } catch (e) { why = String(e && e.message || e); }
  if (data == null) { if (file && file.remove) { try { if ((await file.getFile()).size === 0) await file.remove(); } catch (e) {} }   // an empty file the dialog made is not left behind; a file chosen to be replaced is left as it was
    return toast(why || 'Nothing to save.'); }
  if (!file) return download(name, data, type);
  try { const w = await file.createWritable(); await w.write(data instanceof Blob ? data : new Blob([data], { type })); await w.close(); toast('Saved ' + midCut(file.name, 54)); }
  catch (e) { toast('Could not save ' + midCut(file.name, 40) + '.'); }
}
// A PNG or JPEG is drawn again by R, from the request the server remembers for this figure, so it has the settings the figure was
// drawn with. One at a time, and not while a figure is being generated (R draws one thing at a time): the Download button is
// the draw's Cancel meanwhile, and Generate waits.
let dlBusy = 0, dlWhat = '';
async function raster(h, format, label) {
  if (!h.job || !api.raster) throw new Error('Generate this figure again to save it as ' + label + '.');
  dlBusy++; dlWhat = label; render(); toast('Drawing the ' + label + '\u2026');
  let out; try { out = await api.raster(h.job, format); } catch (e) { out = { error: String(e.message || e) }; } finally { dlBusy--; render(); }
  if (out.error || !out.data) throw new Error(typeof out.error === 'string' ? out.error : 'Nothing came back.');
  const bin = atob(out.data), bytes = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  return new Blob([bytes], { type: KIND[format === 'png' ? 'png' : 'jpg'][1] }); }
function dlMenu(on) { const m = byId('dlMenu'); if (m.hidden === !on) return; m.hidden = !on; byId('dlBtn').setAttribute('aria-expanded', String(on)); if (on) byId('dlSvg').focus(); }
function download(name, text, type) { const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([text], { type })); a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000); toast('Downloaded ' + midCut(name, 54)); }

/* ------------------------------------------------------------------ the search box */
let suggSeq = 0;
async function updateSugg() {
  const box = byId('gene'), text = box.value.trim(), seq = ++suggSeq; let items = [];
  try { items = await api.search(text); } catch (e) { items = []; }
  if (!Array.isArray(items)) items = [];   // an error answer (the connection is lost): no suggestions
  if (seq !== suggSeq) return;
  if (!text) { const rec = S.recent.filter(r => !/:/.test(r)); items = [...rec.map(r => ({ sym: r, recent: true })), ...items.filter(g => g && !rec.includes(g.sym))]; }
  S.sugg = items.filter(Boolean); S.suggFor = text; S.suggAt = Math.min(S.suggAt, Math.max(0, S.sugg.length - 1));
  const ul = byId('sugg'), open = document.activeElement === box && S.sugg.length && !(S.locus && S.locus.label.toLowerCase() === text.toLowerCase());
  ul.hidden = !open; box.setAttribute('aria-expanded', String(!!open)); if (!open) return;
  let head = '';
  ul.innerHTML = S.sugg.map((g, i) => { const h = text ? '' : g.recent ? 'Drawn this session' : 'Suggested', hd = h && h !== head ? `<li class="shead" aria-hidden="true">${h}</li>` : ''; head = h || head;
    return hd + `<li role="option" id="sg${i}" data-i="${i}" class="${i === S.suggAt ? 'active' : ''}" aria-selected="${i === S.suggAt}"><span class="g">${esc(g.sym)}</span><span class="c">${g.chr ? `${esc(g.chr)}:${fmt(g.s)}-${fmt(g.e)} · ${fmtSpan(g.e - g.s)}` : ''}</span></li>`; }).join('');
  box.setAttribute('aria-activedescendant', 'sg' + S.suggAt);
}
const closeSugg = () => { byId('sugg').hidden = true; byId('gene').setAttribute('aria-expanded', 'false'); };
// Use what is in the box. Returns true when there is a locus afterwards.
// What the Region boxes say, as "chr:start-end". A pasted "chr2:1-2" also sets the chromosome. false: it is not a range.
function regionText() { let t = byId('range').value.replace(/\s/g, ''); if (!t) return '';
  const m = t.match(/^([^:]+):(.+)$/), c = m && [...byId('chr').options].find(o => o.value.toLowerCase() === m[1].toLowerCase());
  if (c) t = m[2];
  else if (m && /^[\d,]+-[\d,]+$/.test(m[2])) { S.err = `No chromosome '${m[1]}' in this project.`; toast(S.err); render(); return false; }
  if (!/^[\d,]+-[\d,]+$/.test(t)) { S.err = 'A range looks like 165,430,000-165,620,000.'; toast(S.err); render(); return false; }
  if (c) { byId('chr').value = c.value; byId('range').value = t; }
  return byId('chr').value + ':' + t; }
let commitSeq = 0;
async function commit(text) {
  const inBox = document.activeElement === byId('gene') || document.activeElement === byId('range');
  if (text == null) text = S.region ? regionText() : byId('gene').value.trim();
  if (text === false) return false;
  if (!text) return !!S.locus && !S.err;                      // an empty box: what is on screen stands
  const norm = t => t.replace(/[\s,]/g, '').toLowerCase();
  if (S.locus && !S.err && (norm(text) === norm(S.locus.label) || norm(text) === norm(`${S.locus.chr}:${S.view.s}-${S.view.e}`))) return true;   // nothing new
  const seq = ++commitSeq; let L; try { L = await api.locus(text); } catch (e) { L = { error: 'The server did not answer.' }; }
  if (seq !== commitSeq) return false;                        // a later entry has overtaken this one
  if (!L.error && L.kind === 'region' && L.e - L.s > SPAN_MAX) L = { error: `Up to ${fmtSpan(SPAN_MAX)} at a time.` };
  if (!L.error && L.kind === 'region' && L.e - L.s < SPAN_MIN) L = { error: `At least ${fmtSpan(SPAN_MIN)}.` };
  if (L.error) { if (S.region && S.locus) byId('chr').value = S.locus.chr; S.err = L.error; toast(L.error); render(); return false; }
  setRegion(L.kind === 'region'); if (L.gene) { byId('gene').value = L.label; S.lastGene = L.gene; }
  changed(() => { S.err = null; S.locus = L; S.view = { s: L.view.s, e: L.view.e }; if (L.gene !== (S.info && S.info.gene)) S.info = null; });
  mirror(); loadNearby(); if (inBox && (!document.activeElement || document.activeElement === document.body)) (S.region ? byId('range') : byId('gene')).focus(); return true;
}
// Enter in a box: use what is typed; when that changed nothing, draw.
async function enterIn() { const before = sig(); if (await commit() && sig() === before) generate(); }

/* ------------------------------------------------------------------ wiring */
function setFlip(el, on) { el.setAttribute('aria-checked', String(on)); }
function onFlip(el, fn) { const set = on => { if ((el.getAttribute('aria-checked') === 'true') === on) return; setFlip(el, on); fn(on); }, flip = () => set(el.getAttribute('aria-checked') !== 'true');
  el.addEventListener('click', flip);
  el.addEventListener('keydown', e => { if (e.ctrlKey || e.metaKey || e.altKey) return;
    if (e.key === ' ' || e.key === 'Enter') { e.preventDefault(); flip(); } else if (e.key === 'ArrowLeft') { e.preventDefault(); set(false); } else if (e.key === 'ArrowRight') { e.preventDefault(); set(true); } }); }
// Fit width (the figure as wide as its box, scrolling when it is tall) or fit height (the whole figure inside the box, at the
// top, centred). Only how the figure is shown: downloads and saved figures are the same either way. Kept in the browser, as light and dark are.
let FITH = true; try { FITH = localStorage.getItem('pt-fit') !== 'width'; } catch (e) {}
function applyFit() { const st = byId('stage'), on = FITH && !!S && S.tab === 'figure', f = st.querySelector('svg.fig'); st.classList.toggle('fit-h', on);
  if (f) { if (on) f.setAttribute('preserveAspectRatio', 'xMidYMin meet'); else f.removeAttribute('preserveAspectRatio'); } }
// The switch's first half is the whole figure (fit height), so it is on for fit width.
setFlip(byId('fitFlip'), !FITH);
onFlip(byId('fitFlip'), on => { FITH = !on; try { localStorage.setItem('pt-fit', on ? 'width' : 'height'); } catch (e) {} applyFit(); });
function setRegion(on) { S.region = on; setFlip(byId('locusFlip'), on); byId('gene').hidden = on; byId('chr').hidden = !on; byId('range').hidden = !on; closeSugg(); }
// Light or dark is chosen in the Project settings dialog (R/app-project-actions.R). R draws that dialog anew each time it is
// opened, so its switch is found when it is used, and set to the mode on screen when the dialog shows.
const isDark = () => document.documentElement.getAttribute('data-pt-mode') === 'dark';
function setMode(dark) { document.documentElement.setAttribute('data-pt-mode', dark ? 'dark' : 'light'); try { localStorage.setItem('pt-mode', dark ? 'dark' : 'light'); } catch (e) {}
  const f = byId('modeFlip'); if (f) setFlip(f, dark); }
document.addEventListener('click', e => { if (e.target.closest('#modeFlip')) setMode(!isDark()); });
document.addEventListener('keydown', e => { if (!e.target.closest || !e.target.closest('#modeFlip') || e.ctrlKey || e.metaKey || e.altKey) return;
  if (e.key === ' ' || e.key === 'Enter') { e.preventDefault(); setMode(!isDark()); } else if (e.key === 'ArrowLeft') { e.preventDefault(); setMode(false); } else if (e.key === 'ArrowRight') { e.preventDefault(); setMode(true); } });
if (window.jQuery) window.jQuery(document).on('show.bs.modal', () => { const f = byId('modeFlip'); if (f) setFlip(f, isDark()); });
// Gene or Region only changes which boxes show. Region starts from the window on screen and nothing is redrawn until a range is entered;
// back on Gene the box says the gene on screen (after a typed region: the gene before it, which is gone back to).
onFlip(byId('locusFlip'), async on => { setRegion(on); S.err = null;
  if (on) mirror();
  else if (S.locus && S.locus.gene) byId('gene').value = S.locus.gene;
  else if (S.lastGene) { byId('gene').value = S.lastGene; await commit(S.lastGene); }
  else byId('gene').value = '';
  render(); (S.region ? byId('range') : byId('gene')).focus(); });
onFlip(byId('viewFlip'), on => { S.tab = on ? 'info' : 'figure'; render(); });

byId('gene').addEventListener('input', () => { if (!S) return; S.suggAt = 0; S.err = null; byId('clear').hidden = !byId('gene').value; byId('gene').setAttribute('aria-invalid', 'false'); clearTimeout(updateSugg.t); updateSugg.t = setTimeout(updateSugg, api.wait || 0); });
let blurT = 0;
byId('gene').addEventListener('focus', () => { if (!S) return; clearTimeout(blurT); byId('gene').select(); S.suggAt = 0; updateSugg(); });
byId('gene').addEventListener('click', () => { if (byId('sugg').hidden) updateSugg(); });
byId('gene').addEventListener('blur', () => { blurT = setTimeout(closeSugg, 150); });
byId('gene').addEventListener('change', () => { if (byId('gene').value.trim()) commit(); });   // typed and left: it is used, or marked wrong, never left unsaid
byId('gene').addEventListener('keydown', async e => {
  if (!S) return;   // typed before the project arrived
  const open = !byId('sugg').hidden;
  if (e.key === 'ArrowDown' && S.sugg.length) { e.preventDefault(); S.suggAt = open ? (S.suggAt + 1) % S.sugg.length : 0; updateSugg(); }
  else if (e.key === 'ArrowUp' && open) { e.preventDefault(); S.suggAt = (S.suggAt - 1 + S.sugg.length) % S.sugg.length; updateSugg(); }
  else if (e.key === 'Escape') closeSugg();
  else if (e.key === 'Enter' && !e.ctrlKey && !e.metaKey) { e.preventDefault(); const typed = byId('gene').value.trim().toLowerCase();
    if (open && S.sugg.length && S.suggFor === byId('gene').value.trim()) byId('gene').value = (S.sugg.find(g => g.sym.toLowerCase() === typed) || S.sugg[S.suggAt] || { sym: byId('gene').value }).sym;   // a name typed in full wins over the highlighted row
    closeSugg(); enterIn(); }
});
byId('sugg').addEventListener('mousedown', e => { const li = e.target.closest('li[data-i]'); if (!li) return; e.preventDefault(); byId('gene').value = S.sugg[+li.dataset.i].sym; closeSugg(); commit(); });
byId('range').addEventListener('input', () => { byId('clear').hidden = !byId('range').value; });
byId('range').addEventListener('keydown', e => { if (e.key === 'Enter' && !e.ctrlKey && !e.metaKey) { e.preventDefault(); enterIn(); } });
byId('range').addEventListener('change', () => commit());
byId('chr').addEventListener('change', () => commit());
byId('clear').addEventListener('click', () => { const box = S.region ? byId('range') : byId('gene'); box.value = ''; S.err = null; render(); box.focus(); });
// A newly picked sample joins every kind of per-sample track the figure already shows, next to its own kind.
function follow(n) {
  if (n) for (const T of P.types) { if (T.kind !== 'sample' || !T.samples.includes(n) || inFig(T.id + '|' + n) || (S.rm || []).includes(T.id + '|' + n)) continue;
    if (S.fig.some(x => tk(x).T === T)) place(T.id + '|' + n); }
}
function flipSample(n) { const top = byId('chips').scrollTop;
  changed(() => { const i = S.picked.indexOf(n);
    reranked(() => { if (i >= 0) { S.picked.splice(i, 1); S.gap = S.picked.length === 1 ? i : null;   // remember which side of the pair it was on
        if (S.picked.length === 2 && S.was && S.picked.every(x => S.was.includes(x))) S.picked.splice(0, 2, ...S.was); }   // the same pair again: the same roles
      else { if (S.picked.length === 2) S.was = [...S.picked];
        if (S.picked.length === 1 && S.gap === 0) S.picked.unshift(n); else S.picked.push(n); S.gap = null; if (S.picked.length > 2) S.picked.sort((a, b) => rank(a) - rank(b)); } });
    follow(i < 0 ? n : null); });
  byId('chips').scrollTop = top; }
const setSample = (n, on) => { if (S.picked.includes(n) !== on) flipSample(n); };
const setAll = on => { if (on) P.samples.forEach(x => setSample(x.name, true)); else changed(() => { reranked(() => { S.picked.splice(0); S.gap = null; }); follow(null); }); };
// A click switches one sample. Shift+click switches every sample between the last one clicked and this one.
// Press on a sample and sweep over others: they all take the state the first one got ("painting").
let paint = null, lastChip = null, paintDone = false;
byId('chips').addEventListener('click', e => { const b = e.target.closest('.chip'); if (!b) return; if (paintDone) { paintDone = false; return; } const n = b.dataset.s;
  if (e.shiftKey && lastChip && lastChip !== n) { const ns = P.samples.map(x => x.name), a = ns.indexOf(lastChip), z = ns.indexOf(n), on = !S.picked.includes(n); for (const x of ns.slice(Math.min(a, z), Math.max(a, z) + 1)) setSample(x, on); }
  else flipSample(n);
  lastChip = n; });
byId('chips').addEventListener('pointerdown', e => { const b = e.target.closest('.chip'); if (b && !e.button && !e.shiftKey) paint = { first: b.dataset.s, on: !S.picked.includes(b.dataset.s), done: new Set(), moved: false }; });
window.addEventListener('pointermove', e => { if (!paint) return; const el = document.elementFromPoint(e.clientX, e.clientY), b = el && el.closest && el.closest('#chips .chip'); if (!b) return; const n = b.dataset.s;
  if (!paint.moved) { if (n === paint.first) return; paint.moved = true; paint.done.add(paint.first); setSample(paint.first, paint.on); }
  if (!paint.done.has(n)) { paint.done.add(n); setSample(n, paint.on); lastChip = n; } });
window.addEventListener('pointerup', () => { if (paint && paint.moved) { paintDone = true; setTimeout(() => { paintDone = false; }, 0); } paint = null; });
const hotSample = i => flipSample(P.samples[i].name);   // a number switches its sample on or off, nothing more
byId('allBtn').addEventListener('click', () => { setAll(S.picked.length !== P.samples.length); byId('allBtn').focus(); });
const hidSay = () => toast(byId('hidN').title.split('\n').slice(0, 6).join(' \u00B7 '));
byId('hidN').addEventListener('click', hidSay); byId('hidN').addEventListener('keydown', e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); hidSay(); } });
byId('roles').addEventListener('click', e => { if (e.target.closest('#unpick')) { setAll(false); const c = document.querySelector('#chips .chip'); c && c.focus(); return; }
  if (e.target.closest('#swap')) { changed(() => reranked(() => S.picked.reverse())); byId('swap').focus(); } });
byId('addAll').addEventListener('click', () => changed(() => { for (const k of addable) if (!inFig(k)) place(k); }));
byId('clearAll').addEventListener('click', () => changed(() => { S.fig = []; S.hand = []; }));
byId('addList').addEventListener('click', e => { const x = e.target.closest('[data-x]'); if (x) { S.exp[x.dataset.x] = !S.exp[x.dataset.x]; return render(); }
  const b = e.target.closest('[data-g]'); if (b) toggle(REG[+b.dataset.g]); });
byId('inList').addEventListener('click', e => { const b = e.target.closest('[data-rm]'); if (b) return changed(() => { S.fig = S.fig.filter(k => k !== b.dataset.rm); });
  const d = e.target.closest('[data-dup]'); if (d) changed(() => { const base = d.dataset.dup.replace(/~\d+$/, ''); let n = 2; while (S.fig.includes(base + '~' + n)) n++; S.fig.splice(S.fig.indexOf(d.dataset.dup) + 1, 0, base + '~' + n); }); });
byId('where').addEventListener('click', e => { const b = e.target.closest('button'); if (!b || b.disabled) return; const v = S.view, w = v.e - v.s, c = (v.s + v.e) / 2, d = b.dataset;
  if (d.zoom) { const nw = w * +d.zoom; setView(c - nw / 2, c + nw / 2); } else if (d.pan) setView(v.s + +d.pan * w / 2, v.e + +d.pan * w / 2);
  else if (d.fit) setView(S.locus.view.s, S.locus.view.e); else if (d.gene) { if (S.region) setRegion(false); byId('gene').value = d.gene; commit(d.gene); } });

// Dragging a track: it changes place as soon as the pointer is inside another row.
const inRows = () => [...byId('inList').querySelectorAll('.trow.in')];
byId('inList').addEventListener('pointerdown', e => { const row = e.target.closest('.trow.in'); if (!row || e.target.closest('button') || e.button) return;
  drag = { row, key: row.dataset.k, before: inRows().map(r => r.dataset.k).join('\n') }; row.classList.add('dragging'); row.focus(); e.preventDefault(); });
// The list may have been redrawn under a drag (a figure landed): find the dragged track's new row, or give the drag up.
function dragRow() { if (drag && !drag.row.isConnected) { const r = inRows().find(x => x.dataset.k === drag.key); if (!r) { drag = null; return false; }
    drag.row = r; drag.before = inRows().map(x => x.dataset.k).join('\n'); r.classList.add('dragging'); } return !!drag; }
window.addEventListener('keydown', e => { if (e.key === 'Escape' && drag) { const k = drag.key; drag = null; render(); const el = inRows().find(r => r.dataset.k === k); el && el.focus(); } });   // Escape puts the track back
function dragTo(y) { const rows = inRows(), over = rows.find(r => { const b = r.getBoundingClientRect(); return y >= b.top && y < b.bottom; })
    || (y < rows[0].getBoundingClientRect().top ? rows[0] : y >= rows.at(-1).getBoundingClientRect().bottom ? rows.at(-1) : null);
  if (!over || over === drag.row) return;   // in the gap between two rows nothing moves
  if (rows.indexOf(over) > rows.indexOf(drag.row)) over.after(drag.row); else over.before(drag.row);
  const list = byId('inList'), b = drag.row.getBoundingClientRect(), lb = list.getBoundingClientRect(); if (b.top < lb.top) list.scrollTop -= lb.top - b.top + 4; else if (b.bottom > lb.bottom) list.scrollTop += b.bottom - lb.bottom + 4; }
window.addEventListener('pointermove', e => { if (!drag || !dragRow()) return; drag.y = e.clientY; dragTo(e.clientY); });
window.addEventListener('pointerup', () => { if (!drag || !dragRow()) return; const d = drag; drag = null; d.row.classList.remove('dragging');
  const order = inRows().map(r => r.dataset.k).filter(k => S.fig.includes(k)); if (order.join('\n') === d.before) return;
  const k = d.row.dataset.k; changed(() => { let i = 0; S.fig = S.fig.map(x => order.includes(x) ? order[i++] : x); byHand(k); }); const el = inRows().find(r => r.dataset.k === k); el && el.focus(); });

byId('stage').addEventListener('click', e => { if (e.target.id === 'gen') generate(); });
byId('dlBtn').addEventListener('click', () => { if (dlBusy) { if (api.cancel) api.cancel(); return; } dlMenu(byId('dlMenu').hidden); });
byId('dlMenu').addEventListener('click', e => { const b = e.target.closest('button'), h = S.drawn; if (!b) return;
  if (!byId('dlMenu').hidden) { dlMenu(false); byId('dlBtn').focus(); }
  if (!h) return;   // the figure on screen at the click is the one saved, whatever is chosen meanwhile
  if (b.id === 'dlSvg') return saveFile(h.stem + '.svg', 'svg', () => svgFile(h.svg));
  if (b.id === 'dlRs') return saveFile(h.stem + '.R', 'R', () => h.script);
  const png = b.id === 'dlPng', label = png ? 'PNG' : 'JPEG';
  // What is known already is said now, not after a place for the file was chosen.
  if (dlBusy) return toast('The ' + dlWhat + ' is still being drawn.');
  if (S.busy) return toast('Wait for the figure being drawn.');
  if (!h.job || !api.raster) return toast('Generate this figure again to save it as ' + label + '.');
  saveFile(h.stem + (png ? '.png' : '.jpg'), png ? 'png' : 'jpg', () => raster(h, png ? 'png' : 'jpeg', label)); });
byId('dlMenu').addEventListener('keydown', e => { const items = [...byId('dlMenu').querySelectorAll('button')], i = items.indexOf(document.activeElement);
  if (e.ctrlKey || e.metaKey) return;   // Ctrl+S and Ctrl+Enter work from here too
  e.stopPropagation();   // the page's single-key shortcuts are not for a menu
  if (e.key === 'ArrowDown' || e.key === 'ArrowUp') { e.preventDefault(); items[(i + (e.key === 'ArrowDown' ? 1 : items.length - 1)) % items.length].focus(); }
  else if (e.key === 'Home' || e.key === 'End') { e.preventDefault(); items[e.key === 'Home' ? 0 : items.length - 1].focus(); }
  else if (e.key === 'Escape') { e.preventDefault(); dlMenu(false); byId('dlBtn').focus(); }
  else if (e.key === 'Tab') dlMenu(false); });
document.addEventListener('click', e => { if (!e.target.closest('#dl')) dlMenu(false); });
document.addEventListener('click', e => { if (e.target.closest('#cancelBtn') && api.cancel) api.cancel(); });
byId('hist').addEventListener('click', e => { const b = e.target.closest('button'); if (!b) return; const d = b.dataset;
  if (d.open) { if (S.hist[+d.open] !== S.drawn) reopen(+d.open); }   // the one on screen is open already: a click on it changes nothing
  else if (d.hx) { forget(+d.hx);
    const th = [...document.querySelectorAll('#hist [data-open]')]; (th.find(x => +x.dataset.open === Math.min(+d.hx, S.hist.length - 1)) || byId('setBtn')).focus(); } });
for (const id of ['gene', 'range']) for (const ev of ['focus', 'blur', 'input']) byId(id).addEventListener(ev, () => setTimeout(enterHint, 0));
// Splitters. In figure / Add: neither half goes under 130 px. Tracks / canvas: the left column stays between 280 px and 45% of the window (at most 640 px).
{ const main = document.querySelector('.main'), tr = byId('tracks'), LS = k => { try { return +localStorage.getItem(k) || 0; } catch (e) { return 0; } }, save = (k, v) => { try { localStorage.setItem(k, v); } catch (e) {} };
  const capW = w => Math.round(Math.max(280, Math.min(w, Math.min(640, innerWidth * 0.45)))), capH = h => Math.round(Math.max(130, Math.min(h, tr.clientHeight - 10 - 130)));
  const setW = w => { main.style.setProperty('--left-w', capW(w) + 'px'); save('pt-g-left', capW(w)); }, setH = h => { tr.style.setProperty('--top-h', capH(h) + 'px'); save('pt-g-top', capH(h)); };
  if (LS('pt-g-left')) main.style.setProperty('--left-w', capW(LS('pt-g-left')) + 'px');
  if (LS('pt-g-top')) tr.style.setProperty('--top-h', LS('pt-g-top') + 'px');
  let sp = null;
  for (const [el, vert] of [[byId('splitV'), true], [byId('splitH'), false]]) {
    el.addEventListener('pointerdown', e => { if (e.button) return; sp = { el, vert }; el.classList.add('on'); try { el.setPointerCapture(e.pointerId); } catch (x) {} e.preventDefault(); });
    el.addEventListener('pointermove', e => { if (!sp || sp.el !== el) return; if (vert) setW(e.clientX - main.getBoundingClientRect().left - 7 - parseFloat(getComputedStyle(main).paddingLeft)); else setH(e.clientY - tr.getBoundingClientRect().top - 5); });
    const end = () => { if (sp && sp.el === el) { el.classList.remove('on'); sp = null; } }; el.addEventListener('pointerup', end); el.addEventListener('pointercancel', end);
    el.addEventListener('keydown', e => { const d = { ArrowLeft: -16, ArrowUp: -16, ArrowRight: 16, ArrowDown: 16 }[e.key]; if (!d) return; e.preventDefault(); e.stopPropagation();
      if (vert) setW(tr.getBoundingClientRect().width + d); else setH(tr.querySelector('.half').getBoundingClientRect().height + d); }); }
  addEventListener('resize', () => { if (main.style.getPropertyValue('--left-w')) main.style.setProperty('--left-w', capW(parseFloat(main.style.getPropertyValue('--left-w'))) + 'px');
    if (tr.style.getPropertyValue('--top-h') && !tr.hidden) tr.style.setProperty('--top-h', capH(parseFloat(tr.style.getPropertyValue('--top-h'))) + 'px'); }); }
const setPanel = on => { const p = byId('ptSettings'); if (!p) return; p.hidden = !on; byId('setBtn').setAttribute('aria-expanded', String(on)); if (on) { const f = p.querySelector('input, select, button'); f && f.focus(); } else byId('setBtn').focus(); };
byId('setBtn').addEventListener('click', () => setPanel(byId('ptSettings') && byId('ptSettings').hidden));
document.addEventListener('click', e => { if (e.target.closest('#ptSettingsClose')) setPanel(false); });
document.addEventListener('keydown', e => { if (e.key === 'Escape' && byId('ptSettings') && !byId('ptSettings').hidden && !document.querySelector('dialog[open], .modal.in')) setPanel(false); });
byId('editBtn').addEventListener('click', () => PT.editProject());
byId('keysBtn').addEventListener('click', () => byId('keysDlg').showModal());
byId('keysClose').addEventListener('click', () => byId('keysDlg').close());
document.addEventListener('keydown', async e => {
  const typing = e.target.matches && (e.target.matches('input, select, textarea') || !!e.target.closest('#ptSettings, .modal, .selectize-control')), open = document.querySelector('dialog[open], .modal.in');
  if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) { e.preventDefault(); if (!open) { closeSugg(); if (await commit()) generate(); } return; }
  if (e.key === 's' && (e.ctrlKey || e.metaKey) && !e.shiftKey && !open) { e.preventDefault(); if (S.drawn) byId('dlSvg').click(); return; }   // Ctrl+S saves the figure, not a copy of the page, wherever the focus is
  if (open || typing) return;
  if (e.key === '/') { e.preventDefault(); (S.region ? byId('range') : byId('gene')).focus(); return; }
  if (e.key === '?') return byId('keysDlg').showModal();
  if (!e.ctrlKey && !e.metaKey && !e.altKey && !byId('samples').hidden) {   // samples: 1..9 switch one, A all, 0 none; with Shift, A and 0 are for the tracks
    const d = /^(?:Digit|Numpad)([0-9])$/.exec(e.code), n = P.samples.length, tr = !byId('tracks').hidden;
    if (d && e.shiftKey && d[1] === '0' && tr) { e.preventDefault(); return byId('clearAll').click(); }
    if (d && !e.shiftKey && d[1] === '0') { e.preventDefault(); return setAll(false); }
    if (d && !e.shiftKey && n <= 9 && +d[1] <= n) { e.preventDefault(); return hotSample(+d[1] - 1); }
    const kk = e.key.toLowerCase();   // letters by what is typed, not by where the key sits (French and German keyboards)
    if (kk === 'a' && e.shiftKey && tr) { e.preventDefault(); return byId('addAll').click(); }
    if (kk === 'a' && !e.shiftKey) { e.preventDefault(); return setAll(true); }
  }
  if ((e.key === '[' || e.key === ']') && S.hist.length) { const i = S.hist.indexOf(S.drawn), j = i < 0 ? 0 : i + (e.key === '[' ? 1 : -1); if (S.hist[j]) reopen(j); return; }   // the strip is newest first
});
document.documentElement.dataset.hover = 'select';   // every hover is light terracotta
{ let m = q.get('mode'); if (!m) { try { m = localStorage.getItem('pt-mode'); } catch (e) {} }
  if (m === 'dark') document.documentElement.setAttribute('data-pt-mode', 'dark'); }

/* ------------------------------------------------------------------ the app around the page
   PT is what the Shiny side (R/app-screen-server.R) and a test page call. */
const PT = {
  start(id) { return reset(id || 'app'); },
  // the Figure settings changed on the server: a figure drawn under other settings is another figure
  settings(key) { setSeq++; if (key === SETKEY) return; SETKEY = key; if (S) changed(() => { S.set = { key }; }); },
  editProject() { Shiny.setInputValue('open_settings', Date.now(), { priority: 'event' }); }
};
// Closed or reloaded, as opposed to a lost connection: the page says so on its way out. Over the open connection first (it gets
// there before the connection's own end), and by a plain request as well, for a browser that drops the first.
addEventListener('pagehide', ev => { if (ev.persisted) return;   // kept by the browser to come back to (Back): not a close
  try { if (window.Shiny && Shiny.setInputValue) Shiny.setInputValue('pt_bye', Date.now(), { priority: 'event' }); } catch (e) {}
  if (P && P.bye) { try { fetch(P.bye, { keepalive: true }); } catch (e) {} } });
if (window.Shiny) {
  Shiny.addCustomMessageHandler('pt_rpc_reply', m => { const done = rpcWait[m.id]; if (!done) return; delete rpcWait[m.id]; let j; try { j = JSON.parse(m.json); } catch (e) { j = { error: 'The answer could not be read.' }; } done(j); });
  Shiny.addCustomMessageHandler('pt_settings', m => PT.settings(m.key));
  jQuery(document).on('shiny:sessioninitialized', () => PT.start());
  // The connection is gone (R stopped, the laptop slept): what was asked is answered, later questions are answered at once, and the
  // page stays usable, as the figures drawn so far exist only here (Shiny's grey sheet over the page is switched off in main.css).
  jQuery(document).on('shiny:disconnected', () => { rpcGone = true; for (const id of Object.keys(rpcWait)) { rpcWait[id]({ error: GONE }); delete rpcWait[id]; }
    if (!document.querySelector('.modal.in')) byId('gone').hidden = false; });   // not when the app stopped on purpose (Edit project settings, Reset): its dialog says what to do
}
