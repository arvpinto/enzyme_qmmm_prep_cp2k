/* qm_viewer.js
 *
 * Interactive 3D view of the QM region of the tutorial (Section II), drawn with 3Dmol.js from
 * system_prep/hl_65atoms.mol2 (the MOL2 file written by PyMOL from frame_xray.pdb).
 *
 * Markup (see enzyme_qmmm_prep_II.md):
 *   <div class="qm-viewer" data-mol2="URL of hl_65atoms.mol2" data-fallback="assets/img/qm_region.png">
 *     <div class="qm-viewer__toolbar">  buttons with data-toggle="labels|distances|boundary|hydrogens|spin"
 *                                       and data-action="reset"
 *     <div class="qm-viewer__stage">    the 3Dmol canvas
 *     <p class="qm-viewer__info">       name, residue and serial number of the atom under the pointer
 *   </div>
 *
 * The data below come from the files of system_prep/ and are only used for annotations:
 *   SERIALS   serial number (1-based, as in CP2K and VMD) of each MOL2 atom, in MOL2 order, from amylase.prmtop
 *   BOUNDARY  the five QM atoms bonded to MM atoms (the &LINK sections of forceeval_qmmm.cp2k.inc) and the
 *             position of their MM partner in frame_xray.pdb (the MOL2 coordinates are those of that PDB)
 *   DISTANCES the two distances of the reaction followed in the tutorial
 */
(function () {
  "use strict";

  var SERIALS = [
    3022, 3025, 3026, 3027, 3023, 3024, 3601, 3604, 3605, 3606, 3602, 3603, 3607, 4650, 4653, 4654,
    4655, 4651, 4652, 7678, 7679, 7681, 7683, 7685, 7687, 7689, 7699, 7691, 7698, 7693, 7696, 7680,
    7682, 7686, 7690, 7692, 7684, 7688, 7695, 7694, 7697, 7700, 7702, 7704, 7706, 7708, 7710, 7720,
    7712, 7719, 7714, 7717, 7701, 7703, 7707, 7711, 7713, 7705, 7709, 7716, 7715, 7718, 7783, 7784,
    7785
  ];

  var BOUNDARY = [
    { qm: 3022, mm: { serial: 3020, name: "CA", res: "ASP196", x: 48.394, y: 39.140, z: 46.643 } },
    { qm: 3601, mm: { serial: 3598, name: "CB", res: "GLH232", x: 44.314, y: 40.935, z: 46.063 } },
    { qm: 4650, mm: { serial: 4648, name: "CA", res: "ASP299", x: 46.303, y: 46.178, z: 37.530 } },
    { qm: 7678, mm: { serial: 7668, name: "C4", res: "4GA499", x: 41.471, y: 37.165, z: 37.451 } },
    { qm: 7720, mm: { serial: 7721, name: "C1", res: "4GA502", x: 50.853, y: 40.263, z: 37.295 } }
  ];

  var DISTANCES = [
    { a: 3026, b: 7700, what: "Asp196 OD1 - C1 (nucleophilic attack)" },
    { a: 3607, b: 7699, what: "Glu232 HE2 - O4 (proton transfer)" }
  ];

  var RES_NOTES = {
    ASP196: "ASP196", GLH232: "GLH232", ASP299: "ASP299",
    "4GA499": "4GA499", "4GA500": "4GA500 (+1)", "4GA501": "4GA501 (-1)", WAT511: "WAT511"
  };

  var COLORS = {
    bg: "#0e1014",
    label: "#f3f5f8",
    labelBg: "#1a1d24",
    boundary: "#ff9d42",
    distance: "#f2d24b",
    hover: "#65cbf3"
  };

  var STICK = { radius: 0.13, colorscheme: "cyanCarbon" };
  var BALL = { scale: 0.2, colorscheme: "cyanCarbon" };

  // 3Dmol.js keeps the SYBYL type of MOL2 atoms but not their names and residues, so they are read here.
  // Fields of an ATOM line: id, name, x, y, z, type, substructure id, substructure name, charge.
  // PyMOL writes the substructure name as RESNAME+RESID (e.g. ASP196).
  function parseMol2Atoms(text) {
    var block = text.split("@<TRIPOS>ATOM")[1];
    if (!block) return [];
    block = block.split("@<TRIPOS>")[0];
    return block.split(/\r?\n/).map(function (l) { return l.trim(); }).filter(Boolean).map(function (l) {
      var f = l.split(/\s+/);
      return { id: parseInt(f[0], 10), name: f[1], res: f[7] || "" };
    });
  }

  function residueOf(atom) {
    return atom.qmRes || "";
  }

  function nameOf(atom) {
    return atom.qmName || atom.atom;
  }

  function fallback(root, message) {
    var stage = root.querySelector(".qm-viewer__stage");
    var info = root.querySelector(".qm-viewer__info");
    var bar = root.querySelector(".qm-viewer__toolbar");
    if (bar) bar.hidden = true;
    var src = root.getAttribute("data-fallback");
    if (stage) {
      stage.classList.add("qm-viewer__stage--static");
      stage.innerHTML = src ? '<img src="' + src + '" alt="QM region of the example (static image)">' : "";
    }
    if (info) info.textContent = message;
  }

  function init(root) {
    var stage = root.querySelector(".qm-viewer__stage");
    var info = root.querySelector(".qm-viewer__info");
    var defaultInfo = info ? info.textContent : "";
    var url = root.getAttribute("data-mol2");

    if (typeof window.$3Dmol === "undefined") {
      fallback(root, "The 3D viewer could not be loaded; this is a static image of the QM region.");
      return;
    }

    fetch(url)
      .then(function (r) {
        if (!r.ok) throw new Error("HTTP " + r.status);
        return r.text();
      })
      .then(function (text) { build(text); })
      .catch(function () {
        fallback(root, "The MOL2 file could not be loaded; this is a static image of the QM region.");
      });

    function build(mol2) {
      var viewer = window.$3Dmol.createViewer(stage, { backgroundColor: COLORS.bg, antialias: true });
      var model = viewer.addModel(mol2, "mol2");
      var atoms = model.selectedAtoms({});
      var meta = parseMol2Atoms(mol2);
      var bySerial = {};
      var withSerials = atoms.length === SERIALS.length;
      atoms.forEach(function (atom, i) {
        if (meta.length === atoms.length) {
          atom.qmName = meta[i].name;
          atom.qmRes = meta[i].res;
        }
        atom.topSerial = withSerials ? SERIALS[i] : null;
        if (atom.topSerial) bySerial[atom.topSerial] = atom;
      });

      var state = { labels: true, distances: true, boundary: true, hydrogens: true, spin: false };
      var hoverLabel = null;

      function applyStyle() {
        model.setStyle({}, { stick: STICK, sphere: BALL });
        if (!state.hydrogens) model.setStyle({ elem: "H" }, {});
      }

      function drawAnnotations() {
        viewer.removeAllLabels();
        viewer.removeAllShapes();
        hoverLabel = null;

        if (state.labels) {
          var groups = {};
          atoms.forEach(function (a) {
            if (a.elem === "H") return;
            var r = residueOf(a);
            (groups[r] = groups[r] || []).push(a);
          });
          Object.keys(groups).forEach(function (r) {
            var g = groups[r], c = { x: 0, y: 0, z: 0 };
            g.forEach(function (a) { c.x += a.x; c.y += a.y; c.z += a.z; });
            c.x /= g.length; c.y /= g.length; c.z /= g.length;
            viewer.addLabel(RES_NOTES[r] || r, {
              position: c, fontSize: 13, fontColor: COLORS.label, backgroundColor: COLORS.labelBg,
              backgroundOpacity: 0.8, borderThickness: 0, inFront: true, alignment: "center"
            });
          });
        }

        if (state.distances && withSerials) {
          DISTANCES.forEach(function (d) {
            var a = bySerial[d.a], b = bySerial[d.b];
            if (!a || !b) return;
            var dist = Math.sqrt(Math.pow(a.x - b.x, 2) + Math.pow(a.y - b.y, 2) + Math.pow(a.z - b.z, 2));
            viewer.addCylinder({
              start: { x: a.x, y: a.y, z: a.z }, end: { x: b.x, y: b.y, z: b.z },
              radius: 0.05, color: COLORS.distance, dashed: true, dashLength: 0.18, gapLength: 0.14,
              fromCap: 1, toCap: 1
            });
            viewer.addLabel(dist.toFixed(2) + " \u00c5", {
              position: { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2, z: (a.z + b.z) / 2 },
              fontSize: 12, fontColor: "#14130f", backgroundColor: COLORS.distance, backgroundOpacity: 0.95,
              borderThickness: 0, inFront: true, alignment: "center"
            });
          });
        }

        if (state.boundary && withSerials) {
          BOUNDARY.forEach(function (bd) {
            var q = bySerial[bd.qm];
            if (!q) return;
            var m = bd.mm;
            viewer.addCylinder({
              start: { x: q.x, y: q.y, z: q.z }, end: { x: m.x, y: m.y, z: m.z },
              radius: 0.1, color: COLORS.boundary, fromCap: 1, toCap: 1
            });
            viewer.addSphere({ center: { x: m.x, y: m.y, z: m.z }, radius: 0.32, color: COLORS.boundary });
          });
        }
        viewer.render();
      }

      function describe(atom) {
        var s = nameOf(atom) + " \u00b7 " + residueOf(atom);
        if (atom.topSerial) s += " \u00b7 serial " + atom.topSerial;
        var bd = BOUNDARY.filter(function (b) { return b.qm === atom.topSerial; })[0];
        if (bd) s += " \u00b7 bonded to MM atom " + bd.mm.name + " " + bd.mm.res + " (serial " + bd.mm.serial + ")";
        DISTANCES.forEach(function (d) {
          if (d.a === atom.topSerial || d.b === atom.topSerial) s += " \u00b7 " + d.what;
        });
        return s;
      }

      function showAtom(atom) {
        if (info) info.textContent = describe(atom);
        if (hoverLabel) viewer.removeLabel(hoverLabel);
        hoverLabel = viewer.addLabel(nameOf(atom) + (atom.topSerial ? " " + atom.topSerial : ""), {
          position: { x: atom.x, y: atom.y, z: atom.z }, fontSize: 12, fontColor: "#001520",
          backgroundColor: COLORS.hover, backgroundOpacity: 0.95, borderThickness: 0, inFront: true
        });
        viewer.render();
      }

      function hideAtom() {
        if (hoverLabel) { viewer.removeLabel(hoverLabel); hoverLabel = null; viewer.render(); }
        if (info) info.textContent = defaultInfo;
      }

      viewer.setHoverable({}, true, function (atom) { showAtom(atom); }, function () { hideAtom(); });
      viewer.setClickable({}, true, function (atom) { showAtom(atom); });
      viewer.setHoverDuration(80);

      function fit() {
        viewer.zoomTo();
        viewer.zoom(stage.clientWidth < 560 ? 0.85 : 1.05);
        viewer.render();
      }

      applyStyle();
      drawAnnotations();
      fit();

      root.querySelectorAll("[data-toggle]").forEach(function (btn) {
        var key = btn.getAttribute("data-toggle");
        btn.setAttribute("aria-pressed", String(state[key]));
        btn.addEventListener("click", function () {
          state[key] = !state[key];
          btn.setAttribute("aria-pressed", String(state[key]));
          if (key === "spin") {
            viewer.spin(state.spin ? "y" : false, 0.6);
          } else if (key === "hydrogens") {
            applyStyle();
            viewer.render();
          } else {
            drawAnnotations();
          }
        });
      });
      var reset = root.querySelector('[data-action="reset"]');
      if (reset) {
        reset.addEventListener("click", fit);
      }
      window.addEventListener("resize", function () { viewer.resize(); viewer.render(); });
      root.classList.add("qm-viewer--ready");
    }
  }

  function start() {
    document.querySelectorAll(".qm-viewer").forEach(init);
  }
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", start);
  } else {
    start();
  }
})();
