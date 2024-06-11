


initHydrogens() {
    if (!this.opts.explicitHydrogens)
        for (var t = 0; t < this.graph.vertices.length; t++) {
            let e = this.graph.vertices[t];
            if ("H" !== e.value.element) continue;
            let i = this.graph.vertices[e.neighbours[0]];
            i.value.hasHydrogen = !0, (!i.value.isStereoCenter ||
                i.value.rings.length < 2 && !i.value.bridgedRing ||
                i.value.bridgedRing && i.value.originalRings.length < 2) &&
                (e.value.isDrawn = !1)
        }
}

