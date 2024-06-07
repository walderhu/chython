getMolecularFormula(data = null) {
    let molecularFormula = '';
    let counts = new Map();
    let graph = data === null ? this.graph : new Graph(data, this.opts.isomeric);
    // Initialize element count
    for (var i = 0; i < graph.vertices.length; i++) {
        let atom = graph.vertices[i].value;
        if (counts.has(atom.element)) {
            counts.set(atom.element, counts.get(atom.element) + 1);
        } else {
            counts.set(atom.element, 1);
        }
        // Hydrogens attached to a chiral center were added as vertices,
        // those in non chiral brackets are added here
        if (atom.bracket && !atom.bracket.chirality) {
            if (counts.has('H')) {
                counts.set('H', counts.get('H') + atom.bracket.hcount);
            } else {
                counts.set('H', atom.bracket.hcount);
            }
        }
        // Add the implicit hydrogens according to valency, exclude
        // bracket atoms as they were handled and always have the number
        // of hydrogens specified explicitly
        if (!atom.bracket) {
            let nHydrogens = Atom.maxBonds[atom.element] - atom.bondCount;

            if (atom.isPartOfAromaticRing) {
                nHydrogens--;
            }

            if (counts.has('H')) {
                counts.set('H', counts.get('H') + nHydrogens);
            } else {
                counts.set('H', nHydrogens);
            }
        }
    }
    if (counts.has('C')) {
        let count = counts.get('C');
        molecularFormula += 'C' + (count > 1 ? count : '');
        counts.delete('C');
    }
    if (counts.has('H')) {
        let count = counts.get('H');
        molecularFormula += 'H' + (count > 1 ? count : '');
        counts.delete('H');
    }
    let elements = Object.keys(Atom.atomicNumbers).sort();
    elements.map(e => {
        if (counts.has(e)) {
            let count = counts.get(e);
            molecularFormula += e + (count > 1 ? count : '');
        }
    });
    return molecularFormula;
}