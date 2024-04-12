const DrawerBase = require('./DrawerBase');
const Parser = require('./Parser');

function clean2d(smiles) {
    const drawer = new DrawerBase({});
    const parsed = Parser.parse(smiles);
    drawer.initDraw(parsed, 'light', false);
    drawer.processGraph();

    let vertices = drawer.graph.vertices;
    let xy = Array();
    for (let i = 0; i < vertices.length; i++) {
        let position = vertices[i].position;
        xy.push([position.x, position.y]);
    }
    return xy;
}

module.exports = { clean2d };