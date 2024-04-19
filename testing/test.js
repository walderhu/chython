const { clean2d } = require('./index.js');

const smiles = "CC=CC";
const coordinates = clean2d(smiles);

console.log(coordinates);