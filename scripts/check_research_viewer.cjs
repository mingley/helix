#!/usr/bin/env node
// Optional interaction smoke check using the generated artifact, without network.
// This checks DOM behavior; it is not a browser layout or accessibility audit.
"use strict";
const fs = require("node:fs"), vm = require("node:vm"), assert = require("node:assert/strict");
const path = process.argv[2];
if (!path) throw new Error("usage: node scripts/check_research_viewer.cjs <research.html>");
const html = fs.readFileSync(path, "utf8");
const data = JSON.parse(html.match(/<script id="data" type="application\/json">([\s\S]*?)<\/script>/)[1]);
const scripts = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)];
assert.equal(scripts.length, 1);
class Element {
  constructor(tag) { this.tag = tag; this.children = []; this.attrs = {}; this.events = {}; this.style = {}; this.value = ""; this.textContent = ""; this.namespaceURI = "http://www.w3.org/2000/svg"; this.dataset = {}; }
  appendChild(child) { this.children.push(child); return child; }
  setAttribute(key, value) { this.attrs[key] = value; }
  addEventListener(key, callback) { this.events[key] = callback; }
  replaceChildren() { this.children = []; }
}
const elements = new Map();
const document = {
  getElementById(id) { if (!elements.has(id)) elements.set(id, new Element(id)); return elements.get(id); },
  createElement(tag) { return new Element(tag); },
  createElementNS(_, tag) { return new Element(tag); }
};
document.getElementById("data").textContent = JSON.stringify(data);
vm.runInNewContext(scripts[0][1], { document }, { timeout: 5000 });
const rows = document.getElementById("rows"), locus = document.getElementById("locus");
assert.equal(rows.children.length, data.comparisons.length);
assert.equal(locus.children.filter(c => c.tag === "rect").length, 6);
assert.equal(locus.children.filter(c => c.tag === "circle").length, data.comparisons.length);
const query = document.getElementById("query");
query.value = "synonymous";
query.events.input();
const expected = data.comparisons.filter(c => JSON.stringify(c).toLowerCase().includes("synonymous"));
assert.equal(rows.children.length, expected.length);
assert.ok(expected.length > 0);
rows.children[0].events.click();
assert.equal(JSON.parse(document.getElementById("detail").textContent).id, expected[0].id);
query.value = "no_matching_accession_123";
query.events.input();
assert.equal(rows.children.length, 0);
query.value = "";
query.events.input();
assert.equal(rows.children.length, data.comparisons.length);
assert.ok(document.getElementById("atlas").children.length > 0);
assert.ok(document.getElementById("effects").children.length > 0);
assert.ok(document.getElementById("assay").children.length > 0);
console.log("Research viewer interactions verified: locus, rows, filter, selection and restored table.");
