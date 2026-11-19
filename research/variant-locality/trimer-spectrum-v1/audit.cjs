#!/usr/bin/env node
"use strict";

const {createHash} = require("node:crypto");
const {readFileSync} = require("node:fs");
const {resolve} = require("node:path");
const {reverseComplement, spectrumSummary} = require("./count.cjs");

const ROOT = resolve(__dirname, "../../..");
const SOURCE = "research/context-source-alternatives";
const CROSSCHECK = SOURCE + "/design-coordinate-crosscheck.json";
const PANEL = "datasets/context-panel/panel.json";
const PANEL_SHA256 = "fce6b2a1aaf5aa4d530a2e298dccb840bc7c594a05b6831fdb0c6c0f9e2a9be3";

function sha256(bytes) {
    return createHash("sha256").update(bytes).digest("hex");
}

function read(relativePath) {
    return readFileSync(resolve(ROOT, relativePath));
}

function readExamples() {
    const {checks} = JSON.parse(read(CROSSCHECK));
    return checks.map((item, index) => {
        const source = SOURCE + "/design-coordinate-example-" + (index + 1) + ".json";
        const snapshot = JSON.parse(read(source));
        if (snapshot.genome !== "hg38" || snapshot.chrom !== item.chromosome ||
            snapshot.start !== item.start || snapshot.end !== item.end ||
            !["+", "-"].includes(item.strand) ||
            item.zero_based_half_open_reference_matches !== true ||
            typeof snapshot.dna !== "string" || !/^[ACGTacgt]+$/.test(snapshot.dna)) {
            throw new Error("reference metadata does not match for " + item.id);
        }
        const forward = snapshot.dna.toUpperCase();
        const sequence = item.strand === "+" ? forward : reverseComplement(forward);
        if (sequence.length !== 200 || item.end - item.start !== 200 ||
            sha256(sequence) !== item.core_sha256) {
            throw new Error("reference sequence does not match for " + item.id);
        }
        const stats = spectrumSummary(sequence);
        return {
            id: item.id,
            reference_snapshot: source,
            sequence_sha256: item.core_sha256,
            prefix: stats.prefix,
            suffix: stats.suffix,
            fixed_end_oriented_trimer_class_size: stats.strings,
            reverse_complement_orbits_lower_bound: stats.reverseComplementOrbitsLowerBound,
        };
    });
}

function distribution(values) {
    values.sort((a, b) => a < b ? -1 : Number(a > b));
    const quantile = p => values[Math.ceil(p * values.length) - 1];
    return {
        minimum: values[0],
        percentile_1: quantile(1 / 100),
        percentile_25: quantile(1 / 4),
        percentile_50: quantile(1 / 2),
        percentile_75: quantile(3 / 4),
        percentile_99: quantile(99 / 100),
        maximum: values.at(-1),
    };
}

function summarizePanel() {
    const bytes = read(PANEL);
    if (sha256(bytes) !== PANEL_SHA256) throw new Error("unexpected panel SHA-256");
    const {records} = JSON.parse(bytes);
    if (!Array.isArray(records) || records.length !== 55276) {
        throw new Error("unexpected panel size");
    }
    const strings = [];
    const orbitBounds = [];
    const distinctInputs = new Set();
    let incompatibleEnds = 0;
    let atLeast1e50 = 0;
    let atLeast1e70 = 0;
    for (const {sequence} of records) {
        if (typeof sequence !== "string" || sequence.length !== 200) {
            throw new Error("unexpected panel sequence length");
        }
        distinctInputs.add(sequence);
        const stats = spectrumSummary(sequence);
        const {strings: count, reverseComplementOrbitsLowerBound: lowerBound} = stats;
        strings.push(count);
        orbitBounds.push(lowerBound);
        if (stats.prefix !== reverseComplement(stats.suffix)) incompatibleEnds++;
        if (lowerBound >= 10n ** 50n) atLeast1e50++;
        if (lowerBound >= 10n ** 70n) atLeast1e70++;
    }
    return {
        panel: PANEL,
        panel_sha256: PANEL_SHA256,
        records: records.length,
        distinct_input_sequences: distinctInputs.size,
        records_where_terminal_dinucleotides_exclude_reverse_complement_pairs: incompatibleEnds,
        records_with_at_least_1e50_reverse_complement_orbits: atLeast1e50,
        records_with_at_least_1e70_reverse_complement_orbits: atLeast1e70,
        fixed_end_oriented_trimer_class_sizes: distribution(strings),
        reverse_complement_orbit_lower_bounds: distribution(orbitBounds),
    };
}

function serialize(report) {
    return JSON.stringify(report, (_, value) => typeof value === "bigint" ? value.toString() : value, 2) + "\n";
}

function main(mode) {
    if (mode === "--examples") return serialize({examples: readExamples()});
    if (mode === "--panel") return serialize(summarizePanel());
    if (mode === "--verify") {
        const filename = "research/variant-locality/trimer-spectrum-v1/examples.json";
        if (serialize({examples: readExamples()}) !== read(filename).toString("utf8")) {
            throw new Error("replayed examples differ from " + filename);
        }
        return "Both reference examples match the retained exact counts.\n";
    }
    throw new Error("usage: node audit.cjs --examples | --panel | --verify");
}

if (require.main === module) process.stdout.write(main(process.argv[2]));

module.exports = {main};
