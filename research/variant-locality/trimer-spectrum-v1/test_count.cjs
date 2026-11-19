"use strict";

const assert = require("node:assert/strict");
const {test} = require("node:test");
const {countWithFixedEnds, reverseComplement, spectrumSummary} = require("./count.cjs");

function spectrumKey(sequence) {
    const trimers = [];
    for (let i = 0; i < sequence.length - 2; i++) {
        trimers.push(sequence.slice(i, i + 3));
    }
    return JSON.stringify([sequence.slice(0, 2), sequence.slice(-2), trimers.sort()]);
}

function allSequences(alphabet, length) {
    if (length === 0) return [""];
    return allSequences(alphabet, length - 1).flatMap(s => [...alphabet].map(b => s + b));
}

test("exhaustive short strings agree with direct grouping", () => {
    for (let length = 3; length <= 7; length++) {
        const classes = new Map();
        for (const s of allSequences("ACGT", length)) {
            const key = spectrumKey(s);
            const group = classes.get(key);
            if (group) group.size++;
            else classes.set(key, {sequence: s, size: 1});
        }
        for (const {sequence, size} of classes.values()) {
            assert.equal(countWithFixedEnds(sequence), BigInt(size), sequence);
        }
    }
});

test("open and closed walks account for loops and parallel edges", () => {
    assert.equal(countWithFixedEnds("AAAAAA"), 1n);
    assert.equal(countWithFixedEnds("AAACAA"), 2n);
    assert.equal(countWithFixedEnds("AACAAGAA"), 2n);
    assert.equal(countWithFixedEnds("AACAAGAAC"), 2n);
    assert.equal(countWithFixedEnds("ACACAC"), 1n);
    assert.throws(() => countWithFixedEnds("AC"), TypeError);
    assert.throws(() => countWithFixedEnds("ACNT"), TypeError);
    assert.throws(() => countWithFixedEnds("ACGT\n"), TypeError);
});

test("distinct two-edit strings can preserve oriented trimer counts", () => {
    const left = "ACGT".repeat(24) + "AC";
    const right = "TGCA".repeat(23) + "TG";
    const a = left + "ACGACTAC" + right;
    const b = left + "ACTACGAC" + right;
    assert.equal(a.length, 200);
    assert.equal(b.length, 200);
    const differences = [...a].flatMap((base, i) => base === b[i] ? [] : [i]);
    assert.deepEqual(differences, [100, 103]);
    assert.notEqual(a, reverseComplement(b));
    assert.equal(spectrumKey(a), spectrumKey(b));
    const stats = spectrumSummary(a);
    assert.ok(stats.strings >= 2n);
    assert.ok(stats.reverseComplementOrbitsLowerBound >= 2n);
});
