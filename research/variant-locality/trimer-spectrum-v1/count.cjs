#!/usr/bin/env node
"use strict";

const BASES = "ACGT";
const CONTEXT_COUNT = 16;

function contextCode(sequence, offset) {
    return 4 * BASES.indexOf(sequence[offset]) + BASES.indexOf(sequence[offset + 1]);
}

function determinant(matrix) {
    // Bareiss elimination keeps all intermediate operations integral.
    const n = matrix.length;
    if (n === 0) return 1n;
    let previous = 1n;
    let sign = 1n;
    for (let column = 0; column < n - 1; column++) {
        if (matrix[column][column] === 0n) {
            const swap = matrix.findIndex((row, i) => i > column && row[column] !== 0n);
            if (swap === -1) return 0n;
            [matrix[column], matrix[swap]] = [matrix[swap], matrix[column]];
            sign = -sign;
        }
        const pivot = matrix[column][column];
        for (let row = column + 1; row < n; row++) {
            for (let j = column + 1; j < n; j++) {
                matrix[row][j] =
                    (pivot * matrix[row][j] - matrix[row][column] * matrix[column][j]) / previous;
            }
            matrix[row][column] = 0n;
        }
        previous = pivot;
    }
    return sign * matrix[n - 1][n - 1];
}

const factorials = [1n];
function factorial(n) {
    while (n >= factorials.length) {
        const next = factorials.length;
        factorials.push(factorials[next - 1] * BigInt(next));
    }
    return factorials[n];
}

function countWithFixedEnds(sequence) {
    if (typeof sequence !== "string" || sequence.length < 3 || /[^ACGT]/.test(sequence)) {
        throw new TypeError("expected at least three uppercase A/C/G/T bases");
    }
    const edges = Array.from({length: CONTEXT_COUNT}, () => Array(CONTEXT_COUNT).fill(0));
    const out = Array(CONTEXT_COUNT).fill(0);
    const active = new Set();
    let vertex = contextCode(sequence, 0);
    active.add(vertex);
    for (let i = 1; i < sequence.length - 1; i++) {
        const next = contextCode(sequence, i);
        edges[vertex][next]++;
        out[vertex]++;
        active.add(next);
        vertex = next;
    }
    const root = vertex;

    // The root is the final two bases. Self loops participate in the ordering
    // of outgoing edges, but they cannot be last exits in a tree to the root.
    const nonRoot = [...active].filter(v => v !== root);
    const minor = nonRoot.map(u => nonRoot.map(v =>
        BigInt(u === v ? out[u] - edges[u][u] : -edges[u][v])
    ));
    let numerator = determinant(minor);
    let denominator = 1n;
    for (const u of active) {
        numerator *= factorial(out[u] - (u === root ? 0 : 1));
        for (const multiplicity of edges[u]) {
            if (multiplicity > 1) denominator *= factorial(multiplicity);
        }
    }
    if (numerator <= 0n || numerator % denominator !== 0n) {
        throw new Error("invalid Euler trail count for a witnessed spectrum");
    }
    return numerator / denominator;
}

const COMPLEMENT = {A: "T", C: "G", G: "C", T: "A"};
function reverseComplement(sequence) {
    return [...sequence].reverse().map(base => COMPLEMENT[base]).join("");
}

function spectrumSummary(sequence) {
    const strings = countWithFixedEnds(sequence);
    const prefix = sequence.slice(0, 2);
    const suffix = sequence.slice(-2);
    return {
        prefix,
        suffix,
        strings,
        reverseComplementOrbitsLowerBound:
            prefix === reverseComplement(suffix) ? (strings + 1n) / 2n : strings,
    };
}

module.exports = {countWithFixedEnds, reverseComplement, spectrumSummary};
