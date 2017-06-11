#!/usr/bin/env python3.6
# Author: Marek Marusic
# Purpose: Project for BIF 2017

import datetime
import sys
from operator import itemgetter

seqid, source, typeMethod, start, end, score, strand, phase = range(8)


def printUsage():
    print("Usage: ./resolve_overlaps.py fileName")


def findBestNonOverlaps(inputFile):
    sequenceRegions = readSequences(inputFile)
    print("##gff-version 3\n##date " + datetime.datetime.now().strftime("%Y-%m-%d"))
    for key in sequenceRegions:
        print(key, file=sys.stderr)
        sequences = sorted(sequenceRegions[key], key=itemgetter(end, start))
        bestNonOverlappedIndexes = findBestSolutionInSeq(sequences)

        for index in reversed(bestNonOverlappedIndexes):
            printSequence(sequences[index])


def readSequences(file):
    sequenceRegions = {}
    with open(file) as openedFile:
        for line in openedFile:
            if isComment(line):
                pass
            else:
                parts = line.split()
                key = createIndex(parts)
                if sequenceRegions.get(key) is None:
                    sequenceRegions[key] = []
                    sequenceRegions[key].append(createEmptySequence(parts))
                sequenceRegions[key].append(createSequenceFrom(parts))
    return sequenceRegions


def isComment(line):
    return line.startswith("#")


def createIndex(parts):
    return parts[0] + parts[2] + parts[6]


def createEmptySequence(parts):
    return [parts[0], "##", parts[2], 0, 0, 0, parts[6], "##"]


def createSequenceFrom(parts):
    return [parts[seqid], parts[source], parts[typeMethod], int(parts[start]), int(parts[end]), int(parts[score]), parts[strand], parts[phase]]


def findBestSolutionInSeq(sequences):
    # Algoritmus Weighted Interval Scheduling: Bottom-Up
    closestNonOverlap = [0] * len(sequences)
    memory = [0] * len(sequences)
    solution = [0] * len(sequences)
    resultIndexes = []

    for i in range(1, len(sequences)):
        # Little optimization don't look for closest nonoverlapped if there isn't any
        if sequences[i][start] > sequences[1][end]:
            for j in reversed(range(1, i)):
                if sequences[i][start] > sequences[j][end]:
                    closestNonOverlap[i] = j
                    break

        if sequences[i][score] + memory[closestNonOverlap[i]] > memory[i-1]:
            memory[i] = sequences[i][score] + memory[closestNonOverlap[i]]
            solution[i] = closestNonOverlap[i]
        else:
            memory[i] = memory[i-1]
            solution[i] = i-1

    i = len(sequences) - 1
    while i > 0:
        if solution[i] == closestNonOverlap[i]:
            resultIndexes.append(i)
        i = solution[i]

    print("Result = ", memory[len(sequences)-1], file=sys.stderr)
    return resultIndexes


def printSequence(sequence):
    output = ""
    for part in sequence:
        output += str(part) + "\t"
    print(output[:-1])


# Main program start
if len(sys.argv) < 2:
    printUsage()
else:
    findBestNonOverlaps(sys.argv[1])









