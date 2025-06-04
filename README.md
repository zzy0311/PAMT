


# PAMT: Parallel Motif Transition Discovery

## Overview

This repository contains the implementation of our Parallel Motif Transition Discovery (PAMT) algorithm for temporal graphs. The algorithm employs a Temporal Zone Partitioning (TZP) strategy to divide the temporal graph into multiple independent partitions (growth zones). In each partition, the motif transitions are counted exactly while enabling high parallelism. Finally, the motifs are encoded using a deterministic relabeling scheme. This approach significantly improves efficiency and scalability without sacrificing accuracy.

## Directory Structure

```
.
├── Makefile
├── ptmt.cpp
├── ptmt.hpp
├── partition.hpp
├── encodeMotif.cpp
├── encodeMotif.hpp
└── README.md
```

## Compilation

To compile the project, run:
```bash
make
```
This command will generate the executable `PTMT`.

## Usage

Run the executable using the following command:
```bash
./PAMT input.txt Max_edge Max_memory consecutive omega
```
Where:
- **`input.txt`** is the input file containing the temporal edges.
- **`Max_edge`** is the maximum number of edges allowed per motif transition.
- **`Max_memory`** is the delta parameter used for partitioning.
- **`consecutive`** is a string flag (e.g., `"YES"`) specifying  consecutive transitions.
- **`omega`** is the temporal expansion factor used in the partitioning process.

For example:
```bash
./PAMT dataset.txt 3 600 YES 20
```

## Requirements

- C++11 (or later)
- OpenMP support
- Standard C++ libraries

