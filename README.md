# Enhanced Hashing for Exact String Matching

C implementations of exact online string matching algorithms using enhanced hashing techniques, based on the paper:

**Improving String Matching through Enhanced Hashing**  
Francesco Pio Marino and Jorma Tarhio

## Overview

This repository contains C implementations of hashing-based exact string matching algorithms. The goal is to improve search speed by replacing or extending traditional hash functions with faster software and hardware-assisted hashing methods.

Implemented families include:

- **EPSM variants**
- **WFRq variants**
- **SSB variants**
- Hashing methods based on:
  - rolling hashes
  - CRC32
  - CRC32C on 64-bit keys
  - X-Hash
  - M-Hash using carry-less multiplication
  - SSEF-style hashing

## Algorithms

The repository includes variants of:

- `EPSM`
- `EPSM16C4`
- `EPSM16C8`
- `EPSM08M`
- `WFR`
- `WFRqX`
- `WFRqM`
- `WFRqC4`
- `WFRqC8`
- `WFR16S`
- `SSB13`
- `SSB08X`
- `SSB16S`
- `SSBC8`

Additional baseline hash variants may include:

- `WFRCL`
- `WFRMul`
- `WFRSM`
- `WFRWy`
- `WFRXXH3`

## Requirements

A C compiler with optimization support is required.

Recommended:

```bash
gcc
clang
make
```

## Running Experiments

All algorithms are executed through the SMART framework.

After compilation, experiments can be launched using:

```bash
./bin/smart run -plen 8 32768 \
  -text englishTexts \
  -pre \
  -runs 20 \
  -all \
  -ts 100 \
  -tb 500000 \
  -occ
