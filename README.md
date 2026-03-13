# 🚂 Tracks Game Solver

A Java Swing application that solves the **Tracks logic puzzle** using two algorithmically distinct approaches: **Backtracking DFS** and **Dynamic Programming**.

> Built as part of the Design and Analysis of Algorithms course — March 2026.

---

## 🧩 What is the Tracks Puzzle?

Tracks is a logic puzzle played on an **n × n grid**. The goal is to lay a continuous railway track from a start cell **A** to an end cell **B**, such that:
- Every step is orthogonally adjacent (no diagonals)
- No cell is visited more than once
- Every row and column contains **exactly as many track cells as its border clue**

---

## ⚙️ Algorithms Implemented

### 1. Backtracking DFS
A depth-first search with a **place → recurse → undo** cycle.

**Key features:**
- **Overshoot pruning** — prunes if any row/column exceeds its clue *(O(n) per call)*
- **Capacity pruning** — prunes if remaining unvisited cells can't satisfy a clue *(O(n²) per call)*
- **Merge sort candidate ordering** — neighbours sorted by Manhattan distance to end cell

| Metric | Value |
|---|---|
| Worst-case time | O(bᵏ) exponential |
| Space | O(N) — stack only |
| Completeness | ✅ Guaranteed |

---

### 2. Dynamic Programming (Bottom-Up, 3-Phase)

Builds a **layered state DAG**, eliminating redundant re-exploration via state de-duplication.

**Three phases:**
1. **Forward BFS** — enumerate all distinct reachable states layer by layer
2. **Backward propagation** — label which states can reach a valid solution
3. **Greedy extraction** — walk forward once along solvable states (zero backtracking)

**State key:** `(currentCell, rowUsed[], colUsed[], visitedMask)` hashed into a 64-bit key using multiplicative hashing. Fingerprint updated in **O(1)** per step via incremental XOR.

| Metric | Value |
|---|---|
| Total time | O(S · n) pseudo-polynomial |
| Space | O(S · n) — full memo table |
| Completeness | ✅ Guaranteed |
| Backtracking ops | 0 |

---

## 📊 Algorithm Comparison

| Feature | Backtracking | Dynamic Programming |
|---|---|---|
| Strategy | Try all paths, undo on failure | Store and reuse partial results |
| Memory | Low (stack only) | High (memo table) |
| Best for | Small, complex puzzles | Larger grids with many equivalent states |
| Backtracking ops | Many | **Zero** |

---

## 🎮 Application Features

- **Interactive player mode** with move validation
- **Computer turn** — solver picks the best next move
- **Hint system** — highlights the recommended next cell
- **Animated auto-solve** — step-by-step visualization (100ms delay)
- **5 seeded puzzle patterns** + unlimited random generation
- Grid sizes: **8×8 to 10×10**

---

## 🔗 Classical Problem Connections

| Classical Problem | Connection |
|---|---|
| Hamiltonian Path | Tracks BT is a constrained self-avoiding walk |
| N-Queens | Same increment/decrement counting pattern |
| Subset Sum B&B | Capacity pruning = proactive infeasibility bound |
| DAG Reachability | DP Phase 2 is exactly backward OR-propagation |
| 0/1 Knapsack | Same layered table, fixed-budget architecture |
| Bitmask Hamiltonian DP | Visited set as first-class state component |

---

## 👥 Team

| Name | Roll Number |
|---|---|
| Mithesh G S | CB.SC.U4CSE24715 |
| Naveen Velan S N | CB.SC.U4CSE24734 |
| Niranjan Reddy | CB.SC.U4CSE24735 |
| Sidambarisvar Balamurugan | CB.SC.U4CSE24751 |

---

## 🛠️ Tech Stack

- **Language:** Java 17
- **GUI:** Java Swing
- **Data structures:** `BitSet`, `HashMap`, `ArrayList`
- **Threading:** Solver on background thread; animation on Swing EDT
