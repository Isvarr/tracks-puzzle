import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.util.List;
import java.util.Timer;


public class TracksGame extends JFrame {

    private int gridSize = 8;
    private static final int CELL_SIZE = 50;

    private GamePanel gamePanel;
    private SidePanel sidePanel;
    private GameMode  currentMode = GameMode.PATTERN1;

    enum GameMode { PATTERN1, PATTERN2, PATTERN3, PATTERN4, PATTERN5, RANDOM }

    public TracksGame() {
        setTitle("Tracks Lab  —  Pure DP Solver");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLayout(new BorderLayout());
        gamePanel = new GamePanel();
        sidePanel = new SidePanel(gamePanel);
        add(gamePanel, BorderLayout.CENTER);
        add(sidePanel, BorderLayout.EAST);
        pack();
        setLocationRelativeTo(null);
        setResizable(false);
    }

    public static void main(String[] args) {
        try { UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName()); }
        catch (Exception ignored) {}
        SwingUtilities.invokeLater(() -> new TracksGame().setVisible(true));
    }

    // ════════════════════════════════════════════════════════════════
    //  CELL
    // ════════════════════════════════════════════════════════════════

    static class Cell {
        int r, c;
        boolean hasTrack  = false;
        boolean isCrossed = false;
        boolean isHint    = false;
        Cell    trackParent;

        Cell(int r, int c) { this.r = r; this.c = c; }

        @Override public String toString() { return "(" + r + "," + c + ")"; }
    }

    // ════════════════════════════════════════════════════════════════
    //  DP STATE  — one node in the forward BFS DAG
    // ════════════════════════════════════════════════════════════════

    static class DPState {
        final int     cellId;      // r*gridSize + c
        final int[]   rowUsed;     // track cells placed per row so far
        final int[]   colUsed;     // track cells placed per col so far
        final long    visitedFP;   // fingerprint of the visited-cell set
        final BitSet  visited;     // actual visited set (for child generation)
        final long    key;         // de-duplication hash

        boolean solvable = false;  // set during Phase 2

        DPState(int cellId, int[] rowUsed, int[] colUsed,
                BitSet visited, long visitedFP, long key) {
            this.cellId    = cellId;
            this.rowUsed   = rowUsed;
            this.colUsed   = colUsed;
            this.visited   = visited;
            this.visitedFP = visitedFP;
            this.key       = key;
        }
    }

    // ════════════════════════════════════════════════════════════════
    //  GAME PANEL
    // ════════════════════════════════════════════════════════════════

    class GamePanel extends JPanel {

        private Cell[][]   grid;
        private int[]      rowClues, colClues;
        private Cell       startNode, endNode, currentHead;

        private boolean    isPlayerTurn     = true;
        private boolean    gameOver         = false;
        private boolean    autoSolving      = false;
        private boolean    computerThinking = false;
        private List<Cell> solutionPath     = null;

        private int lastDpStates = 0;
        private int lastDpGoals  = 0;

        // ── constructor ──────────────────────────────────────────

        GamePanel() {
            setBackground(Color.WHITE);
            initGame(GameMode.PATTERN1);

            addMouseListener(new MouseAdapter() {
                @Override public void mouseClicked(MouseEvent e) {
                    if (gameOver || !isPlayerTurn || autoSolving || computerThinking) return;

                    int SX = 50, SY = 50;
                    int mc = (e.getX() - SX) / CELL_SIZE;
                    int mr = (e.getY() - SY) / CELL_SIZE;

                    if (mr < 0 || mc < 0 || mr >= gridSize || mc >= gridSize) return;

                    if (SwingUtilities.isRightMouseButton(e)) {
                        if (!grid[mr][mc].hasTrack && grid[mr][mc] != currentHead) {
                            grid[mr][mc].isCrossed = !grid[mr][mc].isCrossed;
                            repaint();
                        }
                        return;
                    }

                    if (SwingUtilities.isLeftMouseButton(e) && isValidMove(mr, mc)) {
                        clearHints();
                        makeMove(mr, mc);
                        repaint();
                        checkGameStatus();

                        if (!gameOver) {
                            computerThinking = true;
                            isPlayerTurn = false;
                            sidePanel.updateStatus("Computer thinking…");
                            new Timer().schedule(new TimerTask() {
                                @Override public void run() { computerTurn(); }
                            }, 50);
                        }
                    }
                }
            });
        }

        // ── init ────────────────────────────────────────────────

        void initGame(GameMode mode) {
            currentMode = mode;
            switch (mode) {
                case PATTERN1: gridSize = 8;  break;
                case PATTERN2: gridSize = 9;  break;
                case PATTERN3: gridSize = 10; break;
                case PATTERN4: gridSize = 8;  break;
                case PATTERN5: gridSize = 10; break;
                default:       gridSize = 8;  break;
            }

            setPreferredSize(new Dimension(gridSize * CELL_SIZE + 100, gridSize * CELL_SIZE + 100));
            grid = new Cell[gridSize][gridSize];
            for (int r = 0; r < gridSize; r++)
                for (int c = 0; c < gridSize; c++)
                    grid[r][c] = new Cell(r, c);

            if (mode == GameMode.PATTERN1) {
                generateFixedPattern1();
            } else {
                long seed;
                switch (mode) {
                    case PATTERN2: seed = 12399L;  break;
                    case PATTERN3: seed = 44455L;  break;
                    case PATTERN4: seed = 99887L;  break;
                    case PATTERN5: seed = 101010L; break;
                    default:       seed = System.currentTimeMillis(); break;
                }
                generateHighComplexityPuzzle(seed);
            }

            currentHead      = startNode;
            isPlayerTurn     = true;
            gameOver         = false;
            autoSolving      = false;
            computerThinking = false;
            solutionPath     = null;
            lastDpStates     = 0;
            lastDpGoals      = 0;

            TracksGame.this.pack();
            repaint();

            if (sidePanel != null) {
                String nm = (mode == GameMode.RANDOM) ? "Random" : "Pattern " + mode.name().charAt(7);
                sidePanel.updateStatus("Player's Turn");
                sidePanel.updateAnalysis("<html><b>Puzzle:</b> " + nm
                        + "<br/>Size: " + gridSize + "×" + gridSize
                        + "<br/>Track length: " + sumArray(rowClues) + "</html>");
            }
        }

        // ── puzzle generators ────────────────────────────────────

        private void generateFixedPattern1() {
            int[][] coords = {
                {1,1},{1,2},{1,3},{1,4},{2,4},{3,4},
                {3,3},{3,2},{3,1},{4,1},{5,1},
                {5,2},{5,3},{5,4},{5,5},{5,6},
                {4,6},{3,6},{2,6},{1,6}
            };
            List<Cell> path = new ArrayList<>();
            for (int[] p : coords) path.add(grid[p[0]][p[1]]);
            startNode = path.get(0);
            endNode   = path.get(path.size() - 1);
            startNode.hasTrack = true;
            calcClues(path);
        }

        private void generateHighComplexityPuzzle(long seed) {
            Random rand = new Random(seed);
            int target = (int)((gridSize * gridSize) * 0.50);
            while (true) {
                for (int r = 0; r < gridSize; r++)
                    for (int c = 0; c < gridSize; c++)
                        grid[r][c] = new Cell(r, c);
                int sr = rand.nextInt(gridSize - 2) + 1;
                int sc = rand.nextInt(gridSize - 2) + 1;
                boolean[][] vis = new boolean[gridSize][gridSize];
                vis[sr][sc] = true;
                List<Cell> best = new ArrayList<>(), cur = new ArrayList<>();
                cur.add(grid[sr][sc]);
                buildPath(grid[sr][sc], vis, cur, best, target, rand);
                if (best.size() >= target) {
                    startNode = best.get(0);
                    endNode   = best.get(best.size()-1);
                    startNode.hasTrack = true;
                    calcClues(best);
                    break;
                }
            }
        }

        private void buildPath(Cell c, boolean[][] vis, List<Cell> cur,
                                List<Cell> best, int target, Random rand) {
            if (cur.size() > best.size()) { best.clear(); best.addAll(cur); }
            if (best.size() >= target) return;
            List<Cell> nb = getNeighbors(c);
            Collections.shuffle(nb, rand);
            for (Cell n : nb) {
                if (!vis[n.r][n.c] && neighborCount(n, vis) <= 1) {
                    vis[n.r][n.c] = true; cur.add(n);
                    buildPath(n, vis, cur, best, target, rand);
                    if (best.size() >= target) return;
                    cur.remove(cur.size()-1); vis[n.r][n.c] = false;
                }
            }
        }

        private int neighborCount(Cell c, boolean[][] vis) {
            int k = 0;
            for (Cell n : getNeighbors(c)) if (vis[n.r][n.c]) k++;
            return k;
        }

        private void calcClues(List<Cell> path) {
            rowClues = new int[gridSize]; colClues = new int[gridSize];
            for (Cell c : path) { rowClues[c.r]++; colClues[c.c]++; }
        }

        // ── public controls ──────────────────────────────────────

        void resetCurrentPuzzle() { initGame(currentMode); }
        void newRandomPuzzle()    { initGame(GameMode.RANDOM); }
        void loadPattern(int n) {
            GameMode[] m = { GameMode.PATTERN1, GameMode.PATTERN2,
                             GameMode.PATTERN3, GameMode.PATTERN4, GameMode.PATTERN5 };
            if (n >= 1 && n <= 5) initGame(m[n-1]);
        }

        // ── grid helpers ─────────────────────────────────────────

        List<Cell> getNeighbors(Cell c) {
            List<Cell> nb = new ArrayList<>();
            int[] dr = {-1,1,0,0}, dc = {0,0,-1,1};
            for (int i = 0; i < 4; i++) {
                int nr = c.r+dr[i], nc = c.c+dc[i];
                if (nr >= 0 && nr < gridSize && nc >= 0 && nc < gridSize)
                    nb.add(grid[nr][nc]);
            }
            return nb;
        }

        boolean isValidMove(int r, int c) {
            if (r<0||c<0||r>=gridSize||c>=gridSize) return false;
            Cell t = grid[r][c];
            if (t.hasTrack || t.isCrossed) return false;
            if (Math.abs(r-currentHead.r)+Math.abs(c-currentHead.c) != 1) return false;
            if (countTrackInRow(r) >= rowClues[r]) return false;
            if (countTrackInCol(c) >= colClues[c]) return false;
            return true;
        }

        private void makeMove(int r, int c) {
            Cell next = grid[r][c];
            next.hasTrack = true; next.trackParent = currentHead; currentHead = next;
        }

        private int countTrackInRow(int r) {
            int k=0; for (int c=0;c<gridSize;c++) if (grid[r][c].hasTrack) k++; return k;
        }
        private int countTrackInCol(int c) {
            int k=0; for (int r=0;r<gridSize;r++) if (grid[r][c].hasTrack) k++; return k;
        }

        private int idx(int r, int c)          { return r * gridSize + c; }
        private int sumArray(int[] a)           { int s=0; for (int v:a) s+=v; return s; }
        private int manhattan(int id1, int id2) {
            return Math.abs(id1/gridSize - id2/gridSize)
                 + Math.abs(id1%gridSize - id2%gridSize);
        }

        // ════════════════════════════════════════════════════════
        //  PURE DP SOLVER — no backtracking, no A*, no recursion
        // ════════════════════════════════════════════════════════

        /**
         * Solves from the given head cell and returns the list of
         * cells to visit next, or null if no solution exists / timeout.
         */
        private List<Cell> runDPSolver(Cell fromHead) {

            long deadline = System.currentTimeMillis() + timeLimitMs();

            // ── build initial state from already-placed path ──────
            BitSet initVis = new BitSet(gridSize * gridSize);
            int[]  initRow = new int[gridSize];
            int[]  initCol = new int[gridSize];

            for (Cell t = fromHead; t != null; t = t.trackParent) {
                int id = idx(t.r, t.c);
                if (!initVis.get(id)) {
                    initVis.set(id);
                    initRow[t.r]++;
                    initCol[t.c]++;
                }
            }

            int total  = sumArray(rowClues);
            int placed = sumArray(initRow);
            int budget = total - placed;
            if (budget < 0) return null;
            if (budget == 0) return (fromHead == endNode) ? new ArrayList<>() : null;

            int endId = idx(endNode.r, endNode.c);

            // ════════════════════════════════════════════════════
            //  PHASE 1 — Forward BFS
            //  layers[k] : key → DPState  (all states reachable in k steps)
            // ════════════════════════════════════════════════════
            @SuppressWarnings("unchecked")
            HashMap<Long,DPState>[] layers = new HashMap[budget + 1];

            long rootFP  = visFingerprint(initVis);
            long rootKey = stateKey(idx(fromHead.r, fromHead.c), initRow, initCol, rootFP);
            DPState root = new DPState(idx(fromHead.r, fromHead.c),
                                       initRow, initCol, initVis, rootFP, rootKey);

            layers[0] = new HashMap<>();
            layers[0].put(rootKey, root);

            int totalStates = 1;

            for (int depth = 0; depth < budget; depth++) {
                if (System.currentTimeMillis() > deadline) return null;

                HashMap<Long,DPState> curLayer  = layers[depth];
                HashMap<Long,DPState> nextLayer = new HashMap<>();
                layers[depth + 1] = nextLayer;

                int stepsLeft = budget - depth - 1;   // steps remaining after this move

                for (DPState s : curLayer.values()) {
                    Cell curCell = grid[s.cellId / gridSize][s.cellId % gridSize];

                    for (Cell nb : getNeighbors(curCell)) {
                        int nid = idx(nb.r, nb.c);

                        // ── pruning filters ──────────────────────
                        if (grid[nb.r][nb.c].isCrossed)        continue;
                        if (s.visited.get(nid))                 continue;
                        if (s.rowUsed[nb.r] >= rowClues[nb.r]) continue;
                        if (s.colUsed[nb.c] >= colClues[nb.c]) continue;
                        // admissible: must still reach endNode in remaining steps
                        if (stepsLeft < manhattan(nid, endId)) continue;

                        // ── build child arrays ───────────────────
                        int[]  nRow = Arrays.copyOf(s.rowUsed, gridSize);
                        int[]  nCol = Arrays.copyOf(s.colUsed, gridSize);
                        nRow[nb.r]++;
                        nCol[nb.c]++;

                        BitSet nVis = (BitSet) s.visited.clone();
                        nVis.set(nid);
                        long nFP  = updatedFingerprint(s.visitedFP, nid);
                        long nKey = stateKey(nid, nRow, nCol, nFP);

                        // ── de-duplicate: first path wins ────────
                        if (!nextLayer.containsKey(nKey)) {
                            nextLayer.put(nKey, new DPState(nid, nRow, nCol, nVis, nFP, nKey));
                            totalStates++;
                        }
                    }
                }

                if (nextLayer.isEmpty()) return null; // dead end — no surviving states
            }

            // ════════════════════════════════════════════════════
            //  PHASE 2 — Backward DP propagation
            //  Mark goals, then propagate solvability upwards.
            // ════════════════════════════════════════════════════
            int goals = 0;
            for (DPState s : layers[budget].values()) {
                if (s.cellId != endId) continue;
                boolean ok = true;
                for (int i = 0; i < gridSize; i++)
                    if (s.rowUsed[i] != rowClues[i] || s.colUsed[i] != colClues[i]) { ok=false; break; }
                if (ok) { s.solvable = true; goals++; }
            }
            if (goals == 0) return null;

            // Sweep layers from (budget-1) down to 0.
            // A state is solvable iff it has at least one solvable child.
            for (int depth = budget - 1; depth >= 0; depth--) {
                HashMap<Long,DPState> curLayer  = layers[depth];
                HashMap<Long,DPState> nextLayer = layers[depth + 1];
                int stepsLeft = budget - depth - 1;

                for (DPState s : curLayer.values()) {
                    if (s.solvable) continue;   // already marked (won't happen in Phase2 but safe)
                    Cell curCell = grid[s.cellId / gridSize][s.cellId % gridSize];

                    for (Cell nb : getNeighbors(curCell)) {
                        int nid = idx(nb.r, nb.c);
                        if (grid[nb.r][nb.c].isCrossed)        continue;
                        if (s.visited.get(nid))                 continue;
                        if (s.rowUsed[nb.r] >= rowClues[nb.r]) continue;
                        if (s.colUsed[nb.c] >= colClues[nb.c]) continue;
                        if (stepsLeft < manhattan(nid, endId)) continue;

                        int[]  nRow = Arrays.copyOf(s.rowUsed, gridSize);
                        int[]  nCol = Arrays.copyOf(s.colUsed, gridSize);
                        nRow[nb.r]++; nCol[nb.c]++;
                        long nFP  = updatedFingerprint(s.visitedFP, nid);
                        long nKey = stateKey(nid, nRow, nCol, nFP);

                        DPState child = nextLayer.get(nKey);
                        if (child != null && child.solvable) {
                            s.solvable = true;
                            break;   // one solvable child is enough
                        }
                    }
                }
            }

            if (!root.solvable) return null;

            // ════════════════════════════════════════════════════
            //  PHASE 3 — Greedy path extraction
            //  Walk forward; at each layer pick the first solvable child.
            //  No backtracking — the backward pass guarantees at least
            //  one solvable child exists at every step.
            // ════════════════════════════════════════════════════
            List<Cell> path = new ArrayList<>();
            DPState    cur  = root;

            for (int depth = 0; depth < budget; depth++) {
                Cell curCell = grid[cur.cellId / gridSize][cur.cellId % gridSize];
                int stepsLeft = budget - depth - 1;
                HashMap<Long,DPState> nextLayer = layers[depth + 1];

                DPState chosen = null;
                for (Cell nb : getNeighbors(curCell)) {
                    int nid = idx(nb.r, nb.c);
                    if (grid[nb.r][nb.c].isCrossed)        continue;
                    if (cur.visited.get(nid))               continue;
                    if (cur.rowUsed[nb.r] >= rowClues[nb.r]) continue;
                    if (cur.colUsed[nb.c] >= colClues[nb.c]) continue;
                    if (stepsLeft < manhattan(nid, endId)) continue;

                    int[]  nRow = Arrays.copyOf(cur.rowUsed, gridSize);
                    int[]  nCol = Arrays.copyOf(cur.colUsed, gridSize);
                    nRow[nb.r]++; nCol[nb.c]++;
                    long nFP  = updatedFingerprint(cur.visitedFP, nid);
                    long nKey = stateKey(nid, nRow, nCol, nFP);

                    DPState child = nextLayer.get(nKey);
                    if (child != null && child.solvable) { chosen = child; break; }
                }

                if (chosen == null) return null; // should not happen after correct Phase 2
                path.add(grid[chosen.cellId / gridSize][chosen.cellId % gridSize]);
                cur = chosen;
            }

            lastDpStates = totalStates;
            lastDpGoals  = goals;
            return path;
        }

        // ── DP helpers ───────────────────────────────────────────

        /**
         * 64-bit de-duplication key combining cell position,
         * row/col usage counts, and a visited-set fingerprint.
         */
        private long stateKey(int cellId, int[] rowUsed, int[] colUsed, long visFP) {
            long h = (long) cellId * 0x9E3779B97F4A7C15L;
            for (int i = 0; i < gridSize; i++) {
                h ^= Long.rotateLeft((long) rowUsed[i] * 0x6C62272E07BB0142L, i * 5);
                h ^= Long.rotateLeft((long) colUsed[i] * 0xD2B74407B1CE6B31L, i * 7 + 3);
            }
            h ^= visFP * 0xBF58476D1CE4E5B9L;
            return h;
        }

        /** Fingerprint of a full visited BitSet (for the root state). */
        private long visFingerprint(BitSet vis) {
            long fp = 0;
            long[] words = vis.toLongArray();
            for (int i = 0; i < words.length; i++)
                fp ^= Long.rotateLeft(words[i], i * 11);
            return fp;
        }

        /**
         * O(1) incremental update: fingerprint after setting one new bit.
         * Uses the same rotation scheme as visFingerprint so it stays consistent.
         */
        private long updatedFingerprint(long oldFP, int bitId) {
            int word  = bitId / 64;
            int shift = bitId % 64;
            return oldFP ^ Long.rotateLeft(1L << shift, word * 11);
        }

        private long timeLimitMs() {
            if (gridSize <= 8)  return 6000L;
            if (gridSize == 9)  return 10000L;
            return 15000L;
        }

        // ── computer turn ─────────────────────────────────────────

        private void computerTurn() {
            if (gameOver) { computerThinking = false; return; }

            List<Cell> moves = validMovesFrom(currentHead);
            if (moves.isEmpty()) {
                computerThinking = false; gameOver = true;
                sidePanel.updateStatus("Game Over — No moves"); repaint(); return;
            }
            // prefer moves closer to endNode
            moves.sort((a,b) -> Integer.compare(
                    Math.abs(a.r-endNode.r)+Math.abs(a.c-endNode.c),
                    Math.abs(b.r-endNode.r)+Math.abs(b.c-endNode.c)));

            new Thread(() -> {
                Cell best = null;
                for (Cell cand : moves) {
                    // Temporarily register cand as placed to probe the DP
                    Cell op = cand.trackParent;
                    cand.trackParent = currentHead;
                    cand.hasTrack    = true;
                    if (runDPSolver(cand) != null) best = cand;
                    cand.hasTrack    = false;
                    cand.trackParent = op;
                    if (best != null) break;
                }
                final Cell move = (best != null) ? best : moves.get(0);
                SwingUtilities.invokeLater(() -> {
                    computerThinking = false;
                    makeMove(move.r, move.c);
                    sidePanel.updateStatus("Computer → " + move);
                    repaint(); checkGameStatus();
                    if (!gameOver) { isPlayerTurn = true; sidePanel.updateStatus("Your Turn"); }
                });
            }).start();
        }

        private List<Cell> validMovesFrom(Cell head) {
            List<Cell> list = new ArrayList<>();
            for (Cell n : getNeighbors(head)) if (isValidMove(n.r, n.c)) list.add(n);
            return list;
        }

        // ── hint ─────────────────────────────────────────────────

        void showHint() {
            if (gameOver || autoSolving || computerThinking) return;
            clearHints();
            List<Cell> moves = validMovesFrom(currentHead);
            if (moves.isEmpty()) { sidePanel.updateStatus("⚠ No valid moves!"); return; }
            moves.sort((a,b) -> Integer.compare(
                    Math.abs(a.r-endNode.r)+Math.abs(a.c-endNode.c),
                    Math.abs(b.r-endNode.r)+Math.abs(b.c-endNode.c)));

            new Thread(() -> {
                Cell hint = null;
                for (Cell cand : moves) {
                    Cell op = cand.trackParent;
                    cand.trackParent = currentHead; cand.hasTrack = true;
                    if (runDPSolver(cand) != null) hint = cand;
                    cand.hasTrack = false; cand.trackParent = op;
                    if (hint != null) break;
                }
                final Cell fh = hint;
                SwingUtilities.invokeLater(() -> {
                    if (fh != null) {
                        grid[fh.r][fh.c].isHint = true;
                        sidePanel.updateStatus("💡 Hint: " + fh);
                    } else {
                        sidePanel.updateStatus("⚠ No solution found from here");
                    }
                    repaint();
                });
            }).start();
        }

        private void clearHints() {
            for (int r=0; r<gridSize; r++) for (int c=0; c<gridSize; c++) grid[r][c].isHint = false;
        }

        // ── auto-solve ────────────────────────────────────────────

        void autoSolve() {
            if (gameOver || autoSolving || computerThinking) return;
            autoSolving = true; clearHints();
            sidePanel.updateStatus("Solving with Pure DP…");

            new Thread(() -> {
                List<Cell> sol = runDPSolver(currentHead);
                SwingUtilities.invokeLater(() -> {
                    if (sol != null && !sol.isEmpty()) {
                        solutionPath = sol; animateSolution(0);
                    } else {
                        gracefulExit("⚠ Solver timed out",
                                "DP could not finish within the time limit.<br/>Try a smaller puzzle.");
                    }
                });
            }).start();
        }

        private void animateSolution(int i) {
            if (solutionPath == null || i >= solutionPath.size()) {
                autoSolving = false; checkGameStatus(); return;
            }
            Cell next = solutionPath.get(i);
            if (isValidMove(next.r, next.c)) {
                makeMove(next.r, next.c); repaint();
                new Timer().schedule(new TimerTask() {
                    @Override public void run() {
                        SwingUtilities.invokeLater(() -> animateSolution(i + 1));
                    }
                }, 100);
            } else {
                autoSolving = false; gracefulExit("⚠ Stopped", "Stale solution step.");
            }
        }

        private void gracefulExit(String title, String body) {
            autoSolving = false; computerThinking = false;
            isPlayerTurn = false; gameOver = true;
            sidePanel.updateStatus(title);
            sidePanel.updateAnalysis("<html><b>" + body + "</b><br/>Click Reset to try again.</html>");
            repaint();
        }

        // ── game status ───────────────────────────────────────────

        private void checkGameStatus() {
            if (currentHead != endNode) return;
            gameOver = true;
            boolean ok = true;
            for (int i = 0; i < gridSize; i++)
                if (countTrackInRow(i) != rowClues[i] || countTrackInCol(i) != colClues[i])
                    { ok = false; break; }
            if (ok) { sidePanel.updateStatus("🎉 PUZZLE SOLVED!"); showAnalysis(); }
            else      sidePanel.updateStatus("❌ Clues not satisfied");
        }

        private void showAnalysis() {
            int steps = 0; Cell tmp = currentHead;
            while (tmp != null) { steps++; tmp = tmp.trackParent; }
            sidePanel.updateAnalysis(
                "<html><b>✓ Solved!</b><br/>"
                + "<b>Algorithm:</b> Pure DP<br/>"
                + "&nbsp;<b>Phase 1:</b> forward BFS<br/>"
                + "&nbsp;<b>Phase 2:</b> backward propagation<br/>"
                + "&nbsp;<b>Phase 3:</b> greedy extraction<br/>"
                + "<b>DP states:</b> " + lastDpStates + "<br/>"
                + "<b>Goal states:</b> " + lastDpGoals + "<br/>"
                + "<b>Track length:</b> " + steps + "<br/>"
                + "<b>Coverage:</b> " + ((steps * 100) / (gridSize * gridSize)) + "%</html>");
        }

        // ── painting ──────────────────────────────────────────────

        @Override protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            int SX = 50, SY = 50;

            for (int r = 0; r < gridSize; r++) {
                g2.setFont(new Font("SansSerif", Font.BOLD, 14));
                int rc = countTrackInRow(r);
                g2.setColor(rc==rowClues[r] ? new Color(0,150,0) : rc>rowClues[r] ? Color.RED : Color.BLACK);
                g2.drawString(rc + "/" + rowClues[r], SX-45, SY + r*CELL_SIZE + 30);

                for (int c = 0; c < gridSize; c++) {
                    int x = SX + c*CELL_SIZE, y = SY + r*CELL_SIZE;

                    if (r == 0) {
                        int cc = countTrackInCol(c);
                        g2.setColor(cc==colClues[c] ? new Color(0,150,0) : cc>colClues[c] ? Color.RED : Color.BLACK);
                        g2.drawString(cc + "/" + colClues[c], x+12, SY-10);
                    }

                    if (grid[r][c].isHint) {
                        g2.setColor(new Color(255,255,100,200));
                        g2.fillRect(x+1, y+1, CELL_SIZE-2, CELL_SIZE-2);
                    }

                    g2.setColor(Color.LIGHT_GRAY);
                    g2.setStroke(new BasicStroke(1));
                    g2.drawRect(x, y, CELL_SIZE, CELL_SIZE);

                    if (grid[r][c].isCrossed) {
                        g2.setColor(new Color(255,100,100));
                        g2.setStroke(new BasicStroke(3));
                        g2.drawLine(x+10, y+10, x+CELL_SIZE-10, y+CELL_SIZE-10);
                        g2.drawLine(x+CELL_SIZE-10, y+10, x+10, y+CELL_SIZE-10);
                    }

                    if (grid[r][c] == startNode) {
                        g2.setColor(new Color(50,100,255));
                        g2.setFont(new Font("SansSerif", Font.BOLD, 20));
                        g2.drawString("A", x+17, y+33);
                    } else if (grid[r][c] == endNode) {
                        g2.setColor(new Color(220,40,40));
                        g2.setFont(new Font("SansSerif", Font.BOLD, 20));
                        g2.drawString("B", x+17, y+33);
                    }
                }
            }

            // Draw placed track segments
            g2.setColor(new Color(34, 139, 34));
            g2.setStroke(new BasicStroke(6, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
            Cell tmp = currentHead;
            while (tmp != null && tmp.trackParent != null) {
                int x1 = SX + tmp.c*CELL_SIZE + CELL_SIZE/2;
                int y1 = SY + tmp.r*CELL_SIZE + CELL_SIZE/2;
                int x2 = SX + tmp.trackParent.c*CELL_SIZE + CELL_SIZE/2;
                int y2 = SY + tmp.trackParent.r*CELL_SIZE + CELL_SIZE/2;
                g2.drawLine(x1, y1, x2, y2);
                tmp = tmp.trackParent;
            }
        }
    }

    // ════════════════════════════════════════════════════════════════
    //  SIDE PANEL
    // ════════════════════════════════════════════════════════════════

    class SidePanel extends JPanel {
        private JLabel statusLabel, analysisLabel;

        SidePanel(GamePanel gp) {
            setPreferredSize(new Dimension(270, 640));
            setBackground(new Color(240, 240, 245));
            setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
            setBorder(BorderFactory.createEmptyBorder(15, 15, 15, 15));

            JLabel title = new JLabel("Tracks Lab");
            title.setFont(new Font("SansSerif", Font.BOLD, 26));
            title.setForeground(new Color(40, 80, 160));
            title.setAlignmentX(CENTER_ALIGNMENT);

            JLabel algo = new JLabel(
                "<html><center><font color='#555' size='2'>"
                + "<b>Solver: Pure DP (3 phases)</b><br/>"
                + "① forward BFS &nbsp; ② backward propagation<br/>"
                + "③ greedy extraction — zero backtracking"
                + "</font></center></html>");
            algo.setAlignmentX(CENTER_ALIGNMENT);

            JLabel instr = new JLabel(
                "<html><center><b>Controls</b><br/>"
                + "Left-click: place track<br/>"
                + "Right-click: mark ✕<br/><br/>"
                + "<b>Goal:</b> connect A→B &amp; match all clues"
                + "</center></html>");
            instr.setAlignmentX(CENTER_ALIGNMENT);
            instr.setFont(new Font("SansSerif", Font.PLAIN, 12));

            statusLabel = new JLabel("Player's Turn");
            statusLabel.setForeground(new Color(0, 120, 180));
            statusLabel.setFont(new Font("SansSerif", Font.BOLD, 16));
            statusLabel.setAlignmentX(CENTER_ALIGNMENT);

            analysisLabel = new JLabel("");
            analysisLabel.setAlignmentX(CENTER_ALIGNMENT);
            analysisLabel.setVerticalAlignment(SwingConstants.TOP);
            analysisLabel.setFont(new Font("SansSerif", Font.PLAIN, 11));

            JPanel patGrid = new JPanel(new GridLayout(3, 2, 5, 5));
            patGrid.setMaximumSize(new Dimension(230, 120));
            patGrid.setBackground(new Color(240, 240, 245));
            for (int i = 1; i <= 5; i++) {
                final int p = i;
                JButton b = new JButton("Pattern " + i);
                b.setFont(new Font("SansSerif", Font.PLAIN, 11));
                b.addActionListener(e -> gp.loadPattern(p));
                patGrid.add(b);
            }
            JButton rnd = new JButton("Random");
            rnd.setFont(new Font("SansSerif", Font.PLAIN, 11));
            rnd.addActionListener(e -> gp.newRandomPuzzle());
            patGrid.add(rnd);

            JButton hintBtn  = mkBtn("💡 Show Hint",  e -> gp.showHint());
            JButton autoBtn  = mkBtn("🤖 Auto-Solve", e -> gp.autoSolve());
            JButton resetBtn = mkBtn("🔄 Reset",      e -> gp.resetCurrentPuzzle());

            add(title);                             add(gap(6));
            add(algo);                              add(gap(10));
            add(instr);                             add(gap(18));
            add(statusLabel);                       add(gap(18));

            JLabel pl = new JLabel("Select Puzzle:");
            pl.setAlignmentX(CENTER_ALIGNMENT);
            pl.setFont(new Font("SansSerif", Font.BOLD, 14));
            add(pl);                                add(gap(8));
            add(patGrid);                           add(gap(18));
            add(hintBtn);                           add(gap(8));
            add(autoBtn);                           add(gap(8));
            add(resetBtn);                          add(gap(22));
            add(analysisLabel);
        }

        private JButton mkBtn(String t, ActionListener al) {
            JButton b = new JButton(t);
            b.setAlignmentX(CENTER_ALIGNMENT);
            b.setMaximumSize(new Dimension(220, 38));
            b.setFont(new Font("SansSerif", Font.BOLD, 13));
            b.addActionListener(al);
            return b;
        }

        private Component gap(int h) { return Box.createRigidArea(new Dimension(0, h)); }

        void updateStatus(String m)   { statusLabel.setText(m); }
        void updateAnalysis(String h) { analysisLabel.setText(h); }
    }
}
