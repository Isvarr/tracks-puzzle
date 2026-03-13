import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

public class TracksGame extends JFrame {

    private int gridSize = 8;
    private static final int CELL_SIZE = 50;

    private GamePanel gamePanel;
    private SidePanel sidePanel;
    private GameMode currentMode = GameMode.PATTERN1;

    enum GameMode { PATTERN1, PATTERN2, PATTERN3, PATTERN4, PATTERN5, RANDOM }

    public TracksGame() {
        setTitle("Tracks Lab: Backtracking Solver");
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
        try { UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName()); } catch (Exception ignored) {}
        SwingUtilities.invokeLater(() -> new TracksGame().setVisible(true));
    }

    // ══════════════════════════════════════════════════════════════
    //  CELL
    // ══════════════════════════════════════════════════════════════

    static class Cell {
        int r, c;
        boolean hasTrack  = false;
        boolean isCrossed = false;
        boolean isHint    = false;
        Cell    trackParent;

        Cell(int r, int c) { this.r = r; this.c = c; }

        @Override public String toString() { return "(" + r + "," + c + ")"; }
    }

    // ══════════════════════════════════════════════════════════════
    //  GAME PANEL
    // ══════════════════════════════════════════════════════════════

    class GamePanel extends JPanel {

        private Cell[][]   grid;
        private int[]      rowClues, colClues;
        private Cell       startNode, endNode, currentHead;

        private boolean    isPlayerTurn     = true;
        private boolean    gameOver         = false;
        private boolean    autoSolving      = false;
        private boolean    computerThinking = false;
        private List<Cell> solutionPath     = null;

        // Visual search highlight
        private BitSet explored = null;

        // ── constructor ──────────────────────────────────────────

        GamePanel() {
            setBackground(Color.WHITE);
            initGame(GameMode.PATTERN1);

            addMouseListener(new MouseAdapter() {
                @Override public void mouseClicked(MouseEvent e) {
                    if (gameOver || !isPlayerTurn || autoSolving || computerThinking) return;

                    int startX = 50, startY = 50;
                    int c = (e.getX() - startX) / CELL_SIZE;
                    int r = (e.getY() - startY) / CELL_SIZE;

                    if (r >= 0 && c >= 0 && r < gridSize && c < gridSize) {

                        if (SwingUtilities.isRightMouseButton(e)) {
                            if (!grid[r][c].hasTrack && grid[r][c] != currentHead) {
                                grid[r][c].isCrossed = !grid[r][c].isCrossed;
                                repaint();
                            }
                            return;
                        }

                        if (SwingUtilities.isLeftMouseButton(e) && isValidMove(r, c)) {
                            clearHints();
                            makeMove(r, c);
                            repaint();
                            checkGameStatus();

                            if (!gameOver) {
                                computerThinking = true;
                                isPlayerTurn = false;
                                sidePanel.updateStatus("Computer thinking...");
                                new Thread(() -> {
                                    try { Thread.sleep(50); } catch (InterruptedException ignored) {}
                                    computerTurn();
                                }).start();
                            }
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
                case RANDOM:   gridSize = 8;  break;
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
            explored         = new BitSet(gridSize * gridSize);

            TracksGame.this.pack();
            repaint();

            if (sidePanel != null) {
                String modeName = (mode == GameMode.RANDOM) ? "Random"
                        : "Pattern " + mode.name().charAt(7);
                sidePanel.updateStatus("Player's Turn");
                sidePanel.updateAnalysis("<html><b>Puzzle:</b> " + modeName
                        + "<br/>Size: " + gridSize + "x" + gridSize
                        + "<br/>Length: " + sumArray(rowClues) + " cells</html>");
            }
        }

        // ── puzzle generators ────────────────────────────────────

        private void generateFixedPattern1() {
            List<Cell> path = new ArrayList<>();
            int[][] coords = {
                    {1,1},{1,2},{1,3},{1,4},{2,4},{3,4},
                    {3,3},{3,2},{3,1},{4,1},{5,1},
                    {5,2},{5,3},{5,4},{5,5},{5,6},
                    {4,6},{3,6},{2,6},{1,6}
            };
            for (int[] p : coords) path.add(grid[p[0]][p[1]]);
            startNode = path.get(0);
            endNode   = path.get(path.size() - 1);
            startNode.hasTrack = true;
            calculateCluesFromPath(path);
        }

        private void generateHighComplexityPuzzle(long seed) {
            Random rand = new Random(seed);
            int targetLength = (int)((gridSize * gridSize) * 0.50);
            while (true) {
                for (int r = 0; r < gridSize; r++)
                    for (int c = 0; c < gridSize; c++)
                        grid[r][c] = new Cell(r, c);
                int sr = rand.nextInt(gridSize - 2) + 1;
                int sc = rand.nextInt(gridSize - 2) + 1;
                Cell start = grid[sr][sc];
                boolean[][] visited = new boolean[gridSize][gridSize];
                visited[sr][sc] = true;
                List<Cell> bestPath = new ArrayList<>(), currentPath = new ArrayList<>();
                currentPath.add(start);
                generatePathRecursive(start, visited, currentPath, bestPath, targetLength, rand);
                if (bestPath.size() >= targetLength) {
                    startNode = bestPath.get(0);
                    endNode   = bestPath.get(bestPath.size() - 1);
                    startNode.hasTrack = true;
                    calculateCluesFromPath(bestPath);
                    break;
                }
            }
        }

        private void generatePathRecursive(Cell current, boolean[][] visited,
                                            List<Cell> currentPath, List<Cell> bestPath,
                                            int target, Random rand) {
            if (currentPath.size() > bestPath.size()) {
                bestPath.clear(); bestPath.addAll(currentPath);
            }
            if (bestPath.size() >= target) return;
            List<Cell> neighbors = getNeighbors(current);
            Collections.shuffle(neighbors, rand);
            for (Cell n : neighbors) {
                if (!visited[n.r][n.c] && countNeighbors(n, visited) <= 1) {
                    visited[n.r][n.c] = true;
                    currentPath.add(n);
                    generatePathRecursive(n, visited, currentPath, bestPath, target, rand);
                    if (bestPath.size() >= target) return;
                    currentPath.remove(currentPath.size() - 1);
                    visited[n.r][n.c] = false;
                }
            }
        }

        private void calculateCluesFromPath(List<Cell> path) {
            rowClues = new int[gridSize];
            colClues = new int[gridSize];
            for (Cell c : path) { rowClues[c.r]++; colClues[c.c]++; }
        }

        private int countNeighbors(Cell c, boolean[][] visited) {
            int count = 0;
            for (Cell n : getNeighbors(c)) if (visited[n.r][n.c]) count++;
            return count;
        }

        // ── public controls ──────────────────────────────────────

        void resetCurrentPuzzle() { initGame(currentMode); }
        void newRandomPuzzle()    { initGame(GameMode.RANDOM); }
        void loadPattern(int n) {
            switch (n) {
                case 1: initGame(GameMode.PATTERN1); break;
                case 2: initGame(GameMode.PATTERN2); break;
                case 3: initGame(GameMode.PATTERN3); break;
                case 4: initGame(GameMode.PATTERN4); break;
                case 5: initGame(GameMode.PATTERN5); break;
            }
        }

        // ── grid utilities ───────────────────────────────────────

        List<Cell> getNeighbors(Cell c) {
            List<Cell> nb = new ArrayList<>();
            int[] dr = {-1,1,0,0}, dc = {0,0,-1,1};
            for (int i = 0; i < 4; i++) {
                int nr = c.r + dr[i], nc = c.c + dc[i];
                if (nr >= 0 && nr < gridSize && nc >= 0 && nc < gridSize)
                    nb.add(grid[nr][nc]);
            }
            return nb;
        }

        boolean isValidMove(int r, int c) {
            if (r < 0 || c < 0 || r >= gridSize || c >= gridSize) return false;
            Cell target = grid[r][c];
            if (target.hasTrack || target.isCrossed) return false;
            if (Math.abs(r - currentHead.r) + Math.abs(c - currentHead.c) != 1) return false;
            if (countTrackInRow(r) >= rowClues[r]) return false;
            if (countTrackInCol(c) >= colClues[c]) return false;
            return true;
        }

        private void makeMove(int r, int c) {
            Cell next = grid[r][c];
            next.hasTrack    = true;
            next.trackParent = currentHead;
            currentHead      = next;
        }

        private int countTrackInRow(int r) {
            int cnt = 0;
            for (int c = 0; c < gridSize; c++) if (grid[r][c].hasTrack) cnt++;
            return cnt;
        }

        private int countTrackInCol(int c) {
            int cnt = 0;
            for (int r = 0; r < gridSize; r++) if (grid[r][c].hasTrack) cnt++;
            return cnt;
        }

        private int idx(int r, int c) { return r * gridSize + c; }
        private int sumArray(int[] a)  { int s = 0; for (int v : a) s += v; return s; }

        // ── Merge Sort (candidate ordering) ─────────────────────

        private void mergeSort(List<Cell> list, Comparator<Cell> comp) {
            if (list.size() <= 1) return;
            int mid = list.size() / 2;
            List<Cell> L = new ArrayList<>(list.subList(0, mid));
            List<Cell> R = new ArrayList<>(list.subList(mid, list.size()));
            mergeSort(L, comp); mergeSort(R, comp);
            int i = 0, j = 0, k = 0;
            while (i < L.size() && j < R.size())
                list.set(k++, comp.compare(L.get(i), R.get(j)) <= 0 ? L.get(i++) : R.get(j++));
            while (i < L.size()) list.set(k++, L.get(i++));
            while (j < R.size()) list.set(k++, R.get(j++));
        }

        // ══════════════════════════════════════════════════════════
        //  BACKTRACKING DFS SOLVER
        //
        //  Recursively tries every valid neighbour from the current
        //  cell, backtracking when a branch violates row/col clues
        //  or reaches a dead end.  Returns the first complete path
        //  found (start→end with all clues satisfied exactly).
        // ══════════════════════════════════════════════════════════

        private long       btTimeLimit;
        private boolean    btTimedOut;
        private List<Cell> btResult;

        /**
         * Entry point. Returns solution path (cells after fromHead) or null.
         */
        private List<Cell> runBacktracking(Cell fromHead) {
            btTimedOut  = false;
            btResult    = null;
            btTimeLimit = System.currentTimeMillis() + getTimeLimitMs();

            // Reconstruct state from already-placed path
            BitSet visited = new BitSet(gridSize * gridSize);
            int[]  rowUsed = new int[gridSize];
            int[]  colUsed = new int[gridSize];

            Cell t = fromHead;
            while (t != null) {
                int id = idx(t.r, t.c);
                if (!visited.get(id)) {
                    visited.set(id);
                    rowUsed[t.r]++;
                    colUsed[t.c]++;
                }
                t = t.trackParent;
            }

            List<Cell> path = new ArrayList<>();
            btDFS(fromHead, visited, rowUsed, colUsed, path);
            return btTimedOut ? null : btResult;
        }

        private long getTimeLimitMs() {
            return Long.MAX_VALUE; // no timeout
        }

        /**
         * Recursive backtracking DFS.
         *
         * @param cur     current cell (already counted in visited/rowUsed/colUsed)
         * @param visited cells currently on the path
         * @param rowUsed tracks used per row so far
         * @param colUsed tracks used per col so far
         * @param path    cells added after fromHead (for reconstruction)
         * @return true when solution found and btResult is populated
         */
        private boolean btDFS(Cell cur,
                               BitSet visited,
                               int[]  rowUsed,
                               int[]  colUsed,
                               List<Cell> path) {

            if (btTimedOut) return false;

            // ── Goal check ──────────────────────────────────────
            if (cur == endNode) {
                for (int i = 0; i < gridSize; i++) {
                    if (rowUsed[i] != rowClues[i] || colUsed[i] != colClues[i])
                        return false;
                }
                btResult = new ArrayList<>(path);
                return true;
            }

            // ── Pruning: overshoot ───────────────────────────────
            for (int i = 0; i < gridSize; i++) {
                if (rowUsed[i] > rowClues[i] || colUsed[i] > colClues[i]) return false;
            }

            // ── Pruning: capacity (can remaining unvisited cells still satisfy clues?) ──
            for (int i = 0; i < gridSize; i++) {
                int need = rowClues[i] - rowUsed[i];
                if (need > 0) {
                    int avail = 0;
                    for (int j = 0; j < gridSize; j++)
                        if (!visited.get(idx(i, j))) avail++;
                    if (avail < need) return false;
                }
            }
            for (int j = 0; j < gridSize; j++) {
                int need = colClues[j] - colUsed[j];
                if (need > 0) {
                    int avail = 0;
                    for (int i = 0; i < gridSize; i++)
                        if (!visited.get(idx(i, j))) avail++;
                    if (avail < need) return false;
                }
            }

            // ── Visual feedback ──────────────────────────────────
            if (explored != null) explored.set(idx(cur.r, cur.c));

            // ── Expand neighbours (closer to endNode first) ──────
            List<Cell> neighbours = getNeighbors(cur);
            mergeSort(neighbours, (a, b) -> Integer.compare(
                    Math.abs(a.r - endNode.r) + Math.abs(a.c - endNode.c),
                    Math.abs(b.r - endNode.r) + Math.abs(b.c - endNode.c)));

            for (Cell nxt : neighbours) {
                int nid = idx(nxt.r, nxt.c);
                if (grid[nxt.r][nxt.c].isCrossed) continue;
                if (visited.get(nid)) continue;
                if (rowUsed[nxt.r] >= rowClues[nxt.r]) continue;
                if (colUsed[nxt.c] >= colClues[nxt.c]) continue;

                // Place
                visited.set(nid);
                rowUsed[nxt.r]++;
                colUsed[nxt.c]++;
                path.add(nxt);

                if (btDFS(nxt, visited, rowUsed, colUsed, path)) return true;

                // Undo (backtrack)
                path.remove(path.size() - 1);
                visited.clear(nid);
                rowUsed[nxt.r]--;
                colUsed[nxt.c]--;

                if (btTimedOut) return false;
            }

            return false;
        }

        // ── solver entry point ───────────────────────────────────

        private List<Cell> findSolution() {
            List<Cell> sol = runBacktracking(currentHead);
            if (sol != null) {
                System.out.println("[Solver] Solved via Backtracking DFS");
            } else if (btTimedOut) {
                System.out.println("[Solver] Backtracking timed out");
            } else {
                System.out.println("[Solver] No solution exists from this state");
            }
            return sol;
        }

        // ── computer / hint / auto-solve ─────────────────────────

        private void gracefulSolverExit(String title, String reason) {
            autoSolving = false; computerThinking = false;
            isPlayerTurn = false; gameOver = true;
            sidePanel.updateStatus(title);
            sidePanel.updateAnalysis("<html><b>" + reason + "</b><br/><br/>Try another pattern or click <b>Reset</b>.</html>");
            repaint();
        }

        private void computerTurn() {
            if (gameOver) { computerThinking = false; return; }
            List<Cell> validMoves = new ArrayList<>();
            for (Cell n : getNeighbors(currentHead)) if (isValidMove(n.r, n.c)) validMoves.add(n);
            if (validMoves.isEmpty()) {
                computerThinking = false; gameOver = true;
                sidePanel.updateStatus("Game Over - Draw");
                sidePanel.updateAnalysis("<html><font color='orange'><b>No Moves!</b></font></html>");
                repaint(); return;
            }
            mergeSort(validMoves, (a, b) -> Integer.compare(
                    Math.abs(a.r - endNode.r) + Math.abs(a.c - endNode.c),
                    Math.abs(b.r - endNode.r) + Math.abs(b.c - endNode.c)));

            new Thread(() -> {
                Cell best = null;
                for (Cell cand : validMoves) {
                    Cell op = cand.trackParent;
                    cand.trackParent = currentHead;
                    if (runBacktracking(cand) != null) { best = cand; cand.trackParent = op; break; }
                    cand.trackParent = op;
                }
                final Cell move = best;
                SwingUtilities.invokeLater(() -> {
                    computerThinking = false;
                    if (move != null) {
                        makeMove(move.r, move.c);
                        sidePanel.updateStatus("Computer → " + move);
                    } else if (!validMoves.isEmpty()) {
                        makeMove(validMoves.get(0).r, validMoves.get(0).c);
                        sidePanel.updateStatus("Computer → " + validMoves.get(0));
                    } else {
                        gracefulSolverExit("⚠ No Moves", "No valid moves remain."); return;
                    }
                    repaint(); checkGameStatus();
                    if (!gameOver) { isPlayerTurn = true; sidePanel.updateStatus("Your Turn"); }
                });
            }).start();
        }

        void showHint() {
            if (gameOver || autoSolving || computerThinking) return;
            clearHints();
            List<Cell> validMoves = new ArrayList<>();
            for (Cell n : getNeighbors(currentHead)) if (isValidMove(n.r, n.c)) validMoves.add(n);
            if (validMoves.isEmpty()) { sidePanel.updateStatus("⚠ No valid moves!"); return; }
            mergeSort(validMoves, (a, b) -> Integer.compare(
                    Math.abs(a.r - endNode.r) + Math.abs(a.c - endNode.c),
                    Math.abs(b.r - endNode.r) + Math.abs(b.c - endNode.c)));

            new Thread(() -> {
                Cell hint = null;
                for (Cell cand : validMoves) {
                    Cell op = cand.trackParent; cand.trackParent = currentHead;
                    if (runBacktracking(cand) != null) { hint = cand; cand.trackParent = op; break; }
                    cand.trackParent = op;
                }
                final Cell fHint = hint;
                SwingUtilities.invokeLater(() -> {
                    if (fHint != null) { grid[fHint.r][fHint.c].isHint = true; sidePanel.updateStatus("💡 Hint: " + fHint); }
                    else sidePanel.updateStatus("⚠ No solution found!");
                    repaint();
                });
            }).start();
        }

        private void clearHints() {
            for (int r = 0; r < gridSize; r++)
                for (int c = 0; c < gridSize; c++)
                    grid[r][c].isHint = false;
        }

        void autoSolve() {
            if (gameOver || autoSolving || computerThinking) return;
            autoSolving = true; clearHints();
            if (explored != null) explored.clear();
            sidePanel.updateStatus("Auto-solving...");
            new Thread(() -> {
                List<Cell> sol = findSolution();
                SwingUtilities.invokeLater(() -> {
                    if (sol != null && !sol.isEmpty()) { solutionPath = sol; animateSolution(); }
                    else gracefulSolverExit("⚠ Solver Timed Out", "Solver could not finish within the time limit.");
                });
            }).start();
        }

        private void animateSolution() {
            new Thread(() -> {
                for (int i = 0; i < solutionPath.size(); i++) {
                    if (!autoSolving) return;
                    final Cell next = solutionPath.get(i);
                    SwingUtilities.invokeLater(() -> {
                        if (isValidMove(next.r, next.c)) {
                            makeMove(next.r, next.c);
                            repaint();
                        }
                    });
                    try { Thread.sleep(100); } catch (InterruptedException ignored) {}
                }
                SwingUtilities.invokeLater(() -> {
                    autoSolving = false;
                    checkGameStatus();
                });
            }).start();
        }

        private void checkGameStatus() {
            if (currentHead == endNode) {
                gameOver = true;
                boolean valid = true;
                for (int i = 0; i < gridSize; i++)
                    if (countTrackInRow(i) != rowClues[i] || countTrackInCol(i) != colClues[i]) { valid = false; break; }
                if (valid) { sidePanel.updateStatus("🎉 PUZZLE SOLVED! 🎉"); performAnalysis(); }
                else sidePanel.updateStatus("❌ Constraints Failed");
            }
        }

        private void performAnalysis() {
            int steps = 0; Cell tmp = currentHead;
            while (tmp != null) { steps++; tmp = tmp.trackParent; }
            sidePanel.updateAnalysis("<html><b>✓ Solution Verified!</b><br/>"
                    + "<b>Algorithm:</b> Backtracking DFS<br/>"
                    + "<b>Steps:</b> " + steps + "<br/>"
                    + "<b>Coverage:</b> " + ((steps * 100) / (gridSize * gridSize)) + "%</html>");
        }

        // ══════════════════════════════════════════════════════════
        //  PAINT
        // ══════════════════════════════════════════════════════════

        @Override protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            int SX = 50, SY = 50;

            for (int r = 0; r < gridSize; r++) {
                g2.setFont(new Font("SansSerif", Font.BOLD, 14));
                int rCur = countTrackInRow(r);
                g2.setColor(rCur == rowClues[r] ? new Color(0,150,0) : rCur > rowClues[r] ? Color.RED : Color.BLACK);
                g2.drawString(rCur + "/" + rowClues[r], SX - 45, SY + r * CELL_SIZE + 30);

                for (int c = 0; c < gridSize; c++) {
                    int x = SX + c * CELL_SIZE, y = SY + r * CELL_SIZE;
                    if (r == 0) {
                        int cCur = countTrackInCol(c);
                        g2.setColor(cCur == colClues[c] ? new Color(0,150,0) : cCur > colClues[c] ? Color.RED : Color.BLACK);
                        g2.drawString(cCur + "/" + colClues[c], x + 12, SY - 10);
                    }
                    int id = idx(r, c);
                    if (explored != null && explored.get(id)) {
                        g2.setColor(new Color(150, 200, 255, 120));
                        g2.fillRect(x + 1, y + 1, CELL_SIZE - 2, CELL_SIZE - 2);
                    }
                    if (grid[r][c].isHint) {
                        g2.setColor(new Color(255, 255, 100, 200));
                        g2.fillRect(x + 1, y + 1, CELL_SIZE - 2, CELL_SIZE - 2);
                    }
                    g2.setColor(Color.LIGHT_GRAY);
                    g2.setStroke(new BasicStroke(1));
                    g2.drawRect(x, y, CELL_SIZE, CELL_SIZE);
                    if (grid[r][c].isCrossed) {
                        g2.setColor(new Color(255, 100, 100));
                        g2.setStroke(new BasicStroke(3));
                        g2.drawLine(x + 10, y + 10, x + CELL_SIZE - 10, y + CELL_SIZE - 10);
                        g2.drawLine(x + CELL_SIZE - 10, y + 10, x + 10, y + CELL_SIZE - 10);
                    }
                    if (grid[r][c] == startNode) {
                        g2.setColor(new Color(50, 100, 255));
                        g2.setFont(new Font("SansSerif", Font.BOLD, 20));
                        g2.drawString("A", x + 17, y + 33);
                    } else if (grid[r][c] == endNode) {
                        g2.setColor(new Color(255, 50, 50));
                        g2.setFont(new Font("SansSerif", Font.BOLD, 20));
                        g2.drawString("B", x + 17, y + 33);
                    }
                }
            }

            g2.setColor(new Color(34, 139, 34));
            g2.setStroke(new BasicStroke(6, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
            Cell tmp = currentHead;
            while (tmp != null && tmp.trackParent != null) {
                int x1 = SX + tmp.c * CELL_SIZE + CELL_SIZE / 2, y1 = SY + tmp.r * CELL_SIZE + CELL_SIZE / 2;
                int x2 = SX + tmp.trackParent.c * CELL_SIZE + CELL_SIZE / 2, y2 = SY + tmp.trackParent.r * CELL_SIZE + CELL_SIZE / 2;
                g2.drawLine(x1, y1, x2, y2);
                tmp = tmp.trackParent;
            }
        }
    }

    // ══════════════════════════════════════════════════════════════
    //  SIDE PANEL
    // ══════════════════════════════════════════════════════════════

    class SidePanel extends JPanel {
        private JLabel    statusLabel, analysisLabel;
        private GamePanel gamePanel;

        SidePanel(GamePanel gp) {
            this.gamePanel = gp;
            setPreferredSize(new Dimension(270, 600));
            setBackground(new Color(240, 240, 245));
            setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
            setBorder(BorderFactory.createEmptyBorder(15, 15, 15, 15));

            JLabel title = new JLabel("Tracks Lab");
            title.setFont(new Font("SansSerif", Font.BOLD, 26));
            title.setForeground(new Color(40, 80, 160));
            title.setAlignmentX(CENTER_ALIGNMENT);

            JLabel algo = new JLabel("<html><center><font size='2' color='#555'>"
                    + "<b>Solver:</b> Backtracking DFS"
                    + "</font></center></html>");
            algo.setAlignmentX(CENTER_ALIGNMENT);

            JLabel instr = new JLabel("<html><center><b>Controls:</b><br/>"
                    + "Left Click: Add Track<br/>Right Click: Mark ✕<br/><br/>"
                    + "<b>Goal:</b> Connect A→B &amp; match all clues</center></html>");
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

            JPanel patternPanel = new JPanel(new GridLayout(3, 2, 5, 5));
            patternPanel.setMaximumSize(new Dimension(230, 120));
            patternPanel.setBackground(new Color(240, 240, 245));
            for (int i = 1; i <= 5; i++) {
                final int p = i;
                JButton b = new JButton("Pattern " + i);
                b.setFont(new Font("SansSerif", Font.PLAIN, 11));
                b.addActionListener(e -> gp.loadPattern(p));
                patternPanel.add(b);
            }
            JButton rndBtn = new JButton("Random");
            rndBtn.setFont(new Font("SansSerif", Font.PLAIN, 11));
            rndBtn.addActionListener(e -> gp.newRandomPuzzle());
            patternPanel.add(rndBtn);

            JButton hintBtn  = makeBtn("💡 Show Hint",  e -> gp.showHint());
            JButton autoBtn  = makeBtn("🤖 Auto-Solve", e -> gp.autoSolve());
            JButton resetBtn = makeBtn("🔄 Reset",      e -> gp.resetCurrentPuzzle());

            add(title);
            add(Box.createRigidArea(new Dimension(0, 6)));
            add(algo);
            add(Box.createRigidArea(new Dimension(0, 10)));
            add(instr);
            add(Box.createRigidArea(new Dimension(0, 18)));
            add(statusLabel);
            add(Box.createRigidArea(new Dimension(0, 18)));
            JLabel pl = new JLabel("Select Puzzle:");
            pl.setAlignmentX(CENTER_ALIGNMENT);
            pl.setFont(new Font("SansSerif", Font.BOLD, 14));
            add(pl);
            add(Box.createRigidArea(new Dimension(0, 8)));
            add(patternPanel);
            add(Box.createRigidArea(new Dimension(0, 18)));
            add(hintBtn);  add(Box.createRigidArea(new Dimension(0, 10)));
            add(autoBtn);  add(Box.createRigidArea(new Dimension(0, 10)));
            add(resetBtn);
            add(Box.createRigidArea(new Dimension(0, 25)));
            add(analysisLabel);
        }

        private JButton makeBtn(String label, ActionListener al) {
            JButton b = new JButton(label);
            b.setAlignmentX(CENTER_ALIGNMENT);
            b.setMaximumSize(new Dimension(220, 38));
            b.setFont(new Font("SansSerif", Font.BOLD, 13));
            b.addActionListener(al);
            return b;
        }

        void updateStatus(String msg)    { statusLabel.setText(msg); }
        void updateAnalysis(String html) { analysisLabel.setText(html); }
    }
}
